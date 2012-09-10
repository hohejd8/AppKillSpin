#include "ComputeAKV.hpp"
#include "SurfaceFinder/Strahlkorper/StrahlkorperWithMesh.hpp"
#include "Utils/StringParsing/OptionParser.hpp"
#include "AKVsolver.hpp"
#include "Utils/LowLevelUtils/Position.hpp"
#include "SurfaceFinder/StrahlkorperDataSupplier/StrahlkorperParallelInterpolation.hpp"
#include "Dust/Domain/Domain.hpp"
#include "Utils/DataBox/DataBox.hpp"
#include <iomanip>

namespace ComputeItems {

  ComputeAKV::ComputeAKV(const std::string& opts)
    : mResult(0)
  {
    OptionParser p(opts, Help());
    mSkwm       = p.Get<std::string>("StrahlkorperWithMesh");
    mConformalFactor = p.Get<std::string>("ConformalFactor","ConformalFactor");
    mAKVGuess   = p.Get<MyVector<double> >("AKVGuess");//must be three-dimensional
    {
      bool thetap_moved = false;
      REQUIRE(mAKVGuess.Size()==3,"AKVGuess has Size " << mAKVGuess.Size() 
                                  << ", should be 3.");
      if(mAKVGuess[1]<1.e-5){
        mAKVGuess[1] = 1.e-5;
        thetap_moved = true;
      }
      if(mAKVGuess[1]>(M_PI-1.e-5)){
        mAKVGuess[1] = M_PI - 1.e-5;
        thetap_moved = true;
      }
      if(thetap_moved){
        std::cout << "thetap is too close to a pole; the initial guess has been \n"
                     "moved off-axis by an amount 1.e-5." << std::endl;
      }
    }

    mRad        = p.Get<double>("Radius");
    mSolver     = p.Get<std::string>("Solver","Newton");
    mVerbose    = p.Get<bool>("Verbose",false);
    mPrintResiduals = p.Get<bool>("PrintResiduals",false);
    mOutput     = p.Get<std::string>("Output");

    printDiagnostic = MyVector<bool>(MV::Size(6), true);
    if(p.OptionIsDefined("DivNorm")) printDiagnostic[0]=p.Get<bool>("DivNorm");
    if(p.OptionIsDefined("VortNorm")) printDiagnostic[1]=p.Get<bool>("VortNorm");
    if(p.OptionIsDefined("SS")) printDiagnostic[2]=p.Get<bool>("SS");
    if(p.OptionIsDefined("fLNorm")) printDiagnostic[3]=p.Get<bool>("fLNorm");
    if(p.OptionIsDefined("fLambdaNorm")) printDiagnostic[4]=p.Get<bool>("fLambdaNorm");
    if(p.OptionIsDefined("XiDivLNorm")) printDiagnostic[5]=p.Get<bool>("XiDivLNorm");
  }


  void ComputeAKV::print_state (size_t iter, gsl_multiroot_fsolver * s) const
  {
    std::cout << "iter = " << iter
              << " x = " << std::setprecision(8) << std::setw(10) << gsl_vector_get(s->x,0)
              << " " << std::setprecision(8) << std::setw(10) << gsl_vector_get(s->x,1)
              << " " << std::setprecision(8) << std::setw(10) << gsl_vector_get(s->x,2)
              << " f(x) = " << std::setprecision(8) << std::setw(10) << gsl_vector_get(s->f,0)
              << " " << std::setprecision(8) << std::setw(10) << gsl_vector_get(s->f,1)
              << " " << std::setprecision(8) << std::setw(10) << gsl_vector_get(s->f,2)
              << std::endl << std::flush;
  }

  //==========================================================================
  
  void ComputeAKV::RecomputeData(const DataBoxAccess& box) const {
    delete mResult;

    const gsl_multiroot_fsolver_type *T; //solver type
    gsl_multiroot_fsolver *s; //the actual solver itself
    int status;
    size_t iter=0;

    const size_t n = 3; //number of dimensions

    const StrahlkorperWithMesh& skwm = box.Get<StrahlkorperWithMesh>(mSkwm);
    const SurfaceBasis sb(skwm.Grid());

    //Psi needs to be interpolated onto the surface
    const Domain& D=box.Get<Domain>("Domain");
    MyVector<std::string> TensorsToInterp(MV::fill, mConformalFactor);
    std::string mSpatialCoordMap="";
    StrahlkorperDataSuppliers::StrahlkorperParallelInterpolation
          supplier(D, box, box.Get<Domain>("Domain").Communicator(),
                   TensorsToInterp, mSpatialCoordMap,"Spectral");
    DataBox localBox("ComputeAKV DataBox");
    DataBoxInserter localBoxInserter(localBox, POSITION);
    for(int i=0; i<TensorsToInterp.Size(); i++){
      localBoxInserter.AddVolatileItem(TensorsToInterp[i],
                           supplier.Supply(skwm, TensorsToInterp[i]));
    }
    DataBoxAccess lba(localBox, "AKV Recompute");
    const DataMesh& Psi(lba.Get<Tensor<DataMesh> >(mConformalFactor)());

    //initialize L, v
    DataMesh theta = box.Get<StrahlkorperWithMesh>(mSkwm).Grid().SurfaceCoords()(0);
    DataMesh phi = box.Get<StrahlkorperWithMesh>(mSkwm).Grid().SurfaceCoords()(1);
    DataMesh L(DataMesh::Empty);
    DataMesh v(DataMesh::Empty);

    rparams p = {theta,
                 phi,
                 mRad,
                 sb,
                 Psi,
                 L,
                 v,
                 1.e-12,
                 1.e-12,
                 mPrintResiduals};

    gsl_multiroot_function f = {&AKVsolver, n, &p}; //initializes the function

    gsl_vector *x = gsl_vector_alloc(n); //creates initial guess vector
    gsl_vector_set (x, 0, mAKVGuess[0]);
    gsl_vector_set (x, 1, mAKVGuess[1]);
    gsl_vector_set (x, 2, mAKVGuess[2]);

    //Declare the appropriate non-derivative root finder
    if(mSolver=="Hybrids") T = gsl_multiroot_fsolver_hybrids;
    else if(mSolver=="Hybrid") T = gsl_multiroot_fsolver_hybrid;
    else if(mSolver=="Newton") T = gsl_multiroot_fsolver_dnewton;
    else if(mSolver=="Broyden") T = gsl_multiroot_fsolver_broyden;
    else std::cout << "Solver option '" << mSolver << "' not valid." << std::endl;

    s = gsl_multiroot_fsolver_alloc(T, n);
    gsl_multiroot_fsolver_set(s, &f, x);
    if(mVerbose) print_state(iter, s);

    do {
      iter++;

      status = gsl_multiroot_fsolver_iterate(s);

      if(mVerbose) print_state(iter, s);

      if(status){ //if solver is stuck
        std::cout << "GSL multiroot solver is stuck at iter = " << iter << std::endl;
        break;
      }

      status = gsl_multiroot_test_residual(s->f, 1e-11);

    } while(status == GSL_CONTINUE && iter<1000);

    if(iter==1000){
      std::cout << "Iteration was stopped at iter=1000. \n"
                   "You may want to check the solution for validity."
                << std::endl;
    }

    double THETA  = gsl_vector_get(s->x,0);
    double thetap = gsl_vector_get(s->x,1);
    double phip   = gsl_vector_get(s->x,2);

    gsl_vector_free(x); //frees all memory associated with vector x
    gsl_multiroot_fsolver_free(s); //frees all memory associated with solver

    //get thetap, phip within standard bounds
    if(thetap < 0.0){
      thetap = -thetap;
      phip -= M_PI;
    }
    if(thetap > M_PI){
      const int m = (thetap/M_PI);
      thetap -= m*M_PI;
      if(n%2) phip -= M_PI;
    }
    if(phip <= -M_PI){
      const int m = (M_PI - phip)/(2.0*M_PI);
      phip += 2.0*m*M_PI;
    } else if(phip > M_PI) {
      const int m = (phip + M_PI)/(2.0*M_PI);
      phip -= 2.0*m*M_PI;
    }

    if(mVerbose){
      std::cout << "Solution found with : Theta  = " << THETA << "\n"
	      << "                      thetap = " << (180.0/M_PI)*thetap << "\n"
	      << "                        phip = " << (180.0/M_PI)*phip 
	      << std::endl;
    }

    //compute L, v from minimized thetap, phip
    //determine scale factor
    double scale = normalizeKillingVector(sb, Psi, v, mRad);
    if(mVerbose){
      std::cout << "scale factor = " << scale << std::endl;
    }

    //scale L, v
    v *= scale;
    L *= scale;

    //create xi (1-form)
    const SurfaceBasis sbe(box.Get<StrahlkorperWithMesh>(mSkwm).Grid());
    Tensor<DataMesh> tmp_xi = sbe.Gradient(v);
    Tensor<DataMesh> xi(2,"1",DataMesh::Empty);
    xi(0) = tmp_xi(1);
    xi(1) = -tmp_xi(0);

    KillingDiagnostics(sb, L, Psi, xi, mRad, printDiagnostic);

    //approximate Killing vector
    const DataMesh norm = 1.0 / (Psi*Psi*Psi*Psi*mRad);
    Tensor<DataMesh> xi_vec(3,"1",theta);
    xi_vec(0) =  norm * ( cos(theta)*cos(phi)*xi(0) - sin(phi)*xi(1) );
    xi_vec(1) =  norm * ( cos(theta)*sin(phi)*xi(0) + cos(phi)*xi(1) );
    xi_vec(2) = -norm * ( sin(theta)*xi(0) );

    mResult = new Tensor<DataMesh>(xi_vec);

  } //end recomputeData

} // namespace ComputeItems


