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
      REQUIRE(mAKVGuess.Size()==3,"AKVGuess has Size " << mAKVGuess.Size() 
                                  << ", should be 3.");
    mRad        = p.Get<double>("Radius");
    mSolver     = p.Get<std::string>("Solver","Newton");
    mVerbose    = p.Get<bool>("Verbose",false);
    mPrintResiduals = p.Get<bool>("PrintResiduals",false);
    mL_resid_tol = p.Get<double>("L_resid_tol",1.e-12);
    mv_resid_tol = p.Get<double>("v_resid_tol",1.e-12);
    mMin_thetap = p.Get<double>("Min_thetap",1.e-5);
    mOutput     = p.Get<std::string>("Output");

    printDiagnostic = MyVector<bool>(MV::Size(6), true);
    if(p.OptionIsDefined("DivNorm")) printDiagnostic[0]=p.Get<bool>("DivNorm");
    if(p.OptionIsDefined("VortNorm")) printDiagnostic[1]=p.Get<bool>("VortNorm");
    if(p.OptionIsDefined("SS")) printDiagnostic[2]=p.Get<bool>("SS");
    if(p.OptionIsDefined("XiDivLNorm")) printDiagnostic[3]=p.Get<bool>("XiDivLNorm");
    if(p.OptionIsDefined("fLNorm")) printDiagnostic[4]=p.Get<bool>("fLNorm");
    if(p.OptionIsDefined("fLambdaNorm")) printDiagnostic[5]=p.Get<bool>("fLambdaNorm");
  }

  //==========================================================================
  
  void ComputeAKV::RecomputeData(const DataBoxAccess& box) const {
    delete mResult;

    const StrahlkorperWithMesh& skwm = box.Get<StrahlkorperWithMesh>(mSkwm);
    const SurfaceBasis sb(skwm.Grid());
    DataMesh theta = box.Get<StrahlkorperWithMesh>(mSkwm).Grid().SurfaceCoords()(0);
    DataMesh phi = box.Get<StrahlkorperWithMesh>(mSkwm).Grid().SurfaceCoords()(1);
    DataMesh L(DataMesh::Empty);
    DataMesh v(DataMesh::Empty);

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

    //compute some useful quantities
    const DataMesh rp2 = mRad * Psi * Psi;
    const DataMesh llncf = sb.ScalarLaplacian(log(Psi));
    const DataMesh Ricci = (1.0-2.0*llncf) / (rp2*rp2);
    const Tensor<DataMesh> GradRicci = sb.Gradient(Ricci);

    //creating struct rparams p is not necessary here, but it does save
    //findTHETA and findTtp from having long argument lists and
    //creating the struct on their own
    rparams p = {theta,
                 phi,
                 rp2,
                 sb,
                 llncf,
                 GradRicci,
                 L,
                 v,
                 mL_resid_tol,
                 mv_resid_tol,
                 mPrintResiduals};

    //if the initial guess for thetap is close to zero or pi,
    //try solving at thetap = zero
    bool oneDSolutionFound = false;
    bool thetapGuessIsZero = mAKVGuess[1] < mMin_thetap;
    bool thetapGuessIsPi = M_PI-mAKVGuess[1] < mMin_thetap;
    double THETA = mAKVGuess[0];
    double thetap = mAKVGuess[1];
    double phip = mAKVGuess[2];

    if(thetapGuessIsZero){
      oneDSolutionFound = findTHETA(&p,THETA,mVerbose);
      if(oneDSolutionFound){
        thetap = 0.0;
        phip = 0.0;
      } else { //thetap initial guess is bad and too close to zero
        thetap = mMin_thetap;
      }
    }
    if(thetapGuessIsPi){
      oneDSolutionFound = findTHETA(&p,THETA,mVerbose);
      if(oneDSolutionFound){
        thetap = M_PI;
        phip = 0.0;
        //findTHETA tests for thetap=0.  If thetap=Pi, v and L can be
        //corrected by multiplying by -1
        v *= -1.; L *= -1.;
      } else { //thetap initial guess is bad and too close to Pi
        thetap = M_PI - mMin_thetap;
      }
    }

    //if theta=0 was not a good solution
    //or initial guess was not close to zero
    //try the multidimensional root finder
    if(!oneDSolutionFound){
      findTtp(&p, THETA, thetap, phip, mSolver, mVerbose);

      //if thetap solution is close to zero,
      //and we didn't already try thetap=0,
      //try thetap=0 now
      if( thetap < mMin_thetap && !thetapGuessIsZero){
        const double THETA_saved = THETA;
        const DataMesh v_saved = v;
        const DataMesh L_saved = L;
        oneDSolutionFound = findTHETA(&p,THETA,mVerbose);
        if(oneDSolutionFound){
          thetap = 0.0;
          phip = 0.0;
        } else {
          THETA = THETA_saved;
          v = v_saved; L = L_saved;
        } //end if(oneDSolution)
      } // end if(thetap < mMin_thetap)

      //if thetap solution is close to Pi,
      //and we didn't already try thetap=0,
      //try thetap=0 now
      if( M_PI - thetap < mMin_thetap && !thetapGuessIsPi){
        const double THETA_saved = THETA;
        const DataMesh v_saved = v;
        const DataMesh L_saved = L;
        oneDSolutionFound = findTHETA(&p,THETA,mVerbose);
        if(oneDSolutionFound){
          thetap = M_PI;
          phip = 0.0;
        } else {
          THETA = THETA_saved;
          v = v_saved; L = L_saved;
        } //end if(oneDSolution)
      } // end if(thetap < mMin_thetap)

    } // end if(!oneDSolutionFound)

    //get thetap, phip within standard bounds
    if(thetap < 0.0){
      thetap = -thetap;
      phip -= M_PI;
    }
    if(thetap > M_PI){
      const int m = (thetap/M_PI);
      thetap -= m*M_PI;
      if(m%2) phip -= M_PI;
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


