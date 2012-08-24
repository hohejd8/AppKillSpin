#include "ComputeAKV.hpp"
#include "SurfaceFinder/Strahlkorper/StrahlkorperWithMesh.hpp"
#include "Utils/StringParsing/OptionParser.hpp"
#include "AKVsolver.hpp"
//#include "gsl/gsl_multiroots.h"
#include "Utils/LowLevelUtils/Position.hpp"
#include "SurfaceFinder/StrahlkorperDataSupplier/StrahlkorperParallelInterpolation.hpp"
#include "Dust/Domain/Domain.hpp"
#include "Utils/DataBox/DataBox.hpp"
#include <iomanip>

namespace ComputeItems {

  ComputeAKV::ComputeAKV(const std::string& opts)
    : mResult(0)
  {
std::cout << "begin ComputeAKV constructor " << POSITION << std::endl;
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

    mSolver     = p.Get<std::string>("Solver","Newton");
    mVerbose    = p.Get<bool>("Verbose",false);
    mOutput     = p.Get<std::string>("Output");

    printDiagnostic = MyVector<bool>(MV::Size(6), true);
    if(p.OptionIsDefined("DivNorm")) printDiagnostic[0]=p.Get<bool>("DivNorm");
    if(p.OptionIsDefined("VortNorm")) printDiagnostic[0]=p.Get<bool>("VortNorm");
    if(p.OptionIsDefined("SS")) printDiagnostic[0]=p.Get<bool>("SS");
    if(p.OptionIsDefined("fLNorm")) printDiagnostic[0]=p.Get<bool>("fLNorm");
    if(p.OptionIsDefined("fLambdaNorm")) printDiagnostic[0]=p.Get<bool>("fLambdaNorm");
    if(p.OptionIsDefined("XiDivLNorm")) printDiagnostic[0]=p.Get<bool>("XiDivLNorm");
/*
    if(p.OptionIsDefined("DivNorm")){
      if(p.Get<std::string>("DivNorm")=="true") printDiagnostic[0]=true;
    }
    if(p.OptionIsDefined("VortNorm")){
      if(p.Get<std::string>("VortNorm")=="true") printDiagnostic[1]=true;
    }
    if(p.OptionIsDefined("SS")){
      if(p.Get<std::string>("SS")=="true") printDiagnostic[2]=true;
    }
    if(p.OptionIsDefined("fLNorm")){
      if(p.Get<std::string>("fLNorm")=="true") printDiagnostic[3]=true;
    }
    if(p.OptionIsDefined("fLambdaNorm")){
      if(p.Get<std::string>("fLambdaNorm")=="true") printDiagnostic[4]=true;
    }
    if(p.OptionIsDefined("XiDivLNorm")){
      if(p.Get<std::string>("XiDivLNorm")=="true") printDiagnostic[5]=true;
    }
*/
    mCheat=p.Get<bool>("Cheat",false);
    
std::cout << "end ComputeAKV constructor " << POSITION << std::endl;
  }

  void ComputeAKV::print_state (size_t iter, gsl_multiroot_fsolver * s) const
  {
  printf ("iter = %3u x = % .3f % .3f % .3f "
          "f(x) = % .3e % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_vector_get (s->f, 0), 
          gsl_vector_get (s->f, 1),
          gsl_vector_get (s->f, 2));
  }

  //==========================================================================
  
  void ComputeAKV::RecomputeData(const DataBoxAccess& box) const {
std::cout << "begin ComputeAKV RecomputeData " << POSITION << std::endl;
    delete mResult;

    const gsl_multiroot_fsolver_type *T; //solver type
    gsl_multiroot_fsolver *s; //the actual solver itself
std::cout << "created solver " << POSITION << std::endl;
    int status;
    size_t iter=0;

    const size_t n = 3; //number of dimensions
/*
    int subdomain=0;
    for(int sd=0; sd<box.Size(); ++sd) {
      if(box[sd].KeyExists(mConformalFactor)) subdomain=sd;
    }
*/
    const StrahlkorperWithMesh& skwm = box.Get<StrahlkorperWithMesh>(mSkwm);
std::cout << "created skwm " << POSITION << std::endl;
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
std::cout << "before initialize L,v " << POSITION << std::endl;
    //initialize L, v
    DataMesh theta = box.Get<StrahlkorperWithMesh>(mSkwm).Grid().SurfaceCoords()(0);
    DataMesh L = theta;
    DataMesh v = theta;

    struct rparams p = {skwm,
                        Psi,
                        L,
                        v,
                        1.e-12,
                        1.e-12 };

    gsl_multiroot_function f = {&AKVsolver, n, &p}; //initializes the function

    gsl_vector *x = gsl_vector_alloc(n); //creates initial guess vector

    gsl_vector_set (x, 0, mAKVGuess[0]);
    gsl_vector_set (x, 1, mAKVGuess[1]);
    gsl_vector_set (x, 2, mAKVGuess[2]);
std::cout << "after vector set " << POSITION << std::endl;
    //Declare the appropriate non-derivative root finder
    if(mSolver=="Hybrids") T = gsl_multiroot_fsolver_hybrids;
    else if(mSolver=="Hybrid") T = gsl_multiroot_fsolver_hybrid;
    else if(mSolver=="Newton") T = gsl_multiroot_fsolver_dnewton;
    else if(mSolver=="Broyden") T = gsl_multiroot_fsolver_broyden;
    else std::cout << "Solver option '" << mSolver << "' not valid." << std::endl;

    s = gsl_multiroot_fsolver_alloc(T, n); //original
    gsl_multiroot_fsolver_set(s, &f, x);
      print_state(iter, s);
    //the following line may not be necessary
    status = gsl_multiroot_fsolver_iterate(s); //iterates one time, sets status variable

    do {
      iter++;
//diagnostics
/*
std::cout << "iter = " << iter << std::endl;
      if(status){ //if solver is stuck
        std::cout << "solver is stuck at iter = " << iter << std::endl;
        std::cout << "status is " << status << std::endl;
        break;
      }
*/
//end diagnostics
      status = gsl_multiroot_fsolver_iterate(s);

      print_state(iter, s);

      if(status){ //if solver is stuck
        std::cout << "solver is stuck at iter = " << iter << std::endl;
        break;
      }

      //status = gsl_multiroot_test_residual(s->f, 1e-7); //original
      status = gsl_multiroot_test_residual(s->f, 1e-14); //new

    } while(status == GSL_CONTINUE && iter<1000);

    gsl_multiroot_fsolver_free(s); //frees all memory associated with solver

    double THETA  = gsl_vector_get(s->x,0);
    double thetap = gsl_vector_get(s->x,1);
    double phip   = gsl_vector_get(s->x,2);

//---------------------------------
//testing
if(mCheat){
  std::cout << "I AM CHEATING AND PUTTING THE KNOWN SOLUTION IN HERE" << std::endl;
  THETA = -8.38295e-08*(M_PI/180.);
  thetap = 2.20099e-06*(M_PI/180.);
  phip = 112.398*(M_PI/180.);
}
//end testing
//---------------------------------

//diagnostics
std::cout << "after loop " << POSITION << std::endl;
//end diagnostics

    //get thetap, phip within normal bounds
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

    gsl_vector_free(x);

    if(mVerbose){
      std::cout << "Solution found with : Theta  = " << THETA << "\n"
	      << "                      thetap = " << (180.0/M_PI)*thetap << "\n"
	      << "                        phip = " << (180.0/M_PI)*phip 
	      << std::endl;
    }

    //compute L, v from minimized thetap, phip
    DataMesh phi   = box.Get<StrahlkorperWithMesh>(mSkwm).Grid().SurfaceCoords()(1);
    DataMesh rad   = box.Get<StrahlkorperWithMesh>(mSkwm).Radius();

//for testing only ------------------------------
    //const int mNth = skwm.Grid().SurfaceCoords()(0).Extents()[0];
    //const int mNph = skwm.Grid().SurfaceCoords()(0).Extents()[1];
/*
    std::cout << "v before scaled " << POSITION << std::endl;
    for(int i=0; i<mNth; ++i){
      for(int j=0; j<mNph; ++j) {
        std::cout << std::setprecision(10) << v[i*mNph+j] << " " ;
      }
      std::cout << std::endl;
    }
    std::cout << "\n" << std::endl;
*/
//-----------------------------------------------
    //determine scale factor
    double scale = normalizeKillingVector(&p, thetap, phip);
    if(mVerbose){
      std::cout << "scale factor = " << scale << std::endl;
    }

    //scale L, v
    v *= scale;

//diagnostics
/*
    std::cout << "v scaled " << POSITION << std::endl;
    for(int i=0; i<mNth; ++i){
      for(int j=0; j<mNph; ++j) {
        std::cout << std::setprecision(10) << v[i*mNph+j] << " " ;
      }
      std::cout << std::endl;
    }
    std::cout << "\n" << std::endl;
*/
    //L *= scale;
//end diagnostics

    //create xi (1-form)
    const SurfaceBasis sbe(box.Get<StrahlkorperWithMesh>(mSkwm).Grid());
    Tensor<DataMesh> tmp_xi = sbe.Gradient(v);
    Tensor<DataMesh> xi = tmp_xi;
       //initializes with the right structure, but also copies data.
    xi(0) = tmp_xi(1);
    xi(1) = -tmp_xi(0);

//diagnostics
//----------------------------------------
/*
  std::cout << "xi(0) ComputeAKV" << POSITION << std::endl;
  for(int i=0; i<mNth; ++i){
      for(int j=0; j<mNph; ++j) {
            std::cout << std::setprecision(10) << xi(0)[i*mNph+j] << " " ;
      }
      std::cout << std::endl;
  }
  std::cout << "\n" << std::endl;
  std::cout << "xi(1) ComputeAKV" << POSITION << std::endl;
  for(int i=0; i<mNth; ++i){
      for(int j=0; j<mNph; ++j) {
            std::cout << std::setprecision(10) << xi(1)[i*mNph+j] << " " ;
      }
      std::cout << std::endl;
  }
  std::cout << "\n" << std::endl;
*/
//----------------------------------------
//end diagnostics

    KillingDiagnostics(skwm, L, Psi, xi, printDiagnostic);

    //approximate Killing vector
    const DataMesh norm = 1.0 / (Psi*Psi*Psi*Psi*rad);
    Tensor<DataMesh> xi_vec(3,"1",theta);
    xi_vec(0) =  norm * ( cos(theta)*cos(phi)*xi(0) - sin(phi)*xi(1) );
    xi_vec(1) =  norm * ( cos(theta)*sin(phi)*xi(0) + cos(phi)*xi(1) );
    xi_vec(2) = -norm * ( sin(theta)*xi(0) );



    mResult = new Tensor<DataMesh>(xi_vec);

  } //end recomputeData

} // namespace ComputeItems


