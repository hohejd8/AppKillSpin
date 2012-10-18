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
    mResidualSize = p.Get<double>("ResidualSize",1.e-11);
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
    const DataMesh theta = box.Get<StrahlkorperWithMesh>(mSkwm).Grid().SurfaceCoords()(0);
    const DataMesh phi = box.Get<StrahlkorperWithMesh>(mSkwm).Grid().SurfaceCoords()(1);
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
    const DataMesh Ricci = 2.0 * (1.0-2.0*llncf) / (rp2*rp2);
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
    double THETA = mAKVGuess[0];
    double thetap = mAKVGuess[1];
    double phip = mAKVGuess[2];
    bool thetapGuessIsZero = thetap < 1.e-5 || M_PI-thetap < 1.e-5;


    //run the appropriate AKV solvers here
    RunAKVsolvers(THETA, thetap, phip, mMin_thetap,
                  mResidualSize, mVerbose, &p, mSolver);

    if(mVerbose){
      std::cout << "Solution found with : Theta  = " << THETA << "\n"
	      << "                      thetap = " << (180.0/M_PI)*thetap << "\n"
	      << "                        phip = " << (180.0/M_PI)*phip 
	      << std::endl;
    }

    //rotate v, Psi for analysis
    DataMesh rotated_v = RotateOnSphere(v,theta,phi,
                                          sb,thetap,phip);
    DataMesh rotated_Psi = RotateOnSphere(Psi,theta,phi,
                                          sb,thetap,phip);

    //determine scale factor at the equator
    const double scaleAtEquator = normalizeKVAtOnePoint(sb, rotated_Psi, rotated_v, mRad, M_PI/2., 0.0);
    if(mVerbose) std::cout << "scale factor at equator = " << scaleAtEquator << std::endl;
    const double avgScale = normalizeKVAtAllPoints(sb, rotated_Psi, theta, phi, rotated_v, mRad);
    MyVector<double> scaleInnerProduct = AKVInnerProduct(v, v, Ricci, rp2, sb);
    if(mVerbose) std::cout << "scale factor ip1 = " << scaleInnerProduct[0] << std::endl;
    if(mVerbose) std::cout << "scale factor ip2 = " << scaleInnerProduct[1] << std::endl;
    if(mVerbose) std::cout << "scale factor ip3 = " << scaleInnerProduct[2] << std::endl;
    std::cout << std::setprecision(15) 
              << scaleInnerProduct[0] << " " 
              << normalizeKVAtAllPoints(sb, rotated_Psi, theta, phi,
                                        rotated_v*scaleAtEquator, mRad) << std::endl;
    std::cout << std::setprecision(15) 
              << scaleInnerProduct[0] << " " 
              << normalizeKVAtAllPoints(sb, rotated_Psi, theta, phi,
                                        rotated_v*scaleInnerProduct[0], mRad) << std::endl;
    std::cout << std::setprecision(15) 
              << scaleInnerProduct[1] << " " 
              << normalizeKVAtAllPoints(sb, rotated_Psi, theta, phi,
                                        rotated_v*scaleInnerProduct[1], mRad) << std::endl;
    std::cout << std::setprecision(15) 
              << scaleInnerProduct[2] << " " 
              << normalizeKVAtAllPoints(sb, rotated_Psi, theta, phi,
                                        rotated_v*scaleInnerProduct[2], mRad) << std::endl;


    TestScaleFactors(rotated_v, rotated_Psi, mRad, sb, theta,
                     phi, scaleAtEquator, scaleInnerProduct[0]);
    //if(mVerbose) std::cout << "scale factor from inner product = " << scaleInnerProduct << std::endl;

    //compute the inner product
/*
    const double residual_ip_equator = AKVInnerProduct(v*scale, v*scale, Ricci, rp2, sb);
    std::cout << "Residual from the inner product of v scaled by the equator = "
              << residual_ip_equator << std::endl;
    const double residual_ip_average = AKVInnerProduct(v*avgScale, v*avgScale, Ricci, rp2, sb);
    std::cout << "Residual from the inner product of v scaled by average scale factor = "
              << residual_ip_average << std::endl;
*/
    //scale L, v
    v *= scaleAtEquator;
    L *= scaleAtEquator;

    //create xi (1-form)
    Tensor<DataMesh> tmp_xi = sb.Gradient(v);
    Tensor<DataMesh> xi(2,"1",DataMesh::Empty);
    xi(0) = tmp_xi(1);
    xi(1) = -tmp_xi(0);

    //perform diagnostics
    KillingDiagnostics(sb, L, Psi, xi, mRad, printDiagnostic);

    //approximate Killing vector
    const DataMesh norm = 1.0 / (Psi*Psi*Psi*Psi*mRad);
    Tensor<DataMesh> xi_vec(3,"1",DataMesh::Empty);
    xi_vec(0) =  norm * ( cos(theta)*cos(phi)*xi(0) - sin(phi)*xi(1) );
    xi_vec(1) =  norm * ( cos(theta)*sin(phi)*xi(0) + cos(phi)*xi(1) );
    xi_vec(2) = -norm * ( sin(theta)*xi(0) );

    mResult = new Tensor<DataMesh>(xi_vec);

  } //end recomputeData

} // namespace ComputeItems


