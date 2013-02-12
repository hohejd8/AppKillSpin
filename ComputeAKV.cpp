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
    mAKVSolution     = p.Get<std::string>("AKVSolution");
    mWithRicciScaling = p.Get<bool>("WithRicciScaling",true);
    mSkwm       = p.Get<std::string>("StrahlkorperWithMesh");
    mConformalFactor = p.Get<std::string>("ConformalFactor","ConformalFactor");
    mInterpolateConformalFactor = p.Get<bool>("InterpolateConformalFactor",true);
    mAKVGuess   = p.Get<MyVector<double> >("AKVGuess");//must be three-dimensional
      REQUIRE(mAKVGuess.Size()==3,"AKVGuess has Size " << mAKVGuess.Size() 
                                  << ", should be 3.");
    mRad        = p.Get<double>("Radius");
    mSolver     = p.Get<std::string>("Solver","Newton");
    mVerbose    = p.Get<bool>("Verbose",false);
    mPrintResiduals = p.Get<bool>("PrintResiduals",false);
    mResidualSize = p.Get<double>("ResidualSize",1.e-10);
    mL_resid_tol = p.Get<double>("L_resid_tol",1.e-12);
    mv_resid_tol = p.Get<double>("v_resid_tol",1.e-12);
    mMin_thetap = p.Get<double>("Min_thetap",1.e-5);
    mPrintTtpSolution = p.Get<bool>("PrintTtpSolution",true);
    mPrintInnerProducts = p.Get<bool>("PrintInnerProducts",false);
    mScaleFactor = p.Get<std::string>("ScaleFactor","Equator");
      REQUIRE(   mScaleFactor=="Equator"
              || mScaleFactor=="InnerProduct1"
              || mScaleFactor=="InnerProduct2"
              || mScaleFactor=="InnerProduct3"
              || mScaleFactor=="InnerProduct4"
              || mScaleFactor=="InnerProduct5"
              || mScaleFactor=="InnerProduct6"
              || mScaleFactor=="InnerProducts"
              || mScaleFactor=="Optimize",
              "ScaleFactor is '" << mScaleFactor << "'; must be Equator,"
              " InnerProduct1, InnerProduct2, InnerProduct3, or Optimize.")
    mPrintScaleFactor = p.Get<bool>("PrintScaleFactor",false);
    mPrintSurfaceNormalization = p.Get<bool>("PrintSurfaceNormalization",false);
    mPrintSteps = p.Get<bool>("PrintSteps",false);
    mPrintBisectionResults = p.Get<bool>("PrintBisectionResults",false);
    printDiagnostic = MyVector<bool>(MV::Size(6), true);
    if(p.OptionIsDefined("DivNorm")) printDiagnostic[0]=p.Get<bool>("DivNorm");
    if(p.OptionIsDefined("VortNorm")) printDiagnostic[1]=p.Get<bool>("VortNorm");
    if(p.OptionIsDefined("SS")) printDiagnostic[2]=p.Get<bool>("SS");
    if(p.OptionIsDefined("XiDivLNorm")) printDiagnostic[3]=p.Get<bool>("XiDivLNorm");
    if(p.OptionIsDefined("fLNorm")) printDiagnostic[4]=p.Get<bool>("fLNorm");
    if(p.OptionIsDefined("fLambdaNorm")) printDiagnostic[5]=p.Get<bool>("fLambdaNorm");
  }

  //==========================================================================
  
  void ComputeAKV::RecomputeData(const DataBoxAccess& boxa) const {
    delete mResult;
    const StrahlkorperWithMesh& skwm = boxa.Get<StrahlkorperWithMesh>(mSkwm);
    const SurfaceBasis sb(skwm.Grid().Basis());
    const DataMesh theta = boxa.Get<StrahlkorperWithMesh>(mSkwm).Grid().SurfaceCoords()(0);
    const DataMesh phi = boxa.Get<StrahlkorperWithMesh>(mSkwm).Grid().SurfaceCoords()(1);


      DataBox localBox("ComputeAKV DataBox");
      DataBoxInserter localBoxInserter(localBox, POSITION);
      DataBoxAccess lba(localBox, "ComputeAKV");
    //if Psi needs to be interpolated onto the surface
    if(mInterpolateConformalFactor){
      const DataBoxAccess& rootDBA = boxa.Root();
      const Domain& D=rootDBA.Get<Domain>("Domain");
      MyVector<std::string> TensorsToInterp(MV::fill, mConformalFactor);
      std::string mSpatialCoordMap="";
      StrahlkorperDataSuppliers::StrahlkorperParallelInterpolation
          supplier(D, rootDBA, D.Communicator(),
                   TensorsToInterp, mSpatialCoordMap,"Spectral");
      //DataBox localBox("ComputeAKV DataBox");
      //DataBoxInserter localBoxInserter(localBox, POSITION);
      for(int i=0; i<TensorsToInterp.Size(); i++){
        localBoxInserter.AddVolatileItem(TensorsToInterp[i],
                             supplier.Supply(skwm, TensorsToInterp[i]));
      }
      //DataBoxAccess lba(localBox, "ComputeAKV");
      //const DataMesh& Psi(lba.Get<Tensor<DataMesh> >(mConformalFactor)());
    } else {
      localBoxInserter.AddVolatileItem(mConformalFactor,
                  Tensor<DataMesh>(1,"1",boxa.Get<DataMesh>(mConformalFactor)));
    }
    const DataMesh& Psi(lba.Get<Tensor<DataMesh> >(mConformalFactor)());

    //compute some useful quantities
    const DataMesh rp2 = mRad * Psi * Psi;
    const DataMesh r2p4 = rp2*rp2;
    const DataMesh llncf = sb.ScalarLaplacian(log(Psi));
    const DataMesh Ricci = 2.0 * (1.0-2.0*llncf) / r2p4;
    const Tensor<DataMesh> GradRicci = sb.Gradient(Ricci);

    const int axes = 3;

    //if the initial guess for thetap is close to zero or pi,
    //try solving at thetap = zero
    //set the initial guesses
    double THETA[3] = {mAKVGuess[0],0.0,0.0};
    double thetap[3] = {mAKVGuess[1],0.,0.};
    double phip[3] = {mAKVGuess[2],0.,0.};

    //save the v, xi solutions along particular axes
    MyVector<DataMesh> v(MV::Size(3),DataMesh::Empty);
    MyVector<DataMesh> rotated_v(MV::Size(3),DataMesh::Empty);
    MyVector<Tensor<DataMesh> > xi(MV::Size(axes),Tensor<DataMesh>(2,"1",DataMesh::Empty));

    //save the <v_i|v_j> inner product solutions
    double v0v0 = 0.;
    double v1v1 = 0.;
    double v2v2 = 0.;
    double v0v1 = 0.;
    double v0v2 = 0.;
    double v1v2 = 0.;

    //tolerances; add as input later
    const double symmetry_tol = 1.e-11;
    const double min_thetap = 1.e-5;

    //if the diagnostics below decide that there is a bad solution for v[a]
    //(usually a repeated solution), this flag will indicate that the
    //solver should be run again
    bool badAKVSolution = false;

    for(int a=0; a<axes; a++){//index over orthonormal directions to find AKV solutions

      //generate a guess for the next axis of symmetry based on prior solutions.
      //do not run AxisInitialGuess if the previous solution was bad; a new guess
      //has already been provided
      if(!badAKVSolution) AxisInitialGuess(thetap, phip, a);

      badAKVSolution = false;

      //create L
      DataMesh L(DataMesh::Empty);

      //setup struct with all necessary data
      rparams p = {theta,phi,rp2,sb,llncf,GradRicci,L,
                 v[a],mL_resid_tol,mv_resid_tol,mPrintResiduals,mWithRicciScaling};

      RunAKVsolvers(THETA[a], thetap[a], phip[a], min_thetap,
                    mResidualSize, mVerbose, &p, mSolver);

      if(mPrintTtpSolution){
        std::cout << "Solution found with : THETA[" << a << "] = " << THETA[a] << "\n"
  	<< "                     thetap[" << a << "] = " << (180.0/M_PI)*thetap[a] << "\n"
	<< "                       phip[" << a << "] = " << (180.0/M_PI)*phip[a] 
	<< std::endl;
      }

      //check inner products
      // <v_i|v_j> = Integral 0.5 * (Ricci) * Grad(v_i) \cdot Grad(v_j) dA
      switch(a){
        case 0:
          //compute inner product <v_0|v_0>
          v0v0 = AKVInnerProduct(v[0],v[0],Ricci,sb,mWithRicciScaling)*sqrt(2.)*M_PI;
          break;
        case 1:
          //compute inner products <v_1|v_1>, <v_0|v_1>
          v1v1 = AKVInnerProduct(v[1],v[1],Ricci,sb,mWithRicciScaling)*sqrt(2.)*M_PI;
          v0v1 = AKVInnerProduct(v[0],v[1],Ricci,sb,mWithRicciScaling)*sqrt(2.)*M_PI;
          if(fabs(v0v0) == fabs(v0v1)) badAKVSolution = true;
          //if(fabs(v0v1) > 1.e-04) badAKVSolution = true;
          break;
        case 2:
          //compute inner products <v_2|v_2>, <v_0|v_2>, <v_1|v_2>
          v2v2 = AKVInnerProduct(v[2],v[2],Ricci,sb,mWithRicciScaling)*sqrt(2.)*M_PI;
          v0v2 = AKVInnerProduct(v[0],v[2],Ricci,sb,mWithRicciScaling)*sqrt(2.)*M_PI;
          v1v2 = AKVInnerProduct(v[1],v[2],Ricci,sb,mWithRicciScaling)*sqrt(2.)*M_PI;
          if(fabs(v0v0) == fabs(v0v2)) badAKVSolution = true;
          if(fabs(v1v1) == fabs(v1v2)) badAKVSolution = true;
          //if(fabs(v0v2) > 1.e-04 || fabs(v1v2) > 1.e-04) badAKVSolution = true;
          break;
      }

      if(mPrintInnerProducts && a==2){
          std::cout << "<v_0|v_0> = " << v0v0 << std::endl;
          std::cout << "-THETA_0 <v_0|v_0> = " << -THETA[0]*v0v0 << std::endl;
          std::cout << "<v_1|v_1> = " << v1v1 << std::endl;
          std::cout << "<v_0|v_1> = " << v0v1 << std::endl;
          std::cout << "-THETA_1 <v_1|v_1> = " << -THETA[1]*v1v1 << std::endl;
          std::cout << "<v_2|v_2> = " << v2v2 << std::endl;
          std::cout << "<v_0|v_2> = " << v0v2 << std::endl;
          std::cout << "<v_1|v_2> = " << v1v2 << std::endl;
          std::cout << "-THETA_2 <v_2|v_2> = " << -THETA[2]*v2v2 << std::endl;
      }

      //Gram Schmidt orthogonalization
      switch(a){
        case 1:
          if(v0v0<symmetry_tol && v1v1<symmetry_tol){ //two symmetries, v2v2 should also be symmetric
            GramSchmidtOrthogonalization(v[0], v0v0, v[1], v0v1);
          }
          break;
        case 2:
          if(v0v0<symmetry_tol){
            if(v1v1<symmetry_tol){
              REQUIRE(v2v2<symmetry_tol, "Three symmetries required, but only two found.");
              GramSchmidtOrthogonalization(v[0], v0v0, v[2], v0v2);
              GramSchmidtOrthogonalization(v[1], v1v1, v[2], v1v2);
            } else if(v2v2<symmetry_tol){
              REQUIRE(false, "Three symmetries required, but only two found.");
            } else {
              GramSchmidtOrthogonalization(v[1], v1v1, v[2], v1v2);
            }
          } else if(v1v1<symmetry_tol){
            if(v2v2<symmetry_tol){
              REQUIRE(false, "Three symmetries required, but only two found.");
            } else {
              GramSchmidtOrthogonalization(v[0], v0v0, v[2], v0v2);
            }
          } else if(v2v2<symmetry_tol){
              GramSchmidtOrthogonalization(v[0], v0v0, v[1], v0v1);
          }
        break;
      }

      //create xi (1-form)
      xi[a] = ComputeXi(v[a], sb);

      //perform diagnostics
      //Psi and xi are unscaled and unrotated
      KillingDiagnostics(sb, L, Psi, xi[a], mRad, printDiagnostic);

      //rotate v, Psi for analysis
      rotated_v[a] = RotateOnSphere(v[a],theta,phi,
                                          sb,thetap[a],phip[a]);

      DataMesh rotated_Psi = RotateOnSphere(Psi,theta,phi,
                                            sb,thetap[a],phip[a]);

      //determine scale factors
      double scale = 0.0;
      if(mPrintSurfaceNormalization)
        std::cout << "Scale factor:       RMS deviation:" << std::endl;
      if(mScaleFactor=="Equator"){
        scale = NormalizeAKVAtOnePoint(sb, rotated_Psi, rotated_v[a], mRad, M_PI/2., 0.0, mPrintSteps);
      } else if(mScaleFactor=="InnerProduct1"){
        scale = InnerProductScaleFactors(v[a], v[a], Ricci, r2p4, sb,true)[0];
      } else if(mScaleFactor=="InnerProduct2"){
        scale = InnerProductScaleFactors(v[a], v[a], Ricci, r2p4, sb,true)[1];
      } else if(mScaleFactor=="InnerProduct3"){
        scale = InnerProductScaleFactors(v[a], v[a], Ricci, r2p4, sb,true)[2];
      } else if(mScaleFactor=="InnerProduct4"){
        scale = InnerProductScaleFactors(v[a], v[a], Ricci, r2p4, sb,false)[0];
      } else if(mScaleFactor=="InnerProduct5"){
        scale = InnerProductScaleFactors(v[a], v[a], Ricci, r2p4, sb,false)[1];
      } else if(mScaleFactor=="InnerProduct6"){
        scale = InnerProductScaleFactors(v[a], v[a], Ricci, r2p4, sb,false)[2];
      } else if(mScaleFactor=="InnerProducts"){
        std::cout << POSITION << " this is a test bed that can be deleted later" << std::endl;
        const MyVector<double> scaleIP 
              = InnerProductScaleFactors(v[a], v[a], Ricci, r2p4, sb,true);
        const double artificialScale = 0.9;
        std::cout << POSITION << " artificial scaling of rotated_v by " 
                  << artificialScale << std::endl;
        const MyVector<double> artificialScaleIP 
              = InnerProductScaleFactors(v[a]*artificialScale, v[a]*artificialScale,
                                         Ricci, r2p4, sb,true);
        scale = scaleIP[0];
          std::cout << "IP1 " << std::endl;
          PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,rotated_v[a],scaleIP[0],mRad,mPrintSteps);
          std::cout << "IP1 artificial scale" << std::endl;
          PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,rotated_v[a],
                                    artificialScaleIP[0],mRad,mPrintSteps);
          std::cout << "IP2 " << std::endl;
          PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,rotated_v[a],scaleIP[1],mRad,mPrintSteps);
          std::cout << "IP3 " << std::endl;
          PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,rotated_v[a],scaleIP[2],mRad,mPrintSteps);
      } else if(mScaleFactor=="Optimize"){
        const double scaleAtEquator
              = NormalizeAKVAtOnePoint(sb, rotated_Psi, rotated_v[a], mRad, M_PI/2., 0.0, mPrintSteps);
        const MyVector<double> scaleIP 
              = InnerProductScaleFactors(v[a], v[a], Ricci, r2p4, sb,true);
std::cout << POSITION << " Introduce IP4, IP5, IP6" << std::endl;
        const MyVector<double> scaleIPNoRicci
              = InnerProductScaleFactors(v[a], v[a], Ricci, r2p4, sb,false);
        scale = OptimizeScaleFactor(rotated_v[a], rotated_Psi, mRad, sb, theta,
         phi, scaleAtEquator, scaleIP[0], scaleIP[1], scaleIP[2],mPrintSteps,mPrintBisectionResults);
        if(mPrintSurfaceNormalization){
          std::cout << "Equator " << std::endl;
          PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,
                                    rotated_v[a],scaleAtEquator,mRad,mPrintSteps);
          std::cout << "IP1 " << std::endl;
          PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,rotated_v[a],scaleIP[0],mRad,mPrintSteps);
          std::cout << "IP2 " << std::endl;
          PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,rotated_v[a],scaleIP[1],mRad,mPrintSteps);
          std::cout << "IP3 " << std::endl;
          PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,rotated_v[a],scaleIP[2],mRad,mPrintSteps);
          std::cout << "IP4 " << std::endl;
          PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,rotated_v[a],
                                    scaleIPNoRicci[0],mRad,mPrintSteps);
          std::cout << "IP5 " << std::endl;
          PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,rotated_v[a],
                                    scaleIPNoRicci[1],mRad,mPrintSteps);
          std::cout << "IP6 " << std::endl;
          PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,rotated_v[a],
                                    scaleIPNoRicci[2],mRad,mPrintSteps);
        }
      }

      if(mPrintSurfaceNormalization){
        std::cout << mScaleFactor << std::endl;
        PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,rotated_v[a],scale,mRad,mPrintSteps);
      }

      //scale v
      v[a] *= scale;

      //recompute scaled xi (1-form)
      xi[a] = ComputeXi(v[a], sb);

      if(badAKVSolution){
        v[a] = 0.;
        thetap[a] += M_PI/4.; 
        phip[a] += M_PI/4.;
        a--;
        std::cout << "This was a bad / repeated solution, and will be recomputed." << std::endl;
      }
      std::cout << std::endl;
    }//end loop over perpendicular AKV axes

    std::cout << "\n" << std::endl;
  

    //approximate Killing vector
    const DataMesh norm = 1.0 / (Psi*Psi*Psi*Psi*mRad);
    Tensor<DataMesh> xi_vec(3,"1",DataMesh::Empty);
    xi_vec(0) =  norm * ( cos(theta)*cos(phi)*xi[0](0) - sin(phi)*xi[0](1) );
    xi_vec(1) =  norm * ( cos(theta)*sin(phi)*xi[0](0) + cos(phi)*xi[0](1) );
    xi_vec(2) = -norm * ( sin(theta)*xi[0](0) );

    mResult = new Tensor<DataMesh>(xi_vec);

  } //end recomputeData

} // namespace ComputeItems


