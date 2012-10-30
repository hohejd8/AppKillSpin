#include <cmath>
#include <iostream>
#include <iomanip>
#include "Utils/StringParsing/ReadFileIntoString.hpp"
#include "Utils/StringParsing/OptionParser.hpp"
//#include "Utils/Tensor/Tensor.hpp"
//#include "Spectral/BasisFunctions/YlmSpherepack.hpp"
#include "Utils/Math/Gaqd.hpp"
#include "Utils/ErrorHandling/Require.hpp"
//#include "Utils/LowLevelUtils/SimpleNorms.hpp"
//#include "Utils/ErrorHandling/Assert.hpp"
#include "AppKillSpin/AKVsolver.hpp"
#include "Spectral/BasisFunctions/SpherePackIterator.hpp"
#include "Utils/LowLevelUtils/Position.hpp"
//This test file essentially copies and replaces
//pieces of the ComputeAKV ComputeItem. The test
//is actually for the AKVsolver functions, which
//are called by ComputeAKV.

int NumberOfTestsFailed = 0;

DataMesh ConstructConformalFactor(const DataMesh& theta,
                                  const DataMesh& phi,
                                  const int& axisym)
{
  //double D = 1.e-9;
  //const double tpd = 2.0+D;
  //const double tmd = 2.0-D;

  //double B = 0.001*(1.0/(tpd*tpd)-1.0/tpd + 0.25);
  //double C = 0.001*(1.0/(tmd*tmd)-1.0/tmd + 0.25);

  DataMesh Psi(theta); //copy constructor to get the right size
      const double tp = M_PI/4.;
      const double pp = M_PI/4.;
  Psi = 1.28 ;
       //+ 0.001*sin(theta)*sin(theta)*cos(2.0*phi)*cos(theta)
       //+ B*(1.0-3.0*cos(theta)*cos(theta)+3.0*sin(theta)*sin(theta)*cos(2.0*phi))
       //+ C*(1.0-3.0*cos(theta)*cos(theta)-3.0*sin(theta)*sin(theta)*cos(2.0*phi));
  switch(axisym){
    case 0: //complete symmetry, constant
      //do nothing to Psi
      std::cout << "COMPLETELY SYMMETRIC" << std::endl;
      break;
    case 1: //z axisymmetry
      Psi += 0.001*cos(theta)*cos(theta);
      std::cout << "Z AXISYMMETRY" << std::endl;
      break;
    case 2: //x axisymmetry
      Psi += 0.001*(1.0-3.0*cos(theta)*cos(theta)+3.0*sin(theta)*sin(theta)*cos(2.0*phi));
      std::cout << "X AXISYMMETRY" << std::endl;
      break;
    case 3: //y axisymmetry
      Psi += 0.001*(1.0-3.0*cos(theta)*cos(theta)-3.0*sin(theta)*sin(theta)*cos(2.0*phi));
      std::cout << "Y AXISYMMETRY" << std::endl;
      break;
    case 4: //no symmetry
      //Psi+=0.001*(  1.0-3.0*cos(theta+tp)*cos(theta+tp)-3.0*sin(theta+tp)*sin(theta+tp)*cos(2.0*phi+pp)
      //            + 1.0-3.0*cos(theta+tp)*cos(theta)+3.0*sin(theta+tp)*sin(theta)*cos(2.0*phi));
      Psi += 0.001*(phi+theta);
      std::cout << "NO AXISYMMETRY" << std::endl;
      break;
    case 5: //off-axis symmetry
      //Psi += 0.001*(3.0*sqrt(2.0)*sin(2.0*theta)*(cos(phi)+sin(phi))
      //               +3.0*sin(theta)*sin(theta)*sin(2.0*phi));

      //Psi += 0.01*(1.0-3.0*cos(theta)*cos(theta)-3.0*sin(theta)*sin(theta)*(cos(phi)+sin(phi)))//;
      //       *sin(2.0*theta);

      Psi += 0.001*(-1.0+3.0*cos(theta)*cos(theta)
             +3.0*sqrt(2.0)*sin(2.0*theta)*(cos(phi)+sin(phi)));

      //Psi+= 0.001*(-1.0+3.0*cos(theta)*cos(theta) //);
      //       +3.0*sqrt(2.0)*sin(2.0*theta)*(cos(phi)+sin(phi))
      //       +3.0*sin(theta)*sin(theta)*sin(2.0*phi));

      std::cout << "OFF-AXIS AXISYMMETRY" << std::endl;
      break;
  } //end switch

  return Psi;
}

void PrintSurfaceNormalization(const DataMesh& v,
                      const DataMesh& rotated_v,
                      const DataMesh& rotated_Psi,
                      const double& rad,
                      const DataMesh& Ricci,
                      const DataMesh& rp2,
                      const SurfaceBasis& sb,
                      const DataMesh& theta,
                      const DataMesh& phi,
                      const double& scaleFactor)
{
  //std::cout << "Scale factor is : " << scaleFactor << std::endl;

  //v *= scaleFactor;
  //rotated_v *= scaleFactor;

      //const double scaleAboveEquator =
      //          normalizeKVAtOnePoint(sb, rotated_Psi, rotated_v, rad, M_PI/4., 0.0);
      //std::cout << "scale factor at theta=Pi/4   : " 
      //          << std::setprecision(12)
      //          << scaleAboveEquator << std::endl;
      //const double scaleBelowEquator =
      //          normalizeKVAtOnePoint(sb, rotated_Psi, rotated_v, rad, 4.*M_PI/5., 0.0);
      //std::cout << "scale factor at theta=4*Pi/5 : " 
      //          << std::setprecision(12)
      //          << scaleBelowEquator << std::endl;
      const double scaleOverSurface =
                normalizeKVAtAllPoints(sb, rotated_Psi, theta, phi, rotated_v*scaleFactor, rad);
      //std::cout << "scale factor over surface    : " 
      //          << std::setprecision(12)
      //          << scaleOverSurface << std::endl;
      //const double scaleInnerProduct = AKVInnerProduct(v, v, Ricci, rp2, sb);
      std::cout << std::setprecision(15) << scaleFactor << " " << scaleOverSurface << std::endl;
  //std::cout << std::endl;
}

int main(){

  const std::string Help =
    "-------------------------------------------------------------------------\n"
    "TestComputeAKV:                                                          \n"
    "-------------------------------------------------------------------------\n"
    "OPTIONS:                                                                 \n"
    "Nth=<int>  theta resolution [default 3]                                  \n"
    "Nph=<int>  phi resolution [default 4]                                    \n"
    "Radius=<double> radius of sphere. [default 1.0]                          \n"
    "AKVGuess=MyVector<double> a guess for the values of THETA, thetap, phip  \n"
    "         [default (0.,0.,0.)]                                            \n"
    "L_resid_tol=<double> tolerance for L residuals when finding approximate  \n"
    "            Killing vectors.  [default 1.e-12]                           \n"
    "v_resid_tol=<double> tolerance for v residuals when finding approximate  \n"
    "            Killing vectors.  [default 1.e-12]                           \n"
    "min_thetap = for values less than this, thetap is considered close to    \n"
    "               zero. [default 1.e-5]                                     \n"
    "ResidualSize=<double> determines the tolerance for residuals from the    \n"
    "             multidimensional root finder.  [default to 1.e-11]          \n"
    "Solver = <std::string> which gsl multidimensional root finding algorith  \n"
    "        should be used. [default Newton]                                 \n"
    "Verbose=<bool> Print spectral coefficients and warnings if true          \n"
    "        [default true]                                                   \n"
    ; 

  std::string Options = ReadFileIntoString("Test.input");
  OptionParser op(Options,Help);
  const int Nth = op.Get<int>("Nth", 3);
  const int Nph = op.Get<int>("Nph", 4);
  const double rad = op.Get<int>("Radius",1.0);
  MyVector<double> AKVGuess =
                op.Get<MyVector<double> >("AKVGuess",MyVector<double>(MV::Size(3),0.0));
      //must be three-dimensional
      REQUIRE(AKVGuess.Size()==3,"AKVGuess has Size " << AKVGuess.Size() 
                                  << ", should be 3.");
  const double L_resid_tol = op.Get<double>("L_resid_tol", 1.e-12);
  const double v_resid_tol = op.Get<double>("L_resid_tol", 1.e-12);
  const double residualSize = op.Get<double>("ResidualSize", 1.e-11);
  const double min_thetap = op.Get<double>("min_theta",1.e-5);
  const std::string solver = op.Get<std::string>("Solver","Newton");
  const bool verbose = op.Get<bool>("Verbose", false);
  const MyVector<bool> printDiagnostic = MyVector<bool>(MV::Size(6), true);

  //create skm
  const StrahlkorperMesh skm(Nth, Nph);
  //create surface basis
  const SurfaceBasis sb(skm);
  //get theta, phi
  const DataMesh theta(skm.SurfaceCoords()(0));
  const DataMesh phi(skm.SurfaceCoords()(1));

  //set the initial guesses to be along particular axes
  const int axes = 3; //the number of perpendicular axes

  //create conformal factors for every rotation
  const int syms = 5; //the number of axisymmetries we are testing

  for(int s=4; s<5; s++){//index over conformal factor symmetries
  //for(int s=0; s<syms; s++){//index over conformal factor symmetries
    //create conformal factor
    const DataMesh Psi = ConstructConformalFactor(theta, phi, s);

    //set the initial guesses to be along particular axes
    double THETA[3] = {0.,0.,0.};
    double thetap[3] = {0.,M_PI/2.,0.};
    double phip[3] = {0.,0.,0.};

    //save the v solutions along particular axes
    MyVector<DataMesh> v(MV::Size(3),DataMesh::Empty);
    MyVector<DataMesh> rotated_v(MV::Size(3),DataMesh::Empty);

    MyVector<Tensor<DataMesh> > xi(MV::Size(axes),Tensor<DataMesh>(2,"1",DataMesh::Empty));

    //compute some useful quantities
    const DataMesh rp2 = rad * Psi * Psi;
    const DataMesh llncf = sb.ScalarLaplacian(log(Psi));
    const DataMesh Ricci = 2.0 * (1.0-2.0*llncf) / (rp2*rp2);
    const Tensor<DataMesh> GradRicci = sb.Gradient(Ricci);

    for(int a=0; a<axes; a++){//index over perpendicular AKV axes
    //for(int a=0; a<1; a++){//index over perpendicular AKV axes
      //for printing
      switch(a){
        case 0:
          std::cout << thetap[0]*180./M_PI << " " << phip[0]*180./M_PI << std::endl;
          std::cout << "z-axis analysis" << std::endl;
          break;
        case 1:
          thetap[1] = M_PI/2.;
          phip[1] = phip[0]+M_PI/2.;
          std::cout << thetap[1]*180./M_PI << " " << phip[1]*180./M_PI << std::endl;
          std::cout << "x-axis analysis" << std::endl;
          break;
        case 2:
          //perform the cross product of the previous two solutions
          const double alpha = sin(thetap[0])*sin(phip[0])*cos(thetap[1])
                            -sin(thetap[1])*sin(phip[1])*cos(thetap[0]);
          const double beta = cos(thetap[0])*sin(thetap[1])*cos(phi[1])
                            -cos(thetap[1])*sin(thetap[0])*cos(phip[0]);
          const double gamma = sin(thetap[0])*cos(phip[0])*sin(thetap[1])*sin(phip[1])
                            -sin(thetap[1])*cos(phip[1])*sin(thetap[0])*sin(phip[0]);
          thetap[2] = atan2(sqrt(alpha*alpha+beta*beta),gamma);
          phip[2] = atan2(beta, gamma);
          std::cout << thetap[2]*180./M_PI << " " << phip[2]*180./M_PI << std::endl;
          std::cout << "y-axis analysis" << std::endl;
          break;
      }
      //create L, v
      DataMesh L(DataMesh::Empty);

      //setup struct with all necessary data
      rparams p = {theta, phi, rp2, sb, llncf, GradRicci,
                   L, v[a], L_resid_tol, v_resid_tol, verbose};

      RunAKVsolvers(THETA[a], thetap[a], phip[a], min_thetap,
                    residualSize, verbose, &p, solver);

      std::cout << "Solution found with : THETA[" << a << "] = " << THETA[a] << "\n"
  	        << "                     thetap[" << a << "] = " << (180.0/M_PI)*thetap[a] << "\n"
	        << "                       phip[" << a << "] = " << (180.0/M_PI)*phip[a] 
	        << std::endl;

      //rotate v, Psi for analysis
      //DataMesh rotated_v = RotateOnSphere(v[a],theta,phi,
      rotated_v[a] = RotateOnSphere(v[a],theta,phi,
                                          sb,thetap[a],phip[a]);

      DataMesh rotated_Psi = RotateOnSphere(Psi,theta,phi,
                                            sb,thetap[a],phip[a]);

      //determine scale factor
      //std::cout << "\nUnscaled results" << std::endl;
      const double scaleAtEquator =
                normalizeKVAtOnePoint(sb, rotated_Psi, rotated_v[a], rad, M_PI/2., 0.0);
      //std::cout << "scale factor at equator      : " 
      //          << std::setprecision(12)
      //          << scaleAtEquator << std::endl;

      //scale L, v
      //v[a] *= scaleAtEquator;
      //rotated_v[a] *= scaleAtEquator;
      //L *= scaleAtEquator;

      //compare scale factors
      //const double scaleAboveEquator =
      //          normalizeKVAtOnePoint(sb, rotated_Psi, rotated_v[a], rad, M_PI/4., 0.0);
      //std::cout << "scale factor at theta=Pi/4   : " 
      //          << std::setprecision(12)
      //          << scaleAboveEquator << std::endl;
      //const double scaleBelowEquator =
      //          normalizeKVAtOnePoint(sb, rotated_Psi, rotated_v[a], rad, 4.*M_PI/5., 0.0);
      //std::cout << "scale factor at theta=4*Pi/5 : " 
      //          << std::setprecision(12)
      //          << scaleBelowEquator << std::endl;
      //const double scaleOverSurface =
      //          normalizeKVAtAllPoints(sb, rotated_Psi, theta, phi, rotated_v[a], rad);
      //std::cout << "scale factor over surface    : " 
      //          << std::setprecision(12)
      //          << scaleOverSurface << std::endl;
      MyVector<double> scaleInnerProduct = AKVInnerProduct(v[a], v[a], Ricci, rp2, sb);//[0];
      //std::cout << "scale factor over surface    : " 
      //          << std::setprecision(12)
      //          << scaleOverSurface << std::endl;

      //std::cout << "\nUsing the scale factor from the path length at the equator" << std::endl;
      PrintSurfaceNormalization(v[a], rotated_v[a], rotated_Psi,
                       rad, Ricci, rp2, sb, theta, phi, scaleAtEquator);
      //std::cout << "Using the scale factor from the inner product / r2p4" << std::endl;
      PrintSurfaceNormalization(v[a], rotated_v[a], rotated_Psi,
                       rad, Ricci, rp2, sb, theta, phi, scaleInnerProduct[0]);
      PrintSurfaceNormalization(v[a], rotated_v[a], rotated_Psi,
                       rad, Ricci, rp2, sb, theta, phi, scaleInnerProduct[1]);
      PrintSurfaceNormalization(v[a], rotated_v[a], rotated_Psi,
                       rad, Ricci, rp2, sb, theta, phi, scaleInnerProduct[2]);
      TestScaleFactors(rotated_v[a], rotated_Psi, rad, sb, theta,
                       phi, scaleAtEquator*0.99, scaleInnerProduct[0]);

      //scale L, v
      v[a] *= scaleAtEquator;
      L *= scaleAtEquator;

      //create xi (1-form)
      Tensor<DataMesh> tmp_xi = sb.Gradient(v[a]);
      Tensor<DataMesh> xi(2,"1",DataMesh::Empty);
      xi(0) = tmp_xi(1);
      xi(1) = -tmp_xi(0);

      //perform diagnostics
      KillingDiagnostics(sb, L, Psi, xi, rad, printDiagnostic);

      std::cout << std::endl;
    }//end loop over perpendicular AKV axes


    //compute inner product for each individual AKV solution
    //const double zz = AKVInnerProduct(v[0], v[0], Ricci, rp2, sb);
    //std::cout << "z-z inner product = " << zz << std::endl;
    //const double xx = AKVInnerProduct(v[1], v[1], Ricci, rp2, sb);
    //std::cout << "x-x inner product = " << xx << std::endl;
    //const double yy = AKVInnerProduct(v[2], v[2], Ricci, rp2, sb);
    //std::cout << "y-y inner product = " << yy << std::endl;


    //compute inner products between AKV solutions
    std::cout << "z-x inner product : " << std::endl;
    //AKVInnerProduct(v[0], v[1], Ricci, rp2, sb);
    std::cout << "z-y inner product : " << std::endl;
    //AKVInnerProduct(v[0], v[2], Ricci, rp2, sb);
    std::cout << "x-y inner product : " << std::endl;
    //AKVInnerProduct(v[1], v[2], Ricci, rp2, sb);
    std::cout << "\n" << std::endl;
  }


  // Return 0 for success
  return NumberOfTestsFailed;
}
