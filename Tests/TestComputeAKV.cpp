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

  Psi = 1.28 ;
       //+ 0.001*sin(theta)*sin(theta)*cos(2.0*phi)*cos(theta)
       //+ B*(1.0-3.0*cos(theta)*cos(theta)+3.0*sin(theta)*sin(theta)*cos(2.0*phi))
       //+ C*(1.0-3.0*cos(theta)*cos(theta)-3.0*sin(theta)*sin(theta)*cos(2.0*phi));
  switch(axisym){
    case 0: //constant
      //do nothing to Psi
      std::cout << "COMPLETELY SYMMETRIC" << std::endl;
      break;
    case 1: 
      Psi += 0.001*cos(theta)*cos(theta);//z axisymmetry
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
    case 4: //xz axisymmetry
      Psi +=  0.001*(-1.0+3.0*cos(theta)*cos(theta)+3.0*sin(theta)*sin(theta)*cos(2.0*phi)
           +6.0*sin(2.0*theta)*cos(phi));
      std::cout << "xz axisymmetry" << std::endl;
      break;
    case 5: //xyz axisymmetry
      Psi +=  0.001*(-1.0+3.0*cos(theta)*cos(theta)
           +3.0*sqrt(2.0)*sin(2.0*theta)*(cos(phi)+sin(phi))
           +3.0*sin(theta)*sin(theta)*sin(2.0*phi));
      std::cout << "xyz axisymmetry" << std::endl;
      break;
  } //end switch

  return Psi;
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
  const int syms = 4; //the number of axisymmetries we are testing


  for(int s=0; s<syms; s++){//index over conformal factor symmetries
    //create conformal factor
    const DataMesh Psi = ConstructConformalFactor(theta, phi, s);

    //set the initial guesses to be along particular axes
    double THETA[3] = {0.,0.,0.};
    double thetap[3] = {0.,M_PI/2.,M_PI/2.};
    double phip[3] = {0.,0.,M_PI/2.};

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
      //for printing
      switch(a){
        case 0:
          std::cout << "z-axis analysis" << std::endl;
          break;
        case 1:
          std::cout << "x-axis analysis" << std::endl;
          break;
        case 2:
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
      const double scaleAtEquator =
                normalizeKVAtOnePoint(sb, rotated_Psi, rotated_v[a], rad, M_PI/2., 0.0);
      std::cout << "scale factor at equator      : " 

                << scaleAtEquator << std::endl;

      //scale L, v
      v[a] *= scaleAtEquator;
      rotated_v[a] *= scaleAtEquator;
      L *= scaleAtEquator;

      //compare scale factors
      const double scaleAboveEquator =
                normalizeKVAtOnePoint(sb, rotated_Psi, rotated_v[a], rad, M_PI/4., 0.0);
      std::cout << "scale factor at theta=Pi/4   : " 
                << std::setprecision(12)
                << scaleAboveEquator << std::endl;
      const double scaleBelowEquator =
                normalizeKVAtOnePoint(sb, rotated_Psi, rotated_v[a], rad, 4.*M_PI/5., 0.0);
      std::cout << "scale factor at theta=4*Pi/5 : " 
                << std::setprecision(12)
                << scaleBelowEquator << std::endl;
      const double scaleOverSurface =
                normalizeKVAtAllPoints(sb, rotated_Psi, theta, phi, rotated_v[a], rad);
      std::cout << "scale factor over surface    : " 
                << std::setprecision(12)
                << scaleOverSurface << std::endl;
      //scale L, v
      //v[a] *= scale;
      //rotated_v[a] *= scaleAtEquator;
      //L *= scale;

      //create xi (1-form)
      Tensor<DataMesh> tmp_xi = sb.Gradient(v[a]);
      xi[a](0) = tmp_xi(1);
      xi[a](1) = -tmp_xi(0);

      //KillingDiagnostics(sb, L, Psi[s], xi[a], rad, MyVector<bool>(MV::Size(6),true) );

      //print out diagnostics


      AKVInnerProduct(v[a], v[a], Ricci, rp2, sb);
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
    AKVInnerProduct(v[0], v[1], Ricci, rp2, sb);
    std::cout << "z-y inner product : " << std::endl;
    AKVInnerProduct(v[0], v[2], Ricci, rp2, sb);
    std::cout << "x-y inner product : " << std::endl;
    AKVInnerProduct(v[1], v[2], Ricci, rp2, sb);
    std::cout << "\n" << std::endl;
  }


  // Return 0 for success
  return NumberOfTestsFailed;
}
