#include <cmath>
#include <iostream>
#include <iomanip>
#include "Utils/StringParsing/ReadFileIntoString.hpp"
#include "Utils/StringParsing/OptionParser.hpp"
#include "Utils/Math/Gaqd.hpp"
#include "Utils/ErrorHandling/Require.hpp"
#include "AppKillSpin/AKVsolver.hpp"
#include "Spectral/BasisFunctions/SpherePackIterator.hpp"
#include "Utils/LowLevelUtils/Position.hpp"
//This test file copies and replaces some
//pieces of the ComputeAKV ComputeItem. The test
//is actually for the AKVsolver functions, which
//are called by ComputeAKV.  Some additional work
//needs to be done to actually make this a test of
//ComputeAKV.

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
      Psi += 0.001*(phi+theta);
      std::cout << "NO AXISYMMETRY" << std::endl;
      break;
    case 5: //off-axis symmetry
      //based on Merzbacher Eq. 16.60
      //Y(L,0)(thetap,phip) = sqrt(4*PI/(2L+1))\Sum(M=-L..L) Y(L,M)(theta,phi) Y(L,M)*(beta,alpha)
      //assume z-symmetry in the new coordinate system (thetap, phip), where the new z-axis had
      //coordinates (beta, alpha) in the old system

      Psi+= 0.01*(5./(32.*M_PI))*(-1.0+3.0*cos(theta)*cos(theta)
		   +3.0*sqrt(2.0)*sin(2.0*theta)*(cos(phi)+sin(phi))
		   +3.0*sin(theta)*sin(theta)*sin(2.0*phi));
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
  const double scaleOverSurface =
                normalizeKVAtAllPoints(sb, rotated_Psi, theta, phi, rotated_v*scaleFactor, rad);

  std::cout << std::setprecision(15) << scaleFactor << " " << scaleOverSurface << std::endl;
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
    "        [default false]                                                  \n"
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

  for(int s=1; s<2; s++){//index over conformal factor symmetries
  //for(int s=5; s<6; s++){//index over conformal factor symmetries
  //for(int s=0; s<syms; s++){//index over conformal factor symmetries
    //create conformal factor
    const DataMesh Psi = ConstructConformalFactor(theta, phi, s);

    //set the initial guesses to be along particular axes
    double THETA[3] = {0.,0.,0.};
    double thetap[3] = {0.,0.,0.}; //guess for z-axisymmetry
    double phip[3] = {0.,0.,0.}; //guess for z-axisymmetry

    //save the v solutions along particular axes
    MyVector<DataMesh> v(MV::Size(3),DataMesh::Empty);
    MyVector<DataMesh> rotated_v(MV::Size(3),DataMesh::Empty);

    MyVector<Tensor<DataMesh> > xi(MV::Size(axes),Tensor<DataMesh>(2,"1",DataMesh::Empty));

    //compute some useful quantities
    const DataMesh rp2 = rad * Psi * Psi;
    const DataMesh r2p4 = rp2*rp2;
    const DataMesh llncf = sb.ScalarLaplacian(log(Psi));
    const DataMesh Ricci = 2.0 * (1.0-2.0*llncf) / (rp2*rp2);
    const Tensor<DataMesh> GradRicci = sb.Gradient(Ricci);

    for(int a=0; a<axes; a++){//index over perpendicular AKV axes
    //for(int a=2; a<3; a++){//index over a particular AKV axis
      //for printing
      switch(a){
        case 0:
          std::cout << "symmetry axis guess : " 
                    << thetap[0]*180./M_PI << " " << phip[0]*180./M_PI << std::endl;
          std::cout << "z-axis analysis" << std::endl;
          break;
        case 1:
          thetap[1] = M_PI/2.;
          phip[1] = atan2(sin(phip[0]), -cos(phip[0]));
          std::cout << "symmetry axis guess : "
                    << thetap[1]*180./M_PI << " " << phip[1]*180./M_PI << std::endl;
          std::cout << "x-axis analysis" << std::endl;
          break;
        case 2:
          //perform the cross product of the previous two solutions
          const double alpha = sin(thetap[0])*sin(phip[0])*cos(thetap[1])
                            -sin(thetap[1])*sin(phip[1])*cos(thetap[0]);
          const double beta = cos(thetap[0])*sin(thetap[1])*cos(phip[1])
                            -cos(thetap[1])*sin(thetap[0])*cos(phip[0]);
          const double gamma = sin(thetap[0])*cos(phip[0])*sin(thetap[1])*sin(phip[1])
                            -sin(thetap[1])*cos(phip[1])*sin(thetap[0])*sin(phip[0]);
          thetap[2] = atan2(sqrt(alpha*alpha+beta*beta),gamma);
          phip[2] = atan2(beta, alpha);
          std::cout << "symmetry axis guess : "
                    << thetap[2]*180./M_PI << " " << phip[2]*180./M_PI << std::endl;
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
      //create xi (1-form)
      Tensor<DataMesh> tmp_xi = sb.Gradient(v[a]);
      xi[a](0) = tmp_xi(1);
      xi[a](1) = -tmp_xi(0);

      //perform diagnostics
      //note that Psi and xi are unscaled and unrotated at this point
      KillingDiagnostics(sb, L, Psi, xi[a], rad, printDiagnostic);

      //compare scale factors
      std::cout << "\n" << POSITION << std::endl;
      const double scaleAtEquator =
                normalizeKVAtOnePoint(sb, rotated_Psi, rotated_v[a], rad, M_PI/2., 0.0);
                //normalizeKVAtOnePoint(sb, rotated_Psi, rotated_v[a], rad, M_PI/4., 0.0);
      //std::cout<<"scale factor at equator: "<<std::setprecision(12)<<scaleAtEquator<<std::endl;
      MyVector<double> scaleInnerProduct = InnerProductScaleFactors(v[a], v[a], Ricci, r2p4, sb);
      PrintSurfaceNormalization(v[a], rotated_v[a], rotated_Psi,
                       rad, Ricci, rp2, sb, theta, phi, scaleAtEquator);
      PrintSurfaceNormalization(v[a], rotated_v[a], rotated_Psi,
                       rad, Ricci, rp2, sb, theta, phi, scaleInnerProduct[0]);
      PrintSurfaceNormalization(v[a], rotated_v[a], rotated_Psi,
                       rad, Ricci, rp2, sb, theta, phi, scaleInnerProduct[1]);
      PrintSurfaceNormalization(v[a], rotated_v[a], rotated_Psi,
                       rad, Ricci, rp2, sb, theta, phi, scaleInnerProduct[2]);
      //OptimizeScaleFactor(rotated_v[a], rotated_Psi, rad, sb, theta,
      //   phi, scaleAtEquator, scaleInnerProduct[0], scaleInnerProduct[1], scaleInnerProduct[2]);

      //scale L, v
      v[a] *= scaleAtEquator;
      L *= scaleAtEquator;

      //create xi (1-form)
      //tmp_xi = sb.Gradient(v[a]);
      //xi[a](0) = tmp_xi(1);
      //xi[a](1) = -tmp_xi(0);

      std::cout << std::endl;
    }//end loop over perpendicular AKV axes


    //compute inner product for each individual AKV solution
    const double zz = AKVInnerProduct(xi[0], xi[0], Ricci, r2p4, sb);
    std::cout << "z-z inner product = " << zz << std::endl;
    const double xx = AKVInnerProduct(xi[1], xi[1], Ricci, r2p4, sb);
    std::cout << "x-x inner product = " << xx << std::endl;
    const double yy = AKVInnerProduct(xi[2], xi[2], Ricci, r2p4, sb);
    std::cout << "y-y inner product = " << yy << std::endl;


    //compute inner products between AKV solutions
    std::cout << "z-x inner product : "
              << AKVInnerProduct(xi[0], xi[1], Ricci, r2p4, sb) << std::endl;
    std::cout << "x-y inner product : "
              << AKVInnerProduct(xi[1], xi[2], Ricci, r2p4, sb) << std::endl;
    std::cout << "y-z inner product : "
              << AKVInnerProduct(xi[2], xi[0], Ricci, r2p4, sb) << std::endl;
    std::cout << "\n" << std::endl;
  }


  // Return 0 for success
  return NumberOfTestsFailed;
}
