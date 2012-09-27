#include <cmath>
#include <iostream>
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
      //std::cout << Psi << std::endl;
      break;
    case 1: 
      Psi += 0.001*cos(theta)*cos(theta);//z axisymmetry
      //Psi += 0.1*sin(theta)*sin(theta);//fake axisymmetry
      break;
    case 2: //x axisymmetry
      Psi += 0.001*(1.0-3.0*cos(theta)*cos(theta)+3.0*sin(theta)*sin(theta)*cos(2.0*phi));
      break;
    case 3: //y axisymmetry
      Psi += 0.001*(1.0-3.0*cos(theta)*cos(theta)-3.0*sin(theta)*sin(theta)*cos(2.0*phi));
      break;
    case 4: //xz axisymmetry
      Psi +=  0.001*(-1.0+3.0*cos(theta)*cos(theta)+3.0*sin(theta)*sin(theta)*cos(2.0*phi)
           +6.0*sin(2.0*theta)*cos(phi));
      break;
    case 5: //xyz axisymmetry
      Psi +=  0.001*(-1.0+3.0*cos(theta)*cos(theta)
           +3.0*sqrt(2.0)*sin(2.0*theta)*(cos(phi)+sin(phi))
           +3.0*sin(theta)*sin(theta)*sin(2.0*phi));
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

    MyVector<Tensor<DataMesh> > xi(MV::Size(axes),Tensor<DataMesh>(2,"1",DataMesh::Empty));

    //compute some useful quantities
    const DataMesh rp2 = rad * Psi * Psi;
    const DataMesh llncf = sb.ScalarLaplacian(log(Psi));
    const DataMesh Ricci = (1.0-2.0*llncf) / (rp2*rp2);
    const Tensor<DataMesh> GradRicci = sb.Gradient(Ricci);

    for(int a=0; a<axes; a++){//index over perpendicular AKV axes
      //create L, v
      DataMesh L(DataMesh::Empty);
      DataMesh v(DataMesh::Empty);

      //setup struct with all necessary data
      rparams p = {theta, phi, rp2, sb, llncf, GradRicci,
                   L, v, L_resid_tol, v_resid_tol, verbose};

      RunAKVsolvers(THETA[a], thetap[a], phip[a], min_thetap,
                    residualSize, verbose, &p, solver);

      if(true){
        std::cout << "Solution found with : THETA[" << a << "] = " << THETA[a] << "\n"
  	          << "                     thetap[" << a << "] = " << (180.0/M_PI)*thetap[a] << "\n"
	          << "                       phip[" << a << "] = " << (180.0/M_PI)*phip[a] 
	          << std::endl;
      }

      //rotate v, Psi for analysis
      DataMesh rotated_v = RotateOnSphere(v,theta,phi,
                                          sb,thetap[a],phip[a]);

      DataMesh rotated_Psi = RotateOnSphere(Psi,theta,phi,
                                            sb,thetap[a],phip[a]);

      //determine scale factor
      const double scale = normalizeKVAtOnePoint(sb, rotated_Psi, rotated_v, rad, M_PI/2., 0.0);
      //if(true) std::cout << "scale factor = " << scale << std::endl;

      //const double residual_ip_equator = AKVInnerProduct(v*scale, v*scale, Ricci, /*rp2,*/ sb);
      //std::cout << "Residual from the inner product of v scaled by the equator = "
      //          << residual_ip_equator << std::endl;

      //scale L, v
      v *= scale;
      L *= scale;

      //create xi (1-form)
      Tensor<DataMesh> tmp_xi = sb.Gradient(v);
      xi[a](0) = tmp_xi(1);
      xi[a](1) = -tmp_xi(0);

      //KillingDiagnostics(sb, L, Psi[s], xi[a], rad, MyVector<bool>(MV::Size(6),true) );
    }//end loop over perpendicular AKV axes

    //compute inner products between AKV solutions
    const double zx = AKVInnerProduct(xi[0], THETA[0], xi[1], THETA[1], Ricci, sb);
    std::cout << "z-x inner product = " << zx << std::endl;
    const double zy = AKVInnerProduct(xi[0], THETA[0], xi[2], THETA[2], Ricci, sb);
    std::cout << "z-y inner product = " << zy << std::endl;
    const double xy = AKVInnerProduct(xi[1], THETA[1], xi[2], THETA[2], Ricci, sb);
    std::cout << "x-y inner product = " << xy << std::endl;

    std::cout << "\n" << std::endl;
  }


  // Return 0 for success
  return NumberOfTestsFailed;
}
