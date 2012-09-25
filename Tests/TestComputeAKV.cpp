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

//This test file essentially copies and replaces
//pieces of the ComputeAKV ComputeItem. The test
//is actually for the AKVsolver functions, which
//are called by ComputeAKV.

int NumberOfTestsFailed = 0;

//enum Axisymmetry {z, x, y, xz, xyz};

DataMesh ConstructConformalFactor(const DataMesh& theta,
                                  const DataMesh& phi,
                                  const int& axisym)
{
  std::cout << "Inside ConstructConformalFactor; axisym = " << axisym << std::endl;
  double D = 1.e-9;
  const double tpd = 2.0+D;
  const double tmd = 2.0-D;

  double B = 0.001*(1.0/(tpd*tpd)-1.0/tpd + 0.25);
  double C = 0.001*(1.0/(tmd*tmd)-1.0/tmd + 0.25);

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
      //Psi += 0.001*cos(theta)*cos(theta);//z axisymmetry
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
  //std::cout << Psi << std::endl;
  return Psi;
}

/*
void ComputeAKV(rparams p, MyVector<double> AKVGuess, const double rad,
                const double min_thetap, const bool verbose)
{

}
*/

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

  double THETA = AKVGuess[0];
  double thetap = AKVGuess[1];
  double phip = AKVGuess[2];

  //create skm
  const StrahlkorperMesh skm(Nth, Nph);
  //create surface basis
  const SurfaceBasis sb(skm);
  //get theta, phi
  const DataMesh theta(skm.SurfaceCoords()(0));
  const DataMesh phi(skm.SurfaceCoords()(1));
  //create L, v
  DataMesh L(DataMesh::Empty);
  DataMesh v(DataMesh::Empty);

//put all of this in a for loop?
/*
  //create conformal factors for every rotation
  const DataMesh& Psi_z = ConstructConformalFactor(theta, phi, z);
  const DataMesh& Psi_x = ConstructConformalFactor(theta, phi, x);
  const DataMesh& Psi_y = ConstructConformalFactor(theta, phi, y);
  const DataMesh& Psi_xz = ConstructConformalFactor(theta, phi, xz);
  const DataMesh& Psi_xyz = ConstructConformalFactor(theta, phi, xyz);
*/

  for(int i=3; i<4; i++){
    std::cout << "iteration = " << i << std::endl;
    DataMesh Psi_ha(sb.CoefficientMesh());
    SpherePackIterator sit(theta.Extents()[0],theta.Extents()[1]);
    Psi_ha[sit(0,0,SpherePackIterator::a)] = 1.;
    Psi_ha[sit(1,1,SpherePackIterator::a)] = 1.;
    //const DataMesh Psi = sb.Evaluate(Psi_ha);
    const DataMesh& Psi = ConstructConformalFactor(theta, phi, i);
    std::cout << "Psi = " << std::endl;
    //std::cout << Psi << std::endl;
    for(int m=0; m<Nth; m++){
      for(int n=0; n<Nph; n++){
        std::cout << Psi[n*Nth+m] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "\n" << std::endl;
    //compute some useful quantities
    const DataMesh rp2 = rad * Psi * Psi;
    const DataMesh llncf = sb.ScalarLaplacian(log(Psi));
    const DataMesh Ricci = (1.0-2.0*llncf) / (rp2*rp2);
    const Tensor<DataMesh> GradRicci = sb.Gradient(Ricci);

    //setup struct with all necessary data
    rparams p = {theta, phi, rp2, sb, llncf, GradRicci,
                 L, v, L_resid_tol, v_resid_tol, verbose};

    RunAKVsolvers(THETA, thetap, phip, min_thetap,
                  residualSize, verbose, &p, solver);

    if(true){
      std::cout << "Solution found with : THETA  = " << THETA << "\n"
	      << "                      thetap = " << (180.0/M_PI)*thetap << "\n"
	      << "                        phip = " << (180.0/M_PI)*phip 
	      << std::endl;
    }

    //determine scale factor
    const double scale = normalizeKVAtOnePoint(sb, Psi, v, rad, M_PI/2., 0.0);
    if(true) std::cout << "scale factor = " << scale << std::endl;

  }



  // Return 0 for success
  return NumberOfTestsFailed;
}
