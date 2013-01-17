#include <cmath>
#include <iostream>
#include <iomanip>
#include "Utils/StringParsing/ReadFileIntoString.hpp"
#include "Utils/StringParsing/OptionParser.hpp"
#include "Utils/Math/Gaqd.hpp"
#include "Utils/ErrorHandling/Require.hpp"
#include "AppKillSpin/AKVsolver.hpp"
#include "Spectral/BasisFunctions/SpherePackIterator.hpp"//can probably be deleted
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
  DataMesh Psi(theta); //copy constructor to get the right size

  Psi = 1.28 ;

  switch(axisym){
    case 0: //complete symmetry, constant
      //do nothing to Psi
      //Psi = 1.0; //delete this later
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
      //Psi += 0.001*(phi+theta); //too difficult to analize with small l,m
      //Psi += 0.001*(sin(theta)*(cos(phi)+sin(theta)) + cos(theta)*sin(2.*phi));
      Psi += 0.001*(0.//-1.0+3.0*cos(theta)*cos(theta)
		   +3.0*sqrt(2.0)*sin(2.0*theta)*(cos(phi)+sin(phi))
		   +3.0*sin(theta)*sin(theta)*sin(2.0*phi));
      std::cout << "NO AXISYMMETRY" << std::endl;
      break;
    case 5: //off-axis symmetry
      //based on Merzbacher Eq. 16.60
      //Y(L,0)(thetap,phip) = sqrt(4*PI/(2L+1))\Sum(M=-L..L) Y(L,M)(theta,phi) Y(L,M)*(beta,alpha)
      //assume z-symmetry in the new coordinate system (thetap, phip), where the new z-axis had
      //coordinates (beta, alpha) in the old system

      Psi+= 0.001//*(5./(32.*M_PI))
                   *(-1.0+3.0*cos(theta)*cos(theta)
		   +3.0*sqrt(2.0)*sin(2.0*theta)*(cos(phi)+sin(phi))
		   +3.0*sin(theta)*sin(theta)*sin(2.0*phi));
      std::cout << "OFF-AXIS AXISYMMETRY" << std::endl;
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
    "symmetry_tol=<double> abs(THETA) must be less than this value to be      \n"
    "             considered an exact symmetry.  [default 1.e-11]             \n"
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
  const double rad = op.Get<double>("Radius",1.0);
  MyVector<double> AKVGuess =
                op.Get<MyVector<double> >("AKVGuess",MyVector<double>(MV::Size(3),0.0));
      //must be three-dimensional
      REQUIRE(AKVGuess.Size()==3,"AKVGuess has Size " << AKVGuess.Size() 
                                  << ", should be 3.");
  const double L_resid_tol = op.Get<double>("L_resid_tol", 1.e-12);
  const double v_resid_tol = op.Get<double>("L_resid_tol", 1.e-12);
  const double residualSize = op.Get<double>("ResidualSize", 1.e-11);
  const double min_thetap = op.Get<double>("min_theta",1.e-5);
  const double symmetry_tol = op.Get<double>("symmetry_tol",1.e-11);
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

    //set the initial guesses
    double THETA[3] = {AKVGuess[0],0.,0.};
    double thetap[3] = {AKVGuess[1],0.,0.};
    double phip[3] = {AKVGuess[2],0.,0.};

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
    //int symmetries[3] = 0; //counts the number of symmetries

    //compute some useful quantities
    const DataMesh rp2 = rad * Psi * Psi;
    const DataMesh r2p4 = rp2*rp2;
    const DataMesh llncf = sb.ScalarLaplacian(log(Psi));
    const DataMesh Ricci = 2.0 * (1.0-2.0*llncf) / r2p4;
    const Tensor<DataMesh> GradRicci = sb.Gradient(Ricci);

    for(int a=0; a<axes; a++){//index over perpendicular axes to find AKV solutions

      //if the diagnostics below decide that there is a bad solution for v[a]
      //(usually a repeated solution), this flag will indicate that the
      //solver should be run again
      bool badAKVSolution = false;

      //generate a guess for the next axis of symmetry based on prior solutions.
      AxisInitialGuess(thetap, phip, a);

      //create L
      DataMesh L(DataMesh::Empty);

      //setup struct with all necessary data
      rparams p = {theta, phi, rp2, sb, llncf, GradRicci,
                   L, v[a], L_resid_tol, v_resid_tol, verbose, true};

      RunAKVsolvers(THETA[a], thetap[a], phip[a], min_thetap,
                    residualSize, verbose, &p, solver);

      std::cout << "Solution found with : THETA[" << a << "] = " << THETA[a] << "\n"
  	        << "                     thetap[" << a << "] = " << (180.0/M_PI)*thetap[a] << "\n"
	        << "                       phip[" << a << "] = " << (180.0/M_PI)*phip[a] 
	        << std::endl;

      //check inner products
      // <v_i|v_j> = Integral 0.5 * Ricci * Grad(v_i) \cdot Grad(v_j) dA
      switch(a){
        case 0:
          //compute inner product <v_0|v_0>
          v0v0 = AKVInnerProduct(v[0],v[0],Ricci,sb)*sqrt(2.)*M_PI;
          //if(v0v0<symmetry_tol) //symmetries++;
          std::cout << "<v_0|v_0> = " << v0v0 << std::endl;
          std::cout << "-THETA <v_0|v_0> = " << -THETA[a]*v0v0 << std::endl;
          break;
        case 1:
          //compute inner products <v_1|v_1>, <v_0|v_1>
          v1v1 = AKVInnerProduct(v[1],v[1],Ricci,sb)*sqrt(2.)*M_PI;
          v0v1 = AKVInnerProduct(v[0],v[1],Ricci,sb)*sqrt(2.)*M_PI;
          //if(v1v1<symmetry_tol) //symmetries++;
          std::cout << "<v_1|v_1> = " << v1v1 << std::endl;
          std::cout << "<v_0|v_1> = " << v0v1 << std::endl;
          std::cout << "-THETA <v_1|v_1> = " << -THETA[a]*v1v1 << std::endl;
          if(fabs(v0v0) == fabs(v1v2)) badAKVSolution = true;
          break;
        case 2:
          //compute inner products <v_2|v_2>, <v_0|v_2>, <v_1|v_2>
          v2v2 = AKVInnerProduct(v[2],v[2],Ricci,sb)*sqrt(2.)*M_PI;
          v0v2 = AKVInnerProduct(v[0],v[2],Ricci,sb)*sqrt(2.)*M_PI;
          v1v2 = AKVInnerProduct(v[1],v[2],Ricci,sb)*sqrt(2.)*M_PI;
          //if(v2v2<symmetry_tol) //symmetries++;
          std::cout << "<v_2|v_2> = " << v2v2 << std::endl;
          std::cout << "<v_0|v_2> = " << v0v2 << std::endl;
          std::cout << "<v_1|v_2> = " << v1v2 << std::endl;
          std::cout << "-THETA <v_2|v_2> = " << -THETA[a]*v2v2 << std::endl;
          if(fabs(v0v0) == fabs(v0v2)) badAKVSolution = true;
          if(fabs(v1v1) == fabs(v1v2)) badAKVSolution = true;
          break;
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
      KillingDiagnostics(sb, L, Psi, xi[a], rad, printDiagnostic);

      //rotate v, Psi for analysis
      rotated_v[a] = RotateOnSphere(v[a],theta,phi,
                                          sb,thetap[a],phip[a]);

      DataMesh rotated_Psi = RotateOnSphere(Psi,theta,phi,
                                            sb,thetap[a],phip[a]);

      //compare scale factors

      const double scaleAtEquator =
                NormalizeAKVAtOnePoint(sb, rotated_Psi, rotated_v[a], rad, M_PI/2., 0.0);
      PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,rotated_v[a],scaleAtEquator,rad);
      MyVector<double> scaleInnerProduct = InnerProductScaleFactors(v[a], v[a], Ricci, r2p4, sb);
      PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,rotated_v[a],scaleInnerProduct[0],rad);
      PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,rotated_v[a],scaleInnerProduct[1],rad);
      PrintSurfaceNormalization(sb,rotated_Psi,theta,phi,rotated_v[a],scaleInnerProduct[2],rad);
      OptimizeScaleFactor(rotated_v[a], rotated_Psi, rad, sb, theta,
         phi, scaleAtEquator, scaleInnerProduct[0], scaleInnerProduct[1], scaleInnerProduct[2]);

      //scale v
      v[a] *= scaleAtEquator;

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
  }


  // Return 0 for success
  return NumberOfTestsFailed;
}
