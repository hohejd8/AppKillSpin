#include <iostream>
#include <cstdlib>
 
#include "Utils/DataMesh/DataMesh.hpp"
#include "Utils/DataMesh/DataMeshNorms.hpp"
#include "SurfaceFinder/Strahlkorper/SurfaceBasis.hpp"
#include "Utils/ErrorHandling/UtilsForTesting.hpp"
#include "SurfaceFinder/Strahlkorper/StrahlkorperMesh.hpp"
//#include "gsl/gsl_multiroots.h"

//These are the functions we need to test
#include "hohejd8/SurfaceBasisExt.hpp"





int main(int /*argc*/, char** /*argv*/) {


  //create Strahlkorper, SurfaceBasis, SurfaceBasisExt objects
  int l=4;
  //int m=2*l+1;

  StrahlkorperMesh skm(l+1, 2*(l+1));
  //print l,m
  std::cout << "(l,m) should be ("<<l<<","<<l<<"); M->2M+1" << std::endl;
  std::cout << "(l,m) = (" << skm.Basis().L() << "," << skm.Basis().M() << ")"<< std::endl;

  SurfaceBasis sb(skm);
  SurfaceBasisExt sbe(skm);

  //copy surface coords into individual DataMeshes
  DataMesh theta = skm.SurfaceCoords()(0);
  DataMesh phi   = skm.SurfaceCoords()(1);

  //check theta, phi Dim
  std::cout << "theta.Dim() = " << theta.Dim() << std::endl;
  std::cout << "phi.Dim() = " << phi.Dim() << std::endl;

  //check theta, phi extents
  std::cout << "theta.Extents() = " << theta.Extents() << std::endl;
  std::cout << "phi.Extents() = " << phi.Extents() << std::endl;

  //initialize result, solution, and difference DataMeshes
  DataMesh resultDm(theta);
  DataMesh solutionDm(theta);
  DataMesh differenceDm(theta);
  DataMesh coefficientMeshDm(sbe.CoefficientMesh());
    std::cout << "CoefficientMeshDm dimension is "
              << coefficientMeshDm.Dim() << std::endl;

  //check resultDm Dim, Extents
  std::cout << "resultDm.Dim() = " << resultDm.Dim() << std::endl;
  std::cout << "resultDm.Extents() = " << resultDm.Extents() << std::endl;

  //initialize result, solution, and difference Tensor<DataMesh>-es
  Tensor<DataMesh> operandTDm(skm.SurfaceCoords());
  Tensor<DataMesh> resultTDm(skm.SurfaceCoords());
  Tensor<DataMesh> solutionTDm(skm.SurfaceCoords());
  Tensor<DataMesh> differenceTDm(skm.SurfaceCoords());

  //check resultTDm Rank, Dim
  std::cout << "resultTDm.Rank() = " << resultTDm.Rank() << std::endl;
  std::cout << "resultTDm.Dim() = " << resultTDm.Dim() << std::endl;
  std::cout << "\n" << std::endl;

  std::cout << "Test results for various functions.  The printed results \n"
               "show the L1Norm for the difference between the numerical \n"
               "result and analytical solution.  A successful test will  \n"
               "print out 0. \n"
            << std::endl;



  //test VectorPhysToSpec function in SurfaceBasisExt
  std::cout <<
         "VectorPhysToSpec w/ SurfaceBasisExt, \n"
         "f(theta,phi)=-sin(theta) thetaHat + sin(theta) phiHat\n"
         << std::endl;
  operandTDm(0)=-sin(theta);
  operandTDm(1)=sin(theta);
  resultTDm(0)=0.0;
  resultTDm(1)=0.0;
  sbe.VectorPhysToSpec(operandTDm(0).Data(),
                       operandTDm(1).Data(),
                       1,
                       resultTDm(0).Data(),
                       resultTDm(1).Data(),
                       1 );

  std::cout << "Result " << "\n"
            << resultTDm(0) << "\n "
            << resultTDm(1) << "\n "
            << std::endl;


  //test VectorInterpAtPoint function in SurfaceBasisExt
  std::cout <<
      "VectorInterpAtPoint w/ SurfaceBasisExt, \n"
      "f(0,0) = (0,0) \n"
      << std::endl;
  MyVector<double> value1(MV::Size(2),0.0);
  MyVector<double> value2 = sbe.VectorInterpAtPoint(operandTDm, value1, 0.,0.);
  std::cout << "value1 " << value1 << "\n"
            << "value2 " << value2 << std::endl;

  
  std::cout <<
      "VectorInterpAtPoint w/ SurfaceBasisExt, \n"
      "f(pi/2,pi/2) = (-1,1) \n"
      << std::endl;
  value2 = sbe.VectorInterpAtPoint(operandTDm, value1, M_PI/2., M_PI/2.);
  std::cout << "value1 " << value1 << "\n"
            << "value2 " << value2 << std::endl;

  std::cout <<
      "VectorInterpAtPoint w/ SurfaceBasisExt, \n"
      "f(pi/4,pi/4) = (-sqrt(2)/2, sqrt(2)/2) \n"
      << std::endl;
  value2 = sbe.VectorInterpAtPoint(operandTDm, value1, M_PI/4., M_PI/4.);
  std::cout << "value1 " << value1 << "\n"
            << "value2 " << value2 << std::endl;










/**/

  //Return success                                                              
  //return u.NumberOfTestsFailed();
  return EXIT_SUCCESS;
}
