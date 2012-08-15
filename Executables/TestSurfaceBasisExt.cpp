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
  //int l=100;
  int l=13; //compare to Cook routines
  //int m=2*l+1;

  StrahlkorperMesh skm(l+1, 2*(l+1));
  //print l,m
  std::cout << "(l,m) is ("<<l<<","<<l<<"); M->2M+1" << std::endl;
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
  MyVector<Tensor<DataMesh> > resultVTDm(MV::Size(2),skm.SurfaceCoords());
  Tensor<DataMesh> solutionTDm(skm.SurfaceCoords());
  Tensor<DataMesh> differenceTDm(skm.SurfaceCoords());

  //check resultTDm Rank, Dim
  std::cout << "resultTDm.Rank() = " << resultTDm.Rank() << std::endl;
  std::cout << "resultVTDm[0].Rank() = " << resultVTDm[0].Rank() << std::endl;
  std::cout << "resultVTDm[1].Rank() = " << resultVTDm[1].Rank() << std::endl;
  std::cout << "resultTDm.Dim() = " << resultTDm.Dim() << std::endl;
  std::cout << "\n" << std::endl;

  std::cout << "Test results for various functions.  The printed results \n"
               "show the L1Norm for the difference between the numerical \n"
               "result and analytical solution.  A successful test will  \n"
               "print out 0. \n"
            << std::endl;

  //test the grad function in SurfaceBasis
  //UtilsForTesting u;

  std::cout << "Gradient w/ SurfaceBasis f(theta,phi)=cos(theta)" << std::endl;
  resultTDm=sb.Gradient(cos(theta));
  solutionTDm(0)=-sin(theta);
  solutionTDm(1)=0.0;
  differenceTDm(0)=resultTDm(0)-solutionTDm(0);
  differenceTDm(1)=resultTDm(1)-solutionTDm(1);
  if(L1Norm(differenceTDm(0))<1.e-10 && L1Norm(differenceTDm(1))<1.e-10){
    std::cout << "Test passed \n" << std::endl;
  } else {
    std::cout << "L1Norm(0) = " <<L1Norm(differenceTDm(0)) << "\n"
              << "L1Norm(1) = " <<L1Norm(differenceTDm(1)) << "\n"
              << "Difference " << "\n"
              << differenceTDm(0) << "\n "
              << differenceTDm(1) << "\n "
              << "Original Result " << "\n"
              << resultTDm(0) << "\n "
              << resultTDm(1) << "\n "
              << std::endl;
  }

  std::cout << "Gradient w/ SurfaceBasis f(theta,phi)=-sin(theta)*cos(phi)" << std::endl;
  resultTDm=sb.Gradient(-sin(theta)*cos(phi));
  solutionTDm(0)=-cos(theta)*cos(phi);
  solutionTDm(1)=sin(phi);
  differenceTDm(0)=resultTDm(0)-solutionTDm(0);
  differenceTDm(1)=resultTDm(1)-solutionTDm(1);
  if(L1Norm(differenceTDm(0))<1.e-10 && L1Norm(differenceTDm(1))<1.e-10){
    std::cout << "Test passed \n" << std::endl;
  } else {
    std::cout << "L1Norm(0) = " << L1Norm(differenceTDm(0)) << "\n"
              << "L1Norm(1) = " << L1Norm(differenceTDm(1)) << "\n"
              << "Difference " << "\n"
              << differenceTDm(0) << "\n "
              << differenceTDm(1) << "\n "
              << "Original Result " << "\n"
              << resultTDm(0) << "\n "
              << resultTDm(1) << "\n "
              << std::endl;
  }


  //test the grad function in SurfaceBasisExt
  std::cout <<
         "Gradient w/ SurfaceBasisExt, collocation mesh, f(theta,phi)=cos(theta)"
         << std::endl;
  resultTDm=sbe.Gradient(cos(theta));
  solutionTDm(0)=-sin(theta);
  solutionTDm(1)=0.0;
  differenceTDm(0)=resultTDm(0)-solutionTDm(0);
  differenceTDm(1)=resultTDm(1)-solutionTDm(1);
  if(L1Norm(differenceTDm(0))<1.e-10 && L1Norm(differenceTDm(1))<1.e-10){
    std::cout << "Test passed \n" << std::endl;
  } else {
    std::cout << "L1Norm(0) = " <<L1Norm(differenceTDm(0)) << "\n"
              << "L1Norm(1) = " <<L1Norm(differenceTDm(1)) << "\n"
              //<< "Difference " << "\n"
              //<< differenceTDm(0) << "\n "
              //<< differenceTDm(1) << "\n "
              //<< "Original Result " << "\n"
              //<< resultTDm(0) << "\n "
              //<< resultTDm(1) << "\n "
              << std::endl;
  }
  std::cout <<
         "Gradient w/ SurfaceBasisExt, collocation mesh, f(theta,phi)=-sin(theta)*cos(phi)"
         << std::endl;
  resultTDm=sbe.Gradient(-sin(theta)*cos(phi));
  solutionTDm(0)=-cos(theta)*cos(phi);
  solutionTDm(1)=sin(phi);
  differenceTDm(0)=resultTDm(0)-solutionTDm(0);
  differenceTDm(1)=resultTDm(1)-solutionTDm(1);
  if(L1Norm(differenceTDm(0))<1.e-10 && L1Norm(differenceTDm(1))<1.e-10){
    std::cout << "Test passed \n" << std::endl;
  } else {
    std::cout << "L1Norm(0) = " <<L1Norm(differenceTDm(0)) << "\n"
              << "L1Norm(1) = " <<L1Norm(differenceTDm(1)) << "\n"
              //<< "Difference " << "\n"
              //<< differenceTDm(0) << "\n "
              //<< differenceTDm(1) << "\n "
              //<< "Original Result " << "\n"
              //<< resultTDm(0) << "\n "
              //<< resultTDm(1) << "\n "
              << std::endl;
  }
  std::cout <<
         "Gradient w/ SurfaceBasisExt, coefficient mesh, f(theta,phi)=cos(theta)"
         << std::endl;
  coefficientMeshDm=sbe.ComputeCoefficients(cos(theta));
  resultTDm=sbe.Gradient(coefficientMeshDm);
  solutionTDm(0)=-sin(theta);
  solutionTDm(1)=0.0;
  differenceTDm(0)=resultTDm(0)-solutionTDm(0);
  differenceTDm(1)=resultTDm(1)-solutionTDm(1);
  if(L1Norm(differenceTDm(0))<1.e-10 && L1Norm(differenceTDm(1))<1.e-10){
    std::cout << "Test passed \n" << std::endl;
  } else {
    std::cout << "L1Norm(0) = " <<L1Norm(differenceTDm(0)) << "\n"
              << "L1Norm(1) = " <<L1Norm(differenceTDm(1)) << "\n"
              << "Difference " << "\n"
              << differenceTDm(0) << "\n "
              << differenceTDm(1) << "\n "
              << "Original Result " << "\n"
              << resultTDm(0) << "\n "
              << resultTDm(1) << "\n "
              << std::endl;
  }
  std::cout <<
         "Gradient w/ SurfaceBasisExt, coefficient mesh, f(theta,phi)=-sin(theta)*cos(phi)"
         << std::endl;
  coefficientMeshDm=sbe.ComputeCoefficients(-sin(theta)*cos(phi));
  resultTDm=sbe.Gradient(coefficientMeshDm);
  solutionTDm(0)=-cos(theta)*cos(phi);
  solutionTDm(1)=sin(phi);
  differenceTDm(0)=resultTDm(0)-solutionTDm(0);
  differenceTDm(1)=resultTDm(1)-solutionTDm(1);
  if(L1Norm(differenceTDm(0))<1.e-10 && L1Norm(differenceTDm(1))<1.e-10){
    std::cout << "Test passed \n" << std::endl;
  } else {
    std::cout << "L1Norm(0) = " <<L1Norm(differenceTDm(0)) << "\n"
              << "L1Norm(1) = " <<L1Norm(differenceTDm(1)) << "\n"
              << "Difference " << "\n"
              << differenceTDm(0) << "\n "
              << differenceTDm(1) << "\n "
              << "Original Result " << "\n"
              << resultTDm(0) << "\n "
              << resultTDm(1) << "\n "
              << std::endl;
  }






  //test the Laplacian function in SurfaceBasisExt
  std::cout <<
         "Laplacian w/ SurfaceBasisExt, collocation mesh, f(theta,phi)=cos(theta)"
         << std::endl;
  resultDm=sbe.Laplacian(cos(theta));
  solutionDm=-2.0*cos(theta);
  differenceDm=resultDm-solutionDm;
  if(L1Norm(differenceDm)<1.e-10){
    std::cout << "Test passed \n" << std::endl;
  } else {
    std::cout << "L1Norm() = " <<L1Norm(differenceDm) << "\n"
              << "Difference " << "\n"
              << differenceDm << "\n "
              << "Original Result " << "\n"
              << resultDm << "\n "
              << std::endl;
  }
  std::cout <<
         "Laplacian w/ SurfaceBasisExt, collocation mesh, f(theta,phi)=-sin(theta)*cos(phi)"
         << std::endl;
  resultDm=sbe.Laplacian(-sin(theta)*cos(phi));
  solutionDm=(cos(phi)/sin(theta))*2.0*sin(theta)*sin(theta);
  differenceDm=resultDm-solutionDm;
  if(L1Norm(differenceDm)<1.e-10){
    std::cout << "Test passed \n" << std::endl;
  } else {
    std::cout << "L1Norm() = " <<L1Norm(differenceDm) << "\n"
              << "Difference " << "\n"
              << differenceDm << "\n "
              << "Original Result " << "\n"
              << resultDm << "\n "
              << std::endl;
  }

  //test VectorPhysToSpec function in SurfaceBasisExt
  std::cout <<
         "VectorPhysToSpec w/ SurfaceBasisExt, \n"
         "f(theta,phi)=cos(theta) thetaHat + cos(phi) phiHat \n"
         << std::endl;
  operandTDm(0)=cos(theta);
  operandTDm(1)=cos(phi);
  resultTDm(0)=0.0;
  resultTDm(1)=0.0;
  sbe.VectorPhysToSpec(operandTDm(0).Data(),
                       operandTDm(1).Data(),
                       1,
                       resultTDm(0).Data(),
                       resultTDm(1).Data(),
                       1 );

  std::cout << "Original Result " << "\n"
            << resultTDm(0) << "\n "
            << resultTDm(1) << "\n "
            << std::endl;




  //test the Divergence function in SurfaceBasisExt
  std::cout <<
         "Divergence w/ SurfaceBasisExt, collocation mesh, \n"
         "f(theta,phi)=sin(theta) thetaHat"
         << std::endl;
  operandTDm(0)=sin(theta);
  operandTDm(1)=0.0;
  resultDm=sbe.Divergence(operandTDm);
  solutionDm= 2.0*cos(theta);
  differenceDm=resultDm-solutionDm;
  if(L1Norm(differenceDm)<1.e-10){
    std::cout << "Test passed \n" << std::endl;
  } else {
    std::cout << "L1Norm() = " <<L1Norm(differenceDm) << "\n"
              << "Difference " << "\n"
              << differenceDm << "\n "
              << "Original Result " << "\n"
              << resultDm << "\n "
              << std::endl;
  }
  std::cout <<
         "Divergence w/ SurfaceBasisExt, collocation mesh, \n"
         "f(theta,phi)=-sin(theta)*sin(theta)*cos(phi) phiHat"
         << std::endl;
  operandTDm(0)=0.0;
  operandTDm(1)=-sin(theta)*sin(theta)*cos(phi);
  resultDm=sbe.Divergence(operandTDm);
  solutionDm= sin(theta)*sin(phi);
  differenceDm=resultDm-solutionDm;
  if(L1Norm(differenceDm)<1.e-10){
    std::cout << "Test passed \n" << std::endl;
  } else {
    std::cout << "L1Norm() = " <<L1Norm(differenceDm) << "\n"
              << "Difference " << "\n"
              << differenceDm << "\n "
              << "Original Result " << "\n"
              << resultDm << "\n "
              << std::endl;
  }

  //test the Curl function in SurfaceBasisExt
  std::cout <<
         "Curl w/ SurfaceBasisExt, collocation mesh, \n"
         "f(theta,phi)=-sin(theta) phiHat"
         << std::endl;
  operandTDm(0)= 0.0;
  operandTDm(1)= -sin(theta);
  resultDm=sbe.Curl(operandTDm);//look at the result of the curl; not a vector?
  solutionDm = -2.0*cos(theta);
  differenceDm=resultDm-solutionDm;
  if(L1Norm(differenceDm)<1.e-10){
    std::cout << "Test passed \n" << std::endl;
  } else {
    std::cout << "L1Norm() = " <<L1Norm(differenceDm) << "\n"
              << "Difference " << "\n"
              //<< differenceDm << "\n "
              << "Original Result " << "\n"
              //<< resultDm << "\n "
              << std::endl;
  }

  //test the Curl function in SurfaceBasisExt
  std::cout <<
         "Curl w/ SurfaceBasisExt, collocation mesh, \n"
         "f(theta,phi)= sin(phi) thetaHat"
         << std::endl;
  operandTDm(0)= sin(phi);
  operandTDm(1)= 0.;
  resultDm=sbe.Curl(operandTDm);//look at the result of the curl; not a vector?
  solutionDm = -cos(phi)/sin(theta);
  differenceDm=resultDm-solutionDm;
  if(L1Norm(differenceDm)<1.e-10){
    std::cout << "Test passed \n" << std::endl;
  } else {
    std::cout << "L1Norm() = " <<L1Norm(differenceDm) << "\n"
              << "Difference " << "\n"
              //<< differenceDm << "\n "
              << "Original Result " << "\n"
              //<< resultDm << "\n "
              << std::endl;
  }














/**/

  //Return success                                                              
  //return u.NumberOfTestsFailed();
  return EXIT_SUCCESS;
}
