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

  string dataFile = "BHSd8.5_Plus.txt";

  //things I need for this simulation
  //psi, radius, theta, phi, ntheta, nphi

  const int mNth = 14;
  const int mNph = 28;
  //double theta[mNth] = 0.0;
  //double phi[mNph] = 0.0;
  double psi[mNth*mNph] = 0.0;
  double 



  return EXIT_SUCCESS;
}
