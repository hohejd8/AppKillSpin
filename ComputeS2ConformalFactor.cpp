#include "ComputeS2ConformalFactor.hpp"
#include "Utils/StringParsing/OptionParser.hpp"
#include "Utils/ErrorHandling/Require.hpp"
//#include <iomanip>
#include "Utils/LowLevelUtils/Position.hpp"

namespace ComputeItems {

  ComputeS2ConformalFactor::ComputeS2ConformalFactor(const std::string& opts)
    : mResult(0)
  {
    OptionParser p(opts, Help());
    mConformalFactor = p.Get<std::string>("ConformalFactor","ConformalFactor");
    mSkWM = p.Get<std::string>("StrahlkorperWithMesh");
    mSymmetry = p.Get<std::string>("Symmetry");
std::cout << POSITION << std::flush << std::endl;
  }

  void ComputeS2ConformalFactor::RecomputeData(const DataBoxAccess& boxa) const {
    delete mResult;
    const StrahlkorperWithMesh& skwm = boxa.Get<StrahlkorperWithMesh>(mSkWM);
    const SurfaceBasis sb(skwm.Grid().Basis());
    const DataMesh theta = boxa.Get<StrahlkorperWithMesh>(mSkWM).Grid().SurfaceCoords()(0);
    const DataMesh phi = boxa.Get<StrahlkorperWithMesh>(mSkWM).Grid().SurfaceCoords()(1);
std::cout << POSITION << std::flush << std::endl;
    DataMesh Psi(theta); //copy constructor to get the right size

    Psi = 1.28;
std::cout << POSITION << std::flush << std::endl;
    if(mSymmetry=="Spherical"){
      //Do nothing.  Psi is already spherically symmetric.
    } else if(mSymmetry=="Z"){
      Psi += 0.001*cos(theta)*cos(theta);
    } else if(mSymmetry=="X"){
      Psi += 0.001*(1.0-3.0*cos(theta)*cos(theta)+3.0*sin(theta)*sin(theta)*cos(2.0*phi));
    } else if(mSymmetry=="Y"){
      Psi += 0.001*(1.0-3.0*cos(theta)*cos(theta)-3.0*sin(theta)*sin(theta)*cos(2.0*phi));
    } else if(mSymmetry=="OffAxis"){
      //based on Merzbacher Eq. 16.60
      //Y(L,0)(thetap,phip) = sqrt(4*PI/(2L+1))\Sum(M=-L..L) Y(L,M)(theta,phi) Y(L,M)*(beta,alpha)
      //assume z-symmetry in the new coordinate system (thetap, phip), where the new z-axis had
      //coordinates (beta, alpha) in the old system

      Psi+= 0.001
                   *(-1.0+3.0*cos(theta)*cos(theta)
		   +3.0*sqrt(2.0)*sin(2.0*theta)*(cos(phi)+sin(phi))
		   +3.0*sin(theta)*sin(theta)*sin(2.0*phi));
    } else if(mSymmetry=="None"){
      Psi += 0.001*(0.
		   +3.0*sqrt(2.0)*sin(2.0*theta)*(cos(phi)+sin(phi))
		   +3.0*sin(theta)*sin(theta)*sin(2.0*phi));
    } else {
      REQUIRE(false,"Symmetry type not recognized.");
    }
std::cout << POSITION << std::flush << std::endl;
    mResult = new DataMesh(Psi);    
  }
}
