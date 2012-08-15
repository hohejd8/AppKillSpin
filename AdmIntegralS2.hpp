//================================================================
// $Id: AdmIntegralS2.hpp 2011-09-27 hohejd8 $
//================================================================
///
/// \file
/// Defines Observers::AdmIntegralS2.
#ifndef INCLUDED_AdmIntegralS2_hpp
#define INCLUDED_AdmIntegralS2_hpp

#include "Dust/Domain/Observer.hpp"
#include "Utils/MyContainers/MyVector.hpp"
#include "Utils/MyContainers/IPoint.hpp"

//#include "ObservationS2Surface.hpp"
#include "SurfaceFinder/Strahlkorper/SurfaceBasis.hpp"


namespace Observers {

  /// Computes volume integral of scalars on S2 surfaces.
  class AdmIntegralS2 : 
           public Observer, InstantiateObserverAdder<AdmIntegralS2>
{
  public: 
    static std::string ClassID() {return "AdmIntegralS2";}
    static std::string    Help() {
      return "AdmIntegralS2                                            \n"
	"Computes the ADM energy, linear momentum, and angular momentum     \n"
	" on inertial spheres. Not intended for use on Inv-mapped spherical \n"
        " shells. (For that case, use the observer AdmIntegrals instead.)   \n"
	"OPTIONS:                                                           \n"
	"   Metric=string; #name of spatial metric                          \n"
        "   dMetric = string; #FLATTENED (indices at end!!) deriv           \n"
	"                     #of metric                                    \n"
	"   InvMetric=string; #name of inverse spatial metric               \n"
	"   ExtrinsicCurvature=string; #name of extrinsic curvature         \n"
	"   MapPrefix = string; # e.g., 'GridToInertial'                    \n"
	"                                                                   \n"
	"   Radius       = double,double,...;      #Can choose several Rs   \n"
	"                                          #If 2 or more are given, \n"
	"                                          #will extrapolate ADM    \n"
	"                                          #quantities to infinity  \n"
	"                                                                   \n"
	"   Center       = double,double,double;   #Must be 3-D             \n"
	"   L            = int;                                             \n"
	"   S2Surface    = string; #default InertialSphere                  \n"
	"            # Radii, center, maximum Ylm L, and type of the S2     \n"
	"            # surfaces (e.g. InertialSphere).                      \n"
	"            # WARNING: unless you know what you are doing, use     \n"
	"            # the default.                                         \n"
	"                                                                   \n"
	"   BaseFileName = string;  # Will have '*.dat' appended.           \n"
        "                           # default AdmIntegrals                  \n"
	;
    }
  public:
    AdmIntegralS2(const std::string &opts,
		  const DataBox& Box,
                  const MpiCommunicator& Comm,
                  const std::string& Dir);

    void Observe(const DataBox& box,
                 const MpiCommunicator& Comm, 
		 const double TimeForDumpFiles,
		 const std::string& FirstColumnInDatFiles) const;

    MyVector<double> DoIntegrals( const DataBox& box,
                 const MpiCommunicator& Comm ) const;    
      
  private:

    MyVector<std::string> mInput;
    std::string mMetric, mInvMetric, mExtrinsicCurvature;
    MyVector<std::string> mKillingVector;
    std::string mStrahlkorperWithMesh;
    std::string mBaseFileName, mInertialCoordinates, mPsi;
    std::string mInertialXi;
    int mL;
    double mRadius;
    MyVector<double> mCenter;
    IPoint mExtents;
    int mSdIndex;
    int mFLAG;
    bool mTerminateOnPointsOutside;
  };
  
}   

#endif
