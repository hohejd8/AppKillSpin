/// Defines ComputeItems::ComputeS2ConformalFactor.hpp

#ifndef INCLUDED_ComputeConformalFactor_hpp
#define INCLUDED_ComputeConformalFactor_hpp

#include "Utils/DataBox/ComputeItem.hpp"
#include "Utils/DataMesh/DataMesh.hpp"
#include "SurfaceFinder/Strahlkorper/StrahlkorperWithMesh.hpp"

namespace ComputeItems {
  /// Computes a pre-set conformal factor for testing purposes

  class ComputeS2ConformalFactor: public ComputeItem<DataMesh>,
              InstantiateDataBoxAdder<ComputeS2ConformalFactor> {

    public:
      static std::string ClassID() {return "ComputeS2ConformalFactor"; }
      static std::string Help() {
        return
        "ComputeS2ConformalFactor:                              \n"
        "  Computes a pre-set conformal factor on S2 for testing\n"
        "  purposes.  Designed for use with ComputeAKV.         \n"
        "                                                       \n"
        "Options:                                               \n"
        "  ConformalFactor = the name of the conformal factor   \n"
        "                    in the DataBox.                    \n"
        "                    Default ConformalFactor.           \n"
        "  Symmetry = Symmetry type for the conformal factor to \n"
        "             be computed.  Options: Spherical, Z, X, Y,\n"
        "             OffAxis, NoSymmetry.                      \n"
        "  StrahlkorperWithMesh = name of StrahlkorperWithMesh  \n"
        "             in the DataBox to setup the conformal     \n"
        "             factor on the surface.                    \n"
        "                                                       \n"
        "Requires in the DataBox:                               \n"
        "  StrahlkorperWithMesh                                 \n"
        "Presents to DataBox:                                   \n"
        "  DataMesh [ConformalFactor]                           \n";
      };

      ComputeS2ConformalFactor(const std::string& opts);
      std::string Output()         const {return mConformalFactor;}
      const result_type& GetData() const {return *mResult;}
      void RecomputeData(const DataBoxAccess& boxa) const;
    private:
      std::string mConformalFactor, mSkWM, mSymmetry;
      mutable result_type* mResult;

  }; //class ComputeS2ConformalFactor
} //namespace ComputeItems



#endif
