//=======================================
// $Id: ComputeSKWM.hpp 2011-09-09 hohejd8 $
//=======================================

///
/// /file
/// Defines ComputeItems::ComputeSKWM.hpp

#ifndef INCLUDED_ComputeSKWM_hpp
#define INCLUDED_ComputeSKWM_hpp

#include "Utils/DataBox/ComputeItem.hpp"
#include "SurfaceFinder/Strahlkorper/StrahlkorperWithMesh.hpp"
#include "Utils/MyContainers/MyVector.hpp"

namespace ComputeItems {
  /// Adds a pre-defined Strahlkorper surface to the 
  /// DataBox so that routines can hold exact copies
  /// of a given Strahlkorper surface.  May in some ways
  /// have similar functionality to AddStrahlkorperDataBoxes.

  class ComputeSKWM: public ComputeItem<StrahlkorperWithMesh>,
              InstantiateDataBoxAdder<ComputeSKWM> {

    public:
      static std::string ClassID() {return "ComputeSKWM"; }
      static std::string Help() {
        return
        "ComputeSKWM:                                           \n"
        "  Computes a StrahlkorperWithMesh from a Strahlkorper  \n"
        "  and StrahlkorperMesh present in the DataBox.         \n"
        "                                                       \n"
        "Options:                                               \n"
        "  Surface = Strahlkorper; # Strahlkorper in DataBox    \n"
        "  Mesh = StrahlkorperMesh; # StrahlkorperMesh in DataBox \n"
        "  Output = string;                                     \n"
        "         # name of result                              \n"
        "                                                       \n"
        "Requires in the DataBox:                               \n"
        "  Strahlkorper                                         \n"
        "  StrahlkorperMesh                                     \n"
        "Presents to DataBox:                                   \n"
        "  StrahlkorperWithMesh [Output]                        \n";
      };

      ComputeSKWM(const std::string& opts);
      std::string Output()         const {return mOutput;}
      const result_type& GetData() const {return *mResult;}
      void RecomputeData(const DataBoxAccess& box) const;
    private:
      std::string mSurface, mMesh, mOutput;
      mutable result_type* mResult;

  }; //class ComputeSKWM
} //namespace ComputeItems



#endif
