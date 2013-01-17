//=======================================
// $Id: ComputeStrahlkorperWithMesh.hpp 2011-09-09 hohejd8 $
//=======================================

///
/// /file
/// Defines ComputeItems::ComputeStrahlkorperWithMesh.hpp

#ifndef INCLUDED_ComputeStrahlkorperWithMesh_hpp
#define INCLUDED_ComputeStrahlkorperWithMesh_hpp

#include "Utils/DataBox/ComputeItem.hpp"
#include "SurfaceFinder/Strahlkorper/StrahlkorperWithMesh.hpp"
#include "Utils/MyContainers/MyVector.hpp"

namespace ComputeItems {
  /// Adds a pre-defined Strahlkorper surface to the 
  /// DataBox so that routines can hold exact copies
  /// of a given Strahlkorper surface.  May in some ways
  /// have similar functionality to AddStrahlkorperDataBoxes.

  class ComputeStrahlkorperWithMesh: public ComputeItem<StrahlkorperWithMesh>,
              InstantiateDataBoxAdder<ComputeStrahlkorperWithMesh> {

    public:
      static std::string ClassID() {return "ComputeStrahlkorperWithMesh"; }
      static std::string Help() {
        return
        "ComputeStrahlkorperWithMesh:                           \n"
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

      ComputeStrahlkorperWithMesh(const std::string& opts);
      std::string Output()         const {return mOutput;}
      const result_type& GetData() const {return *mResult;}
      void RecomputeData(const DataBoxAccess& boxa) const;
    private:
      std::string mSurface, mMesh, mOutput;
      mutable result_type* mResult;

  }; //class ComputeStrahlkorperWithMesh
} //namespace ComputeItems



#endif
