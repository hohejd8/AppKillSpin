//=======================================
// $Id: FlatspaceCKV.hpp 2011-09-09 hohejd8 $
//=======================================

///
/// /file
/// Defines ComputeItem::FlatspaceCKV

#ifndef INCLUDED_FlatspaceCKV_hpp
#define INCLUDED_FlatspaceCKV_hpp

#include "Utils/DataBox/ComputeItem.hpp"
#include "Utils/Tensor/Tensor.hpp"
#include "Utils/DataMesh/DataMesh.hpp"
#include "Utils/MyContainers/MyVector.hpp"

namespace ComputeItems {
  /// Adds a conformally flat Killing vector to the DataBox.

  class FlatspaceCKV: public ComputeItem<Tensor<DataMesh> >,
              InstantiateDataBoxAdder<FlatspaceCKV> {

    public:
      static std::string ClassID() {return "FlatspaceCKV"; }
      static std::string Help() {
        return
        "FlatspaceCKV:                                \n"
        "  Places a copy of a conformally flat Killing vector   \n"
        "  in the DataBox.                                      \n"
        "Options:                                               \n"
        "  Metric = Tensor<DataMesh>>; # metric in DataBox      \n"
        "  SurfaceWithMesh = StrahlkorperWithMesh;              \n"
        "       #the surface on which the Killing vectors lie   \n"
        "  Output = MyVector<string>;                           \n"
        "         # name of Killing vector in DataBox           \n"
        "         must choose either translational              \n"
        "         'T_x', 'T_y', 'T_z' or rotational (spin)      \n"
        "         'S_x', 'S_y', 'S_z'                           \n"
        "                                                       \n"
        "Requires in the DataBox:                               \n"
        "  Tensor<DataMesh> [Metric]                            \n"
        "  StrahlkorperWithMesh                                 \n"
        "                                                       \n"
        "Presents to DataBox:                                   \n"
        "  Strahlkorper (surface) of input parameters.          \n";
      };

      FlatspaceCKV(const std::string& opts);
      std::string Output()         const {return mOutput;}
      const result_type& GetData() const {return *mResult;}
      void RecomputeData(const DataBoxAccess& box) const;
    private:
      std::string mMetric, mSurface;
      std::string mOutput;
      MyVector<std::string> mlistKV;
      mutable result_type* mResult;

  }; //class FlatspaceCKV
} //namespace ComputeItems



#endif
