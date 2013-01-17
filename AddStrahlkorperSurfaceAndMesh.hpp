/// \file
/// Defines DataBoxAdders::AddStrahlkorperSurfaceAndMesh.
#ifndef AddStrahlkorperSurfaceAndMesh_hpp
#define AddStrahlkorperSurfaceAndMesh_hpp

#include "Utils/DataBox/DataBoxAdder.hpp"

namespace DataBoxAdders {
  /// Adds a pre-defined Strahlkorper surface to the 
  /// DataBox so that routines can hold exact copies
  /// of a given Strahlkorper surface.  May in some ways
  /// have similar functionality to AddStrahlkorperDataBoxes.

  class AddStrahlkorperSurfaceAndMesh: public DataBoxAdder,
              private Factory::Register<AddStrahlkorperSurfaceAndMesh> 
  {
    public:
      static std::string ClassID() {return "AddStrahlkorperSurfaceAndMesh";}
      static std::string Help() {
        return ClassID()+"\n"
        "  Adds Strahlkorper Surface, Mesh, and SurfaceWithMesh \n"
        "  so that routines can have exact copies.              \n"
        "Options:                                               \n"
        "  SurfaceName = string; # name of Surface in DataBox   \n"
        "  MeshName = string; # name of Mesh in DataBox         \n"
        "  L = int; # spherical harmonic                        \n"
        "  M = int; # spherical harmonic                        \n"
        "  Radius = double; # radius of spherical surface       \n"
        "  Center = MyVector<double>; # location of center      \n"
        "                                                       \n"
        "Output quantities:                                     \n"
        "  StrahlkorperSurface [SurfaceName]                    \n"
        "  StrahlkorperMesh [MeshName]                          \n";
      };

      AddStrahlkorperSurfaceAndMesh() {};
      void AddToDataBoxImpl(const DataBoxInserter& boxi,
			    const std::string& Opts) const;
  };
}

#endif // AddStrahlkorperSurfaceAndMesh_hpp
