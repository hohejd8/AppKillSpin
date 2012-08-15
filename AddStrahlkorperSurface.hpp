//=======================================
// $Id: AddStrahlkorperSurface.hpp 2011-09-09 hohejd8 $
//=======================================


///
/// /file
/// Defines DataBoxAdders::AddStrahlkorperSurface

#ifndef INCLUDED_AddStrahlkorperSurface_hpp
#define INCLUDED_AddStrahlkorperSurface_hpp

#include "Utils/DataBox/DataBoxAdder.hpp"

namespace DataBoxAdders {
  /// Adds a pre-defined Strahlkorper surface to the 
  /// DataBox so that routines can hold exact copies
  /// of a given Strahlkorper surface.  May in some ways
  /// have similar functionality to AddStrahlkorperDataBoxes.

  class AddStrahlkorperSurface: public DataBoxAdder,
              private Factory::Register<AddStrahlkorperSurface> {

    public:
      static std::string ClassID() {return "AddStrahlkorperSurface"; }
      static std::string Help() {
        return
        "AddStrahlkorperSurface:                                \n"
        "  Adds a Strahlkorper surface so that routines can     \n"
        "  have exact copies.                                   \n"
        "Options:                                               \n"
        "  Output = string; # name of this surface in DataBox   \n"
        "  L = int; # spherical harmonic                        \n"
        "  M = int; # spherical harmonic                        \n"
        "  Radius = double; # radius of spherical surface       \n"
        "  Center = MyVector<double>; # location of center      \n"
        "                                                       \n"
        "Output quantities:                                     \n"
        "  Strahlkorper (surface) of input parameters.          \n";
      };

      AddStrahlkorperSurface(): DataBoxAdder(ClassID())  {};
      void AddToDataBoxImpl(const DataBoxInserter& box,
                            const std::string& opts) const;
  };
}



#endif
