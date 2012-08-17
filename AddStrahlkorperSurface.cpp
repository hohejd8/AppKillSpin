#include "AddStrahlkorperSurface.hpp"
//#include "Utils/DataBox/AddToDataBox.hpp"
#include "Utils/StringParsing/OptionParser.hpp"
#include "SurfaceFinder/Strahlkorper/Strahlkorper.hpp"
#include "SurfaceFinder/Strahlkorper/StrahlkorperMesh.hpp"

namespace DataBoxAdders {

  void AddStrahlkorperSurface::
  AddToDataBoxImpl(const DataBoxInserter& box,const std::string& Opts) const
  {
    OptionParser p(Opts, Help());
    const int mNth = p.Get<int>("Nth");
    const int mNph = p.Get<int>("Nph");
    const double mRadius = p.Get<double>("Radius");
    const MyVector<double> mCenter = p.Get<MyVector<double> >("Center");
    const std::string mSurface = p.Get<std::string>("SurfaceName");
    const std::string mMesh = p.Get<std::string>("MeshName");

    const Strahlkorper theStrahlkorper(2, 2, mRadius, mCenter);
       //Strahlkorpers must be created with L>=2, despite what the direct documentation
       //states.  A Strahlkorper(L,M,radius,center) creates a SurfaceBasis(L,M),
       //which in turn creates a YlmSpherePack(L+1,2*M), with restrictions that
       //(L+1)>=3 and (2*M)>=4.  See Spectral/BasisFunctions/YlmSpherePack.cpp
       //for more information.
    const StrahlkorperMesh sm(mNth, mNph);

    box.AddVolatileItem(mSurface, theStrahlkorper);
    box.AddVolatileItem(mMesh, sm);
  }
}

