//============================================
// $Id: ComputeSKWM.cpp 2011-09-09 hohejd8
//============================================

#include "ComputeSKWM.hpp"
#include "SurfaceFinder/Strahlkorper/Strahlkorper.hpp"
#include "SurfaceFinder/Strahlkorper/StrahlkorperMesh.hpp"
#include "Utils/StringParsing/OptionParser.hpp"

namespace ComputeItems {

  ComputeSKWM::ComputeSKWM(const std::string& opts)
    : mResult(0)
  {
    OptionParser p(opts, Help());
    mSurface       = p.Get<std::string>("Surface");
    mMesh      = p.Get<std::string>("Mesh");
    mOutput     = p.Get<std::string>("SKWMName");

  }

  //==========================================================================
  
  void ComputeSKWM::RecomputeData(const DataBoxAccess& box) const {
    delete mResult;
    mResult = new StrahlkorperWithMesh(box.Get<StrahlkorperMesh>(mMesh),
                                       box.Get<Strahlkorper>(mSurface));

  } //RecomputeData

} // namespace ComputeItems

