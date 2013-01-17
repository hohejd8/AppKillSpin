//============================================
// $Id: ComputeStrahlkorperWithMesh.cpp 2011-09-09 hohejd8
//============================================

#include "ComputeStrahlkorperWithMesh.hpp"
#include "SurfaceFinder/Strahlkorper/Strahlkorper.hpp"
#include "SurfaceFinder/Strahlkorper/StrahlkorperMesh.hpp"
#include "Utils/StringParsing/OptionParser.hpp"

namespace ComputeItems {

  ComputeStrahlkorperWithMesh::ComputeStrahlkorperWithMesh(const std::string& opts)
    : mResult(0)
  {
    OptionParser p(opts, Help());
    mSurface       = p.Get<std::string>("Surface");
    mMesh      = p.Get<std::string>("Mesh");
    mOutput     = p.Get<std::string>("SKWMName");

  }

  //==========================================================================
  
  void ComputeStrahlkorperWithMesh::RecomputeData(const DataBoxAccess& boxa) const {
    delete mResult;
    mResult = new StrahlkorperWithMesh(boxa.Get<StrahlkorperMesh>(mMesh),
                                       boxa.Get<Strahlkorper>(mSurface));

  } //RecomputeData

} // namespace ComputeItems

