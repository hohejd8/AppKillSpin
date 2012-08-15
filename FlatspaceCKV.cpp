//============================================
// $Id: FlatspaceCKV.cpp 2011-09-09 hohejd8
//============================================

#include "FlatspaceCKV.hpp"
#include "Utils/ErrorHandling/Require.hpp"
#include "SurfaceFinder/Strahlkorper/StrahlkorperWithMesh.hpp"
#include "SurfaceFinder/StrahlkorperDataSupplier/StrahlkorperParallelInterpolation.hpp"
#include "Utils/StringParsing/OptionParser.hpp"
#include "Dust/Domain/Domain.hpp"

namespace ComputeItems {

  FlatspaceCKV::FlatspaceCKV(const std::string& opts)
    //: mOutput(0)
  {
    OptionParser p(opts, Help());
    mMetric      = p.Get<std::string>("Metric");
    mSurface     = p.Get<std::string>("SurfaceWithMesh");
    mOutput      = p.Get<std::string>("Output");

    bool mGoodInput = false;
    MyVector<std::string> listKV(MV::fill, "T_x","T_y","T_z","S_x","S_y","S_z");
    mlistKV = listKV;

    for(int i=0; i<listKV.Size(); i++){
      if(mOutput==listKV[i]) mGoodInput=true;
    }

    REQUIRE(mGoodInput, "Killing vector type (Output) not a match." << Help());

  }

  //==========================================================================
  
  void FlatspaceCKV::RecomputeData(const DataBoxAccess& box) const {
    typedef Tensor<DataMesh> TDm;

    const StrahlkorperWithMesh& swm = box.Get<StrahlkorperWithMesh>(mSurface);

    //get the interpolated metric
    std::string mSpatialCoordMap="";
    std::string mInterpolator="Spectral";
    MyVector<std::string> mTensorsToInterp(MV::fill, mMetric);
    const StrahlkorperDataSuppliers::StrahlkorperParallelInterpolation
           supplier(box.Get<Domain>("Domain"), box,
                    box.Get<Domain>("Domain").Communicator(), mTensorsToInterp,
                    mSpatialCoordMap, mInterpolator);
    const TDm& g = supplier.Supply(swm, mMetric);

    //construct the conformal Killing vectors from the interpolated metric
    const int dim = g.Dim();
    const Mesh mesh(g(0,0));
    const TDm zero(dim, "1", mesh, 0.);

    //inertialX with respect to center of excision sphere
    TDm inertialXC(zero);
    TDm inertialX(zero);
    const MyVector<DataMesh>& X = swm.GlobalCoords();
    for(int i=0; i<dim; i++){
      inertialX(i) = X[i];
      inertialXC(i) = inertialX(i) - swm.Center()[i];
      //inertialXC(i) = swm.GlobalCoords()[i] - swm.Center();
    }

    //Xi (translational Killing vectors) and XiRot (rotational Killing vectors)
    //The ith 1-form is \partial_j x^i
    MyVector<TDm> inertialXi(MV::fill,zero,zero,zero);
    MyVector<TDm> inertialRotXi(MV::fill,zero,zero,zero);

      inertialXi[0](0) = 1.;  inertialXi[0](1) = 0.;  inertialXi[0](2) = 0.;
      inertialXi[1](0) = 0.;  inertialXi[1](1) = 1.;  inertialXi[1](2) = 0.;
      inertialXi[2](0) = 0.;  inertialXi[2](1) = 0.;  inertialXi[2](2) = 1.;

      for(int i=0;i<dim;++i){
        inertialRotXi[0](i)=
            inertialXC(1)*inertialXi[2](i)-inertialXC(2)*inertialXi[1](i);
        inertialRotXi[1](i)=
            inertialXC(2)*inertialXi[0](i)-inertialXC(0)*inertialXi[2](i);
        inertialRotXi[2](i)=
            inertialXC(0)*inertialXi[1](i)-inertialXC(1)*inertialXi[0](i);
      }

    //delete mResult;
    if(mOutput==mlistKV[0]) mResult = new TDm(inertialXi[0]);
    if(mOutput==mlistKV[1]) mResult = new TDm(inertialXi[1]);
    if(mOutput==mlistKV[2]) mResult = new TDm(inertialXi[2]);
    if(mOutput==mlistKV[3]) mResult = new TDm(inertialRotXi[0]);
    if(mOutput==mlistKV[4]) mResult = new TDm(inertialRotXi[1]);
    if(mOutput==mlistKV[5]) mResult = new TDm(inertialRotXi[2]);


  } //RecomputeData

} // namespace ComputeItems

