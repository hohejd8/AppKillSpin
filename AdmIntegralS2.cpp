//================================================================
// $Id: AdmIntegralS2.cpp 2011-06-30 hohejd8 based on AdmIntegralsOverS2 $
//================================================================
#include "AdmIntegralS2.hpp"
#include "Utils/Tensor/Tensor.hpp"
#include "Utils/DataMesh/DataMesh.hpp"
#include "Utils/IO/PrependDirectoryName.hpp"
#include "Utils/StringParsing/OptionParser.hpp"
#include "Utils/DataBox/DataBoxAccess.hpp"
#include "Dust/Domain/Domain.hpp"
#include "SurfaceFinder/Strahlkorper/Strahlkorper.hpp"
#include "SurfaceFinder/Strahlkorper/StrahlkorperWithMesh.hpp"
#include "SurfaceFinder/StrahlkorperDataSupplier/StrahlkorperParallelInterpolation.hpp"
#include "Utils/DataBox/DataBox.hpp"
#include "Utils/DataBox/DataBoxInserter.hpp"
#include "Utils/LowLevelUtils/Position.hpp"
#include "Observers/DumpTensors.hpp"
#include "Utils/IO/FileExists.hpp"
#include "Utils/LowLevelUtils/ConvertNumberToString.hpp"
#include "Utils/MiscUtils/SimpleProfiler.hpp"
#include "stdio.h"
#include "Utils/MpiWrappers/MpiCommunicator.hpp"
#include "Utils/IO/CachedOfStream.hpp"
#include <iomanip>
//TESTING
//#include "SurfaceBasisExt.hpp"
//TESTING

namespace Observers {

  AdmIntegralS2::AdmIntegralS2(const std::string &opts,
			       const DataBox& Box,
                               const MpiCommunicator& Comm,
                               const std::string& Dir)
  {
    typedef MyVector<double>      MvD;
    typedef Tensor<DataMesh> TDm;
    mBaseFileName            = "AdmIntegralS2";

    PrependDirectoryInPlace(mBaseFileName, Dir);

std::cout << "AdmIntegralS2 constructor" << std::endl;
    OptionParser p(opts,Help());
    mMetric                  = p.Get<std::string>("Metric", "g");
    mInvMetric               = p.Get<std::string>("InvMetric", "Invg");
    mExtrinsicCurvature      = p.Get<std::string>("ExtrinsicCurvature", "K");
    mPsi                     = p.Get<std::string>("ConformalFactor", "ConformalFactor");
    mCenter                  = p.Get<MyVector<double> >("Center");
    mRadius                  = p.Get<double>("Radius");
    mL                       = p.Get<int>("L");
    mKillingVector           = p.Get<MyVector<std::string> >("KillingVector");
    mStrahlkorperWithMesh    = p.Get<std::string>("StrahlkorperWithMesh");
    mInertialXi = p.Get<std::string>("InertialXi");

    //set mInput, which is the list of integrals to be output
    //MyVector<std::string> listOfIntegrals(MV::fill,"Sx_Adm","Sy_Adm","Sz_Adm");
    //mInput = listOfIntegrals;
    mInput = mKillingVector; //this needs to be traced down and cleaned up
  }

//*****************************************************//
//*****************************************************//
  MyVector<double> AdmIntegralS2::
  DoIntegrals( const DataBox& Box, const MpiCommunicator& Comm ) const  {   
 
    typedef Tensor<DataMesh> TDm;
    DataBoxAccess box(Box,"AdmIntegralS2::DoIntegrals");

    //---------------------------------------
    // Get StrahlkorperWithMesh, interpolate
    //---------------------------------------
    const StrahlkorperWithMesh&
         surfaceWithMesh(box.Get<StrahlkorperWithMesh>(mStrahlkorperWithMesh));//new

    const Domain& D = box.Get<Domain>("Domain");
    // List tensors to interpolate onto surface, then interpolate them
    MyVector<std::string> TensorsToInterp(MV::fill,mMetric,
                           mInvMetric,
                           mExtrinsicCurvature,
                           mPsi);
    std::string mSpatialCoordMap="";
    const StrahlkorperDataSuppliers::StrahlkorperParallelInterpolation
           supplier(D, box, Comm, TensorsToInterp,
                    mSpatialCoordMap, "Spectral");//new
    REQUIRE(supplier.IsDataAvailable(surfaceWithMesh),
             "StrahlkorperDataSupplier: Can't interpolate "
             "3d data onto 'surfaceWithMesh'.\n");

    //---------------------------------------
    // Setup local DataBox to hold interpolated values
    // (see StrahlkorperMetaObserver)
    //---------------------------------------
    //making a DataBox might be a bit overkill for now, but will allow
    //easy transition to parallel
    DataBox localBox("surfaceWithMesh DataBox");
    DataBoxInserter localBoxInserter(localBox, POSITION);
    for(int i=0; i<TensorsToInterp.Size(); i++){
      localBoxInserter.AddVolatileItem(TensorsToInterp[i],
                           supplier.Supply(surfaceWithMesh, TensorsToInterp[i]));
    }

    //for computation/reading, make local (referenced) copies of interpolated tensors
    DataBoxAccess lba(localBox, "AdmIntegralS2::DoIntegral");
    const TDm& g(lba.Get<TDm>(mMetric));
    const TDm& Invg(lba.Get<TDm>(mMetric));
    const TDm& K(lba.Get<TDm>(mExtrinsicCurvature));
    const DataMesh& Psi(lba.Get<TDm>(mPsi)());

    //---------------------------------
    //Setup non-interpolated values
    //---------------------------------
    //Set template dimension, mesh, DataMesh, and TDm
    const int dim=3;
    const Mesh mesh(g(0,0));
    const DataMesh zeroMesh(mesh, 0.);
    const TDm zero(dim, "1", mesh, 0.);

    //alternative method for getting the normalized normal vectors
    TDm radBasis(surfaceWithMesh.RadialBasisVector());

    //ADM angular momentum integrand and integral
    MyVector<DataMesh>
      admAngMomentumIntegrand(MV::Size(mKillingVector.Size()), zeroMesh);//new



//----------------------------------------------
    //SPIN EXPRESSION 9/14/2011
    // S_{(k)} = \frac{1}{8 \pi} K_{ij} ( s_m  \gamma^{m,i} ) \xi_{(k)n}  \Psi^{-6}
    
    for(int k=0;k<mKillingVector.Size();++k){
      const TDm& KV = box.Get<TDm>(mKillingVector[k]);
      for(int i=0;i<dim;++i){
        for(int j=0;j<dim;++j){
          for(int m=0;m<dim;++m){
            admAngMomentumIntegrand[k] += (K(i,j) )
              * radBasis(m) * Invg(m,i) * 1./(Psi*Psi*Psi*Psi)
              * KV(j) * 1./(Psi*Psi);
          }
        }	      
      }
      admAngMomentumIntegrand[k] /= (8.*M_PI);
    }
//-----------------------------------------------


    //-----------------------------------
    //Integration method taken from OnStrahlkorperComputeSurfaceIntegral.cpp
    //-----------------------------------
    MyVector<double> admAngularMomentum(MV::fill,0.,0.,0.);
    const DataMesh& dA = surfaceWithMesh.SurfaceAreaElement(g);
    //TEST FOR SurfaceBasisExt
    //SurfaceBasisExt sbc(surfaceWithMesh.Grid());
    SurfaceBasis sbc(surfaceWithMesh.Grid());
    for(int i=0; i<mKillingVector.Size(); i++){ 
      admAngularMomentum[i] =
           sbc.Integrate(admAngMomentumIntegrand[i]*dA);
    }
//testing only
std::cout << "spin " << admAngularMomentum << std::endl;
//end testing
    return admAngularMomentum;
  }
 
//***************************************************************
//***************************************************************
 
  void AdmIntegralS2::Observe(const DataBox& Box,
         const MpiCommunicator& Comm, 
	 const double /* TimeForDumpFiles */,
	 const std::string& FirstColumnInDatFiles) const {
    SimpleProfiler bla("Observers::"+ClassID());
    DataBoxAccess box(Box,"AdmIntegralS2::Observe");
    const Domain& mDomain=box.Get<Domain>("Domain");


    //Returns the KV integrals from the surface
    MyVector<double> Integrals = DoIntegrals(Box, Comm);

    // Done with integrals. Now output (only on processor zero).
    if(mDomain.Communicator().Rank()==0) {
      std::string FileName = mBaseFileName+".dat";
      bool newfile = !FileExists(FileName);
      CachedOfStream& Out1=GetCachedOfStream(FileName);
      if(newfile) {
	int i = 1;
	Out1 << "# AdmIntegralS2\n";
	Out1 << "# [" << i++ << "] = time \n";
	for(int var=0;var<mInput.Size();++var) {
	    Out1 << "# [" << i++ << "] = " 
		<< mInput[var] << "(R=" 
		<< DoubleToString(mRadius,3) << ")\n";
//for testing only
	    std::cout << "# [" << i++ << "] = " 
		<< mInput[var] << "(R=" 
		<< DoubleToString(mRadius,3) << ")\n";
//end testing
	}
      }
      Out1 << FirstColumnInDatFiles;
      for(int var=0;var<mInput.Size();++var) {
	  Out1 << " " << DoubleToString(Integrals[var],16);
//for testing only
std::cout << " " << DoubleToString(Integrals[var],16);
//end testing
      }
      Out1 << "\n";
//for testing only
std::cout << "\n";
//end testing
    }
  }
}


