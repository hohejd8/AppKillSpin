//=======================================
// $Id: ComputeSKWM.hpp 2011-10-23 hohejd8 $
//=======================================

///
/// /file
/// Defines ComputeItems::ComputeAKV.hpp

#ifndef INCLUDED_ComputeAKV_hpp
#define INCLUDED_ComputeAKV_hpp

#include "Utils/DataBox/ComputeItem.hpp"
//#include "hohejd8/SurfaceBasisExt.hpp"
#include "SurfaceFinder/Strahlkorper/SurfaceBasis.hpp"
#include "Utils/DataMesh/DataMesh.hpp"
#include "gsl/gsl_multiroots.h"

namespace ComputeItems {
  /// Computes the approximate Killing vector when given
  /// a StrahlkorperWithMesh surface, the conformal factor,
  /// and an initial guess for \Theta, \theta ', and \phi '. 

  class ComputeAKV: public ComputeItem<Tensor<DataMesh> >,
              InstantiateDataBoxAdder<ComputeAKV> {

    public:
      static std::string ClassID() {return "ComputeAKV"; }
      static std::string Help() {
        return
        "ComputeAKV:                                            \n"
        "  Computes the approximate Killing vector which is     \n"
        "  necessary for computing the quasi-local angular      \n"
        "  momentum of a single object.                         \n"
        "                                                       \n"
        "Options:                                               \n"
        "  StrahlkorperWithMesh = name of surface in DataBox    \n"
        "  ConformalFactor = conformal factor in DataBox        \n"
        "  AKVGuess = initial guess for AKV                     \n"
        "  Output = name of solution in DataBox                 \n"
        "  DivNorm = print the L2 norm of the divergence of the \n"
        "            approximate Killing vector.  Default false.\n"
        "  VortNorm = print the L2 norm of the vorticity of the \n"
        "            approximate Killing vector.  Default false.\n"
        "  SS = print the surface average S_{ij}S^{ij} value.   \n"
        "            Default false.                             \n"
        "  fLNorm = print the L2 norm of div(grad(L))+L.        \n"
        "           Default false.                              \n"
        "  fLambdaNorm = print the L2 norm of                   \n"
        "              grad(Psi^4/(1-2*lap(ln(Psi)) x Grad(L)   \n"
        "              Default false.                           \n"
        "  XiDivLNorm = print the norm of the approximate       \n"
        "             Killing vector time div(L). Default false.\n"
        "                                                       \n"
        "Requires in the DataBox:                               \n"
        "  StrahlkorperWithMesh surface                         \n"
        "  Conformal Factor                                     \n"
        "                                                       \n"
        "Presents to DataBox:                                   \n"
        "  Tensor<DataMesh> [Output]                            \n";
      };

      ComputeAKV(const std::string& opts);
      std::string Output()         const {return mOutput;}
      const result_type& GetData() const {return *mResult;}
      void RecomputeData(const DataBoxAccess& box) const;

    private:
      void print_state (size_t iter, gsl_multiroot_fsolver * s) const;
      std::string mSkwm;
      std::string mConformalFactor;
      MyVector<double> mAKVGuess;
      std::string mDivNorm, mVortNorm, mSS, mfLNorm, mfLambdaNorm, mXiDivLNorm;
      MyVector<bool> printDiagnostic;
      bool mVerbose;
      std::string mOutput;
      mutable result_type* mResult;


  }; //class ComputeAKV
} //namespace ComputeItems



#endif
