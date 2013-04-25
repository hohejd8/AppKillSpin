//=======================================
// $Id: ComputeAKV.hpp
//=======================================

///
/// /file
/// Defines ComputeItems::ComputeAKV.hpp

#ifndef INCLUDED_ComputeAKV_hpp
#define INCLUDED_ComputeAKV_hpp

#include "Utils/DataBox/ComputeItem.hpp"
#include "SurfaceFinder/Strahlkorper/SurfaceBasis.hpp"
#include "Utils/DataMesh/DataMesh.hpp"

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
        "  AKVSolution = name of approximate Killing vector     \n"
        "                solution in DataBox.                   \n"
        "  WithRicciScaling = bool, true for using Ricci scalar \n"
        "                scaling with the Lagrange multiplier   \n"
        "                THETA.  Default true.                  \n"
        "  StrahlkorperWithMesh = name of StrahlkorperWithMesh  \n"
        "                         surface in DataBox.  Required.\n"
        "  ConformalFactor = conformal factor in DataBox.       \n"
        "                    Required.                          \n"
        "  InterpolateConformalFactor = boolean, true if the    \n"
        "         conformal factor needs to be interpolated     \n"
        "         onto the Strahlkorper surface.  False for     \n"
        "         testing purposes.  Default true.              \n"
        "  AKVGuess = initial guess for THETA, thetap   and     \n"
        "             phip of primary symmetry axis.            \n"
        "             Default {0.0,0.0,0.0}.                    \n"
        "  Radius = the radius of the surface.  Required.       \n"
        "  Solver = which multidimensional root-finding solver  \n"
        "           from the GSL library will be used.          \n"
        "           Options: Hybrids, Hybrid, Newton, Broyden.  \n"
        "           Default Newton.                             \n"
        "  Verbose = Prints THETA, thetap and phip solutions, as\n"
        "            well as all other diagnostics from the     \n"
        "            solution finder.  Default false.           \n"
        "  PrintResiduals = Prints the ic10, ic1p, ic1m         \n"
        "            residuals at each iteration of the solver. \n"
        "            Default to false.                          \n"
        "  ResidualSize = determines the tolerance for ic10,    \n"
        "            ic1p, and ic1m residuals from the          \n"
        "            multidimensional root finder.              \n"
        "            Default to 1.e-10.                         \n"
        "  L_resid_tol = the tolerance for the variation in L   \n"
        "                in the multidimensional root finder.   \n"
        "                Default 1.e-12.                        \n"
        "  v_resid_tol = the tolerance for the variation in v   \n"
        "                in the multidimensional root finder.   \n"
        "                Default 1.e-12.                        \n"
        "  min_thetap = for values less than this, thetap is    \n"
        "               considered close to zero. Default 1.e-5.\n"
        "  PrintTtpSolution = print the solution(s) for THETA,  \n"
        "                     thetap, and phip.  Default true.  \n"
        "  PrintInnerProducts = boolean to print all the inner  \n"
        "               product values for the v solutions in   \n"
        "               each direction.  Default false.         \n"
        "  ScaleFactor = determine which method for finding the \n"
        "                scale factor for v is to be used.      \n"
        "                Options: Equator, InnerProduct1,       \n"
        "                InnerProduct2, InnerProduct3, Optimize.\n"
        "                See AKVSolver::InnerProductScaleFactors\n"
        "                for details on InnerProduct methods.   \n"
        "                Default Equator.                       \n"
        "  PrintScaleFactor = prints the chosen scale factor.   \n"
        "                     Default false.                    \n"
        "  PrintSurfaceNormalization = prints the RMS deviation \n"
        "                   from being properly scaled.  (Best  \n"
        "                   return value is zero.)              \n"
        "                   Default false.                      \n"
        "  PrintSteps = Prints each step along the Killing path.\n"
        "               Default false.                          \n"
        "  PrintBisectionResults = Prints each step of the      \n"
        "            minimization routine for finding the       \n"
        "            optimal scaling factor.  Default false.    \n"
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
        "  FindPoles = find the maximum and minimum values of   \n"
        "             the v DataMesh. Necessary for computing   \n"
        "             the Mobius transform. When false, the     \n"
        "             corresponding rotation uses (theta', phi')\n"
        "             and its antipode. Default false.          \n"
        "                                                       \n"
        "Requires in the DataBox:                               \n"
        "  StrahlkorperWithMesh surface                         \n"
        "  Conformal Factor                                     \n"
        "                                                       \n"
        "Presents to DataBox:                                   \n"
        "  Tensor<DataMesh> [AKVSolution]                       \n";
      };

      ComputeAKV(const std::string& opts);
      std::string Output()         const {return mAKVSolution;}
      const result_type& GetData() const {return *mResult;}
      void RecomputeData(const DataBoxAccess& box) const;

    private:
      std::string mAKVSolution;
      std::string mSkwm, mConformalFactor, mSolver, mScaleFactor;


      MyVector<double> mAKVGuess;
      double mRad, mL_resid_tol, mv_resid_tol, mMin_thetap, mResidualSize;
      double mTestTheta, mTestPhi, mTestEqTheta, mTestEqPhi;
      std::string mDivNorm, mVortNorm, mSS, mfLNorm, mfLambdaNorm, mXiDivLNorm;

      MyVector<bool> printDiagnostic;
      bool mVerbose, mInterpolateConformalFactor, mPrintResiduals,
           mPrintTtpSolution, mPrintInnerProducts, mWithRicciScaling,
           mPrintScaleFactor, mPrintSurfaceNormalization, mPrintSteps,
           mPrintBisectionResults, mFindPoles, mPrintAllFormsOfv;

      mutable result_type* mResult;
  }; //class ComputeAKV
} //namespace ComputeItems



#endif
