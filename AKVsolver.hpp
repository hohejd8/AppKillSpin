#ifndef INCLUDED_AKVsolver_hpp
#define INCLUDED_AKVsolver_hpp

#include "SurfaceFinder/Strahlkorper/SurfaceBasis.hpp"
#include "SurfaceFinder/Strahlkorper/StrahlkorperWithMesh.hpp"
#include "Utils/DataMesh/DataMesh.hpp"
#include "Utils/DataMesh/ComplexDataMesh.hpp"
#include "Utils/DataBox/DataBoxAccess.hpp"
//#include "gsl/gsl_vector.h"
#include "gsl/gsl_multiroots.h"
#include "gsl/gsl_multimin.h"

//Is the Killing path centered on the axis after the appropriate (theta, phi) rotation?
//If not, return false *and* return the relative position of the axis
//that the Killing path seems to be centered on.
bool IsKillingPathCentered(const SurfaceBasis& sb,
                           double& thetaOffAxis,
                           double& phiOffAxis,
                           const DataMesh& theta,
                           const DataMesh& phi,
                           const DataMesh& rotated_Psi,
                           const DataMesh& rotated_v,
                           const double& rad,
                           const bool& printSteps=false);

//prints the RMS deviation from a perfectly scaled surface
void PrintSurfaceNormalization(const SurfaceBasis& sb,
                      const DataMesh& rotated_Psi,
                      const DataMesh& theta,
                      const DataMesh& phi,
                      const DataMesh& rotated_v,
                      const double& scaleFactor,
                      const double& rad,
                      const bool& printSteps);

//this function determines which AKVsolver to run based on initial guess
void RunAKVsolvers(double& THETA,
                   double& thetap,
                   double& phip,
                   const double& min_thetap,
                   const double& residualSize,
                   const bool& verbose,
                   struct rparams * p,
                   const std::string solver);

//performs a 1D root finder for THETA at thetap=0 (phip=0)
//returns true if THETA root is found such that
//(l,m)=(1,m) residuals are within tolerance
bool FindTHETA(struct rparams * p,
               double& THETA_root,
               const double& residual_size,
               //const bool verbose);
               const bool verbose,
               const double thetap=0.,
               const double phip=0.);

//uses the gsl multidimensional root finder to find
//values for THETA, thetap, phip such that
//(l,m)=(1,m) residuals are within tolerance
void FindTtp(struct rparams * p,
             double& THETA,
             double& thetap,
             double& phip,
             const std::string solver,
             const double& residual_size,
             const bool verbose);

//function header for use with gsl 1D root finder
//runs AKVsolver at thetap, phip = 0
double AKVsolver1D(double THETA,
                   void *params);
                   //void *params,
                   //const double thetap=0.,
                   //const double phip=0.);

//similar to AKVsolver1D, this function runs AKVsolver once
//at thetap, phip != 0
void EvaluateAKV(double THETA, double thetap, double phip, void *params);

//function header for use with gsl multidimensional root finder
int AKVsolver(const gsl_vector * x,
                void *params,
                gsl_vector * f);

//prints the state of the multidimensional root finder
void print_state (size_t iter, gsl_multiroot_fsolver * s);

//helper function required for GSL multidimensional minimizer routine
double InterpolateDataMesh(const gsl_vector *v,
                           void *params);

//returns the theta, phi components of the extrema for a given DataMesh
void DataMeshExtrema(const DataMesh& collocationvalues,
                     const DataMesh& thetaGrid,
                     const DataMesh& phiGrid,
                     const SurfaceBasis& sb,
                     MyVector<double>& minPoint,
                     MyVector<double>& maxPoint);

//create the set of three complex z values that determine the transform
MyVector<std::complex<double> > MobiusTransformPoints
(
const double nPoleTheta,
const double nPolePhi,
const double sPoleTheta,
const double sPolePhi,
double eqTheta=-1., 
double eqPhi=-1.
);

//perform Mobius transform at one point
std::complex<double> Mobius(const MyVector<std::complex<double> > z,
                            const std::complex<double> z4);

//compute the Mobius conformal factor
double MobiusConformalFactor(const MyVector<std::complex<double> > z,
                             const std::complex<double> z4);

//Mobius transformation on a sphere
DataMesh MobiusTransform(const DataMesh& collocationvalues,//unused right now
                         const DataMesh& thetaGrid,
                         const DataMesh& phiGrid,
                         const SurfaceBasis& sb,
                         const MyVector<std::complex<double> > z,
                         DataMesh& thetaMobius,
                         DataMesh& phiMobius);

//rotates a DataMesh by an amount (Theta,Phi)
DataMesh RotateOnSphere
         (const DataMesh& collocationvalues,
          const DataMesh& thetaGrid,
          const DataMesh& phiGrid,
          const SurfaceBasis& sb,
          const double Theta,
          const double Phi);

//provides a simple mapping between (theta, phi) and (thetapp, phipp) coordinates
//through axis rotation of (thetap, phip).  Similar to RotateOnSphere
void CoordinateRotationMapping(const DataMesh& thetaGrid,
                               const DataMesh& phiGrid,
                               const double& thetap,
                               const double& phip);

//This function computes the approximate Killing vector xi (1-form) given
//the scalar quantity v on the surface
Tensor<DataMesh> ComputeXi(const DataMesh& v, const SurfaceBasis& sb);

//calls normalizeKVAtOnePoint for every point in the mesh,
//then returns the average of all the scale factors
double NormalizeAKVAtAllPoints(const SurfaceBasis& sb,
                              const DataMesh& Psi,
                              const DataMesh& theta,
                              const DataMesh& phi,
                              const DataMesh& v,
                              const double& rad,
                              const bool& printSteps);

//determines the path length of following the approximate
//Killing vector around the sphere starting at (theta, phi)
//and returns the ratio of that path to the expected
//value of 2*Pi (the normalization scale factor)
double NormalizeAKVAtOnePoint(const SurfaceBasis& sb,
                              const DataMesh& Psi,
                              const DataMesh& v,
                              const double& rad,
                              const double& thetap=M_PI/2.,
                              const double& phip=0.0,
                              const bool& printSteps=false);
double NormalizeAKVAtOnePoint(const SurfaceBasis& sb,
                              const DataMesh& Psi,
                              const Tensor<DataMesh>& xi,
                              const double& rad,
                              const double& thetap=M_PI/2.0,
                              const double& phip=0.0,
                              const bool& printSteps=false);

//determines the optimal scale factor using a modified bisection
//routine.  "Optimal" means the chosen scale factor returns the
//closest value to 2*pi for all paths along the surface
double OptimizeScaleFactor(const DataMesh& rotated_v,
                      const DataMesh& rotated_Psi,
                      const double& rad,
                      const SurfaceBasis& sb,
                      const DataMesh& theta,
                      const DataMesh& phi,
                      const double& scaleFactor1,
                      const double& scaleFactor2,
                      const double& scaleFactor3,
                      const double& scaleFactor4,
                      const bool& printSteps,
                      const bool& printBisectionResults);

//integrates along a particular Killing path on the surface
bool KillingPath(const SurfaceBasis& sb,
                    const DataMesh& Psi_r,
                    const Tensor<DataMesh>& xi,
                    const double& rad,
                    double& t,
                    const double& theta,
                    const double& phi,
                    double& thetaOffAxis,
                    double& phiOffAxis,
                    const bool printSteps=false);
int PathDerivs(double t_required_by_solver, 
               const double y[],
               double f[],
               void *params);

//determines value of the integral
// \frac{1}{2} \oint ^2R((_s \vec \nabla v_1) \cdot (_s \vec \nabla v_2) d\Omega
double AKVInnerProduct(const DataMesh& v1,
                       const DataMesh& v2,
                       const DataMesh& Ricci,
                       const SurfaceBasis& sb,
                       const bool& withRicciScaling);

//returns the scaling factors related to various forms of the AKVInnerProduct
MyVector<double> InnerProductScaleFactors(const DataMesh& v1,
                       const DataMesh& v2,
                       const DataMesh& Ricci,
                       const DataMesh& r2p4,
                       const SurfaceBasis& sb,
                       const bool& withRicciScaling);

//returns the proper area integral
//   \frac{1}{\sqrt{2} \pi} \oint 1 dA
// = \frac{1}{\sqrt{2} \pi} \oint r^2 Psi^4 d\Omega
double SurfaceArea(const DataMesh& r2p4,
                   const SurfaceBasis& sb);

//function to perform G-S orthogonalization on a (flexible) DataMesh
//based on its inner product with a (fixed) DataMesh
void GramSchmidtOrthogonalization(const DataMesh& fixedMesh,
                                  const double fixedInnerProduct,
                                  DataMesh& flexibleMesh,
                                  const double crossTermInnerProduct);

//this function will create an initial guess for the next axis of
//symmetry based on previous solutions.  It requires theta, phi for
//prior axes solutions, and an index which indicates whether this is
//the first, second, or third guess
void AxisInitialGuess(double theta[], double phi[], const int index);

//performs diagnostics on the approximate Killing vector solution
void KillingDiagnostics(const SurfaceBasis& sb,
                        const DataMesh& L,
                        const DataMesh& Psi,
                        const Tensor<DataMesh>& xi,
                        const double& rad,
                        const MyVector<bool>& printDiagnostic);

//a structure required to do the gsl multidimensional root finding
struct rparams{
  const DataMesh& theta;
  const DataMesh& phi;
  const DataMesh& rp2;
  const SurfaceBasis& sb;
  const DataMesh& llncf;
  const Tensor<DataMesh>& GradRicci;
  DataMesh& L;
  DataMesh& v;
  const double & L_resid_tol;
  const double & v_resid_tol;
  const bool & printResiduals;
  const bool & ricciScaling;
};

//a test structure required to do 1D root finding with thetap, phip != 0
struct rparam1D{
  struct rparams& p;
  const double& thetap;
  const double& phip;
};

//a structure required to follow the Killing path around the surface
struct ODEparams{
  const SurfaceBasis& sb;
  const DataMesh& Psi_ha;
  const Tensor<DataMesh>& xi_ha;
  const double& rad;
};

//a structure to perform the GSL multidimensional minimizing routine
struct interpparams{
  const SurfaceBasis& sb;
  const DataMesh& collocationValues;
};

#endif
