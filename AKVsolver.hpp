#ifndef INCLUDED_AKVsolver_hpp
#define INCLUDED_AKVsolver_hpp

#include "SurfaceFinder/Strahlkorper/SurfaceBasis.hpp"
#include "SurfaceFinder/Strahlkorper/StrahlkorperWithMesh.hpp"
#include "Utils/DataMesh/DataMesh.hpp"
#include "Utils/DataBox/DataBoxAccess.hpp"
//#include "gsl/gsl_vector.h"
#include "gsl/gsl_multiroots.h"


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
//(l,m)=(1,m) residuals are < 1.e-12
bool findTHETA(struct rparams * p,
               double& THETA_root,
               const double& residual_size,
               const bool verbose);

//uses the gsl multidimensional root finder to find
//values for THETA, thetap, phip such that
//(l,m)=(1,m) residuals are < 1.e-12
void findTtp(struct rparams * p,
             double& THETA,
             double& thetap,
             double& phip,
             const std::string solver,
             const double& residual_size,
             const bool verbose);

//function header for use with gsl 1D root finder
double AKVsolver1D(double THETA, void *params);

//function header for use with gsl multidimensional root finder
int AKVsolver(const gsl_vector * x,
                void *params,
                gsl_vector * f);


void print_state (size_t iter, gsl_multiroot_fsolver * s);

//rotates a DataMesh by an amount (Theta,Phi)
DataMesh RotateOnSphere
         (const DataMesh& collocationvalues,
          const DataMesh& thetaGrid,
          const DataMesh& phiGrid,
          const SurfaceBasis& sb,
          const double Theta,
          const double Phi);

//calls normalizeKVAtOnePoint for every point in the mesh,
//then returns the average of all the scale factors
double normalizeKVAtAllPoints(const SurfaceBasis& sb,
                              const DataMesh& Psi,
                              const DataMesh& theta,
                              const DataMesh& phi,
                              const DataMesh& v,
                              const double& rad);

//determines the path length of following the approximate
//Killing vector around the sphere starting at (theta, phi)
//and returns the ratio of that path to the expected
//value of 2*Pi (the normalization scale factor)
double normalizeKVAtOnePoint(const SurfaceBasis& sb,
                              const DataMesh& Psi,
                              //const DataMesh& theta,
                              //const DataMesh& phi,
                              const DataMesh& v,
                              const double& rad,
                              const double& thetap=M_PI/2.,
                              const double& phip=0.0);
double normalizeKVAtOnePoint(const SurfaceBasis& sb,
                              const DataMesh& Psi,
                              //const DataMesh& theta,
                              //const DataMesh& phi,
                              const Tensor<DataMesh>& xi,
                              const double& rad,
                              const double& thetap=M_PI/2.0,
                              const double& phip=0.0);

//This routine will print the value of normalizeKVAtAllPoints for 
//a range of scale factors in the neighborhood of scaleFactor1, scaleFactor2
void TestScaleFactors(const DataMesh& rotated_v,
                      const DataMesh& rotated_Psi,
                      const double& rad,
                      const SurfaceBasis& sb,
                      const DataMesh& theta,
                      const DataMesh& phi,
                      const double& scaleFactor1,
                      const double& scaleFactor2);

bool KillingPath(const SurfaceBasis& sb,
                    const DataMesh& Psi_r,
                    const Tensor<DataMesh>& xi,
                    const double& rad,
                    double& t,
                    const double& theta,
                    const double& phi=0.0,
                    const bool printSteps=false);

int PathDerivs(double t_required_by_solver, 
               const double y[],
               double f[],
               void *params);

//determines value of the integral
// \frac{1}{2} \oint ^2R((_s \vec \nabla v_1) \cdot (_s \vec \nabla v_2) d\Omega
//and returns the ratio of this result to (8*\pi / 3)
double AKVInnerProduct(const DataMesh& v1,
                       const DataMesh& v2,
                       const DataMesh& Ricci,
                       const DataMesh& rp2,
                       const SurfaceBasis& sb);

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
};

//a structure required to follow the Killing path around the surface
struct ODEparams{
  const SurfaceBasis& sb;
  const DataMesh& Psi_ha;
  const Tensor<DataMesh>& xi_ha;
  const double& rad;
};

#endif
