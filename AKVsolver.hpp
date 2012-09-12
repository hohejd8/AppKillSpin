#ifndef INCLUDED_AKVsolver_hpp
#define INCLUDED_AKVsolver_hpp

#include "SurfaceFinder/Strahlkorper/SurfaceBasis.hpp"
#include "SurfaceFinder/Strahlkorper/StrahlkorperWithMesh.hpp"
#include "Utils/DataMesh/DataMesh.hpp"
#include "Utils/DataBox/DataBoxAccess.hpp"
#include "gsl/gsl_vector.h"

//function header for use with gsl 1D root finder
double AKVsolver1D(double THETA, void *params);

//function header for use with gsl multidimensional root finder
int AKVsolver(const gsl_vector * x,
                void *params,
                gsl_vector * f);

//performs a 1D root finder for THETA at thetap=0
bool findTHETA_root(struct rparams * p_compare,
                           double& THETA_root,
                           const bool verbose);

void findTtp(struct rparams * p,
             double& THETA,
             double& thetap,
             double& phip,
             const std::string solver,
             const bool verbose);

//determines normalization factor for approximate Killing vector
double normalizeKillingVector(const SurfaceBasis& sb,
                              const DataMesh& Psi,
                              DataMesh& v,
                              const double& rad);

//rotates a DataMesh by an amount (Theta,Phi)
DataMesh RotateOnSphere
         (const DataMesh& collocationvalues,
          const DataMesh& thetaGrid,
          const DataMesh& phiGrid,
          const SurfaceBasis& sbe,
          const double Theta,
          const double Phi);

bool KillingPath(const SurfaceBasis& sb,
                    const DataMesh& Psi_r,
                    const Tensor<DataMesh>& xi,
                    const double& rad,
                    double& t,
                    const double theta);

int PathDerivs(double t, 
         const double y[],
         double f[],
         void *params);

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
  const double& rad;
  const SurfaceBasis& sb;
  const DataMesh& Psi;
  DataMesh& L;
  DataMesh& v;
  const double & L_resid_tol;
  const double & v_resid_tol;
  const bool printResiduals;
};

//a structure required to follow the Killing path around the surface
struct ODEparams{
  const SurfaceBasis& sb;
  const DataMesh& Psi_ha;
  const Tensor<DataMesh>& xi_ha;
  const double& rad;
};

#endif
