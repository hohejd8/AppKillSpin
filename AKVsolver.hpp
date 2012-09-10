#ifndef INCLUDED_AKVsolver_hpp
#define INCLUDED_AKVsolver_hpp

#include "SurfaceFinder/Strahlkorper/SurfaceBasis.hpp"
#include "SurfaceFinder/Strahlkorper/StrahlkorperWithMesh.hpp"
#include "Utils/DataMesh/DataMesh.hpp"
#include "Utils/DataBox/DataBoxAccess.hpp"
#include "gsl/gsl_vector.h"


//function header for use with gsl root finder
int AKVsolver(const gsl_vector * x,
                void *params,
                gsl_vector * f);

double normalizeKillingVector(const SurfaceBasis& sb,
                              const DataMesh& Psi,
                              DataMesh& v,
                              const double& rad);

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

void KillingDiagnostics(const SurfaceBasis& sb,
                        const DataMesh& L,
                        const DataMesh& Psi,
                        const Tensor<DataMesh>& xi,
                        const double& rad,
                        const MyVector<bool>& printDiagnostic);

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

struct ODEparams{
  const SurfaceBasis& sb;
  const DataMesh& Psi_ha;
  //const DataMesh& Psi;
  const Tensor<DataMesh>& xi_ha;
  //const Tensor<DataMesh>& xi;
  const double& rad;
};

#endif
