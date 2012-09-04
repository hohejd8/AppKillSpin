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


double normalizeKillingVector(void *params,
                              const double thetap,
                              const double phip);

DataMesh RotateOnSphere
         (const DataMesh& collocationvalues,
          const DataMesh& thetaGrid,
          const DataMesh& phiGrid,
          const SurfaceBasis& sbe,
          const double Theta,
          const double Phi);

bool KillingPathNew(void *params);

int func(double t, 
         const double y[],
         double f[],
         void *params);

bool KillingPath(void *params,
                 const DataMesh& Psi, //rotated Psi
                 double& t,
                 double theta,
                 const Tensor<DataMesh>& xi);

MyVector<double> PathDerivs(void *params,
                const DataMesh& Psi, //rotated Psi
                const MyVector<double>& Vin,
                const Tensor<DataMesh>& xi);

void PathRKQC(const double& ds,
              double& h,
              const double& hmax,
              double& hdid,
              MyVector<double>& Vin,
              MyVector<double>& Vp,
              void *params,
              const DataMesh& Psi, //rotated Psi
              const Tensor<DataMesh>& xi,
              MyVector<double>& scale,
              const double epsilon);

void PathRKCK(const double& h,
              const MyVector<double>& Vin,
              MyVector<double>& Vout,
              const MyVector<double>& Vp,
              void *params,
              const DataMesh& Psi, //rotated Psi
              const Tensor<DataMesh>& xi,
              MyVector<double>& error);

void KillingDiagnostics(const StrahlkorperWithMesh& skwm,
                        const DataMesh& L,
                        const DataMesh& Psi,
                        const Tensor<DataMesh>& xi,
                        const MyVector<bool>& printDiagnostic);

struct rparams{
  const StrahlkorperWithMesh& skwm;
  const DataMesh& Psi;
  DataMesh& L;
  DataMesh& v;
  const double L_resid;
  const double v_resid;
  const bool PrintResiduals;
};

struct ODEparams{
  const SurfaceBasis& sb;
  const DataMesh& Psi_r;
};

#endif
