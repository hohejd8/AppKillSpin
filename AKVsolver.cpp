#include "AKVsolver.hpp"
#include "Utils/StringParsing/OptionParser.hpp"
#include "Utils/LowLevelUtils/Position.hpp"
#include "Spectral/BasisFunctions/SpherePackIterator.hpp"
#include "gsl/gsl_odeiv2.h"
#include <iomanip>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


//this function determines which AKVsolver to run based on initial guess
void RunAKVsolvers(double& THETA,
                   double& thetap,
                   double& phip,
                   const double& min_thetap,
                   const double& residualSize,
                   const bool& verbose,
                   struct rparams * p,
                   const std::string solver)
{
  //if the initial guess for thetap is close to zero or pi,
  //try solving at thetap = zero
  bool oneDSolutionFound = false;
  const bool thetapGuessIsZero = thetap < min_thetap;
  const bool thetapGuessIsPi = M_PI-thetap < min_thetap;

  if(thetapGuessIsZero){
    oneDSolutionFound = findTHETA(p,THETA,residualSize,verbose);
    if(oneDSolutionFound){
      thetap = 0.0;
      phip = 0.0;
    } else { //thetap initial guess is bad and too close to zero
      thetap = min_thetap;
    }
  }
  if(thetapGuessIsPi){
    oneDSolutionFound = findTHETA(p,THETA,residualSize,verbose);
    if(oneDSolutionFound){
      thetap = M_PI;
      phip = 0.0;
      //findTHETA tests for thetap=0.  If thetap=Pi, v and L can be
      //corrected by multiplying by -1
      p->v *= -1.; p->L *= -1.;
    } else { //thetap initial guess is bad and too close to Pi
      thetap = M_PI - min_thetap;
    }
  }

  //if theta=0 was not a good solution
  //or initial guess was not close to zero
  //try the multidimensional root finder
  if(!oneDSolutionFound){
    findTtp(p, THETA, thetap, phip, solver,residualSize, verbose);

    //if thetap solution is close to zero,
    //and we didn't already try thetap=0,
    //try thetap=0 now
    if( thetap < min_thetap && !thetapGuessIsZero){
      const double THETA_saved = THETA;
      const DataMesh v_saved = p->v;
      const DataMesh L_saved = p->L;
      oneDSolutionFound = findTHETA(p,THETA,residualSize,verbose);
      if(oneDSolutionFound){
        thetap = 0.0;
        phip = 0.0;
      } else {
        THETA = THETA_saved;
        p->v = v_saved; p->L = L_saved;
      } //end if(oneDSolution)
    } // end if(thetap < min_thetap)

    //if thetap solution is close to Pi,
    //and we didn't already try thetap=0,
    //try thetap=0 now
    if( M_PI - thetap < min_thetap && !thetapGuessIsPi){
      const double THETA_saved = THETA;
      const DataMesh v_saved = p->v;
      const DataMesh L_saved = p->L;
      oneDSolutionFound = findTHETA(p,THETA,residualSize,verbose);
      if(oneDSolutionFound){
        thetap = M_PI;
        phip = 0.0;
      } else {
        THETA = THETA_saved;
        p->v = v_saved; p->L = L_saved;
      } //end if(oneDSolution)
    } // end if(thetap < min_thetap)
  } // end if(!oneDSolutionFound)

  //get thetap, phip within standard bounds
  if(thetap < 0.0){
    thetap = -thetap;
    phip -= M_PI;
  }
  if(thetap > M_PI){
    const int m = (thetap/M_PI);
    thetap -= m*M_PI;
    if(m%2) phip -= M_PI;
  }
  if(phip <= -M_PI){
    const int m = (M_PI - phip)/(2.0*M_PI);
    phip += 2.0*m*M_PI;
  } else if(phip > M_PI) {
    const int m = (phip + M_PI)/(2.0*M_PI);
    phip -= 2.0*m*M_PI;
  }

}


//this function attempts to find THETA assuming thetap=0
bool findTHETA(struct rparams * p,
               double& THETA_root,
               const double& residual_size,
               const bool verbose)
{
  std::cout << "Staring the gsl 1D root finder at thetap = 0.0." << std::endl;
  bool goodSolution = false;
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double THETA_lo = -1.0, THETA_hi = 1.0;
  gsl_function F;

  F.function = &AKVsolver1D;
  F.params = p;

  T = gsl_root_fsolver_brent; //bisection, falsepos; should probably just keep brent
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, THETA_lo, THETA_hi);

  do{
    iter++;
    status = gsl_root_fsolver_iterate (s);
    THETA_root = gsl_root_fsolver_root (s);
    THETA_lo = gsl_root_fsolver_x_lower (s);
    THETA_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (THETA_lo, THETA_hi, 0.0, 0.001);

    if (status == GSL_SUCCESS && verbose)
      std::cout << "Converged" << std::endl;

  } while (status == GSL_CONTINUE && iter < max_iter);
     
  gsl_root_fsolver_free (s);

  //compare the residual from multidimensional root finder to 1D root finder
  //run AKVsolver one more time to get residual (f) output
  gsl_vector * x = gsl_vector_alloc(3);
  gsl_vector_set(x,0,THETA_root);
  gsl_vector_set(x,1,0.0);
  gsl_vector_set(x,2,0.0);
  gsl_vector * f = gsl_vector_alloc(3);
  gsl_vector_set(f,0,0.0);
  gsl_vector_set(f,1,0.0);
  gsl_vector_set(f,2,0.0);

  AKVsolver(x, p, f);

  if(verbose)
    std::cout <<  "f0 = " << gsl_vector_get(f,0) 
              << " fp = " << gsl_vector_get(f,1) 
              << " fm = " << gsl_vector_get(f,2) << std::endl;

  //if residuals are small, this is a good solution
  if(   fabs(gsl_vector_get(f,0)) < residual_size
     && fabs(gsl_vector_get(f,1)) < residual_size
     && fabs(gsl_vector_get(f,2)) < residual_size ) goodSolution = true;

  return goodSolution;
}

//this function attempts to find THETA, thetap and phip
void findTtp(struct rparams * p,
             double& THETA,
             double& thetap,
             double& phip,
             const std::string solver,
             const double& residual_size,
             const bool verbose)
{
  std::cout << "Starting the gsl multidimensional root finder." << std::endl;
  //setup the gsl_multiroot finder
  const gsl_multiroot_fsolver_type *T; //solver type
  gsl_multiroot_fsolver *s; //the actual solver itself

  int status;
  size_t iter=0;
  const size_t n = 3; //number of dimensions
  gsl_multiroot_function f = {&AKVsolver, n, p}; //initializes the function
  gsl_vector *x = gsl_vector_alloc(n); //creates initial guess vector
  gsl_vector_set (x, 0, THETA);
  gsl_vector_set (x, 1, thetap);
  gsl_vector_set (x, 2, phip);

  //Declare the appropriate non-derivative root finder
  if(solver=="Hybrids") T = gsl_multiroot_fsolver_hybrids;
  else if(solver=="Hybrid") T = gsl_multiroot_fsolver_hybrid;
  else if(solver=="Newton") T = gsl_multiroot_fsolver_dnewton;
  else if(solver=="Broyden") T = gsl_multiroot_fsolver_broyden;
  else std::cout << "Solver option '" << solver << "' not valid." << std::endl;

  s = gsl_multiroot_fsolver_alloc(T, n);

  gsl_multiroot_fsolver_set(s, &f, x);
  if(verbose) print_state(iter, s);

  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate(s);
    if(verbose) print_state(iter, s);
    if(status){ //if solver is stuck
      std::cout << "GSL multiroot solver is stuck at iter = " << iter << std::endl;
      break;
    }
    status = gsl_multiroot_test_residual(s->f, residual_size);
  } while(status == GSL_CONTINUE && iter<1000);

  if(iter==1000){
    std::cout << "Iteration was stopped at iter=1000. \n"
                 "You may want to check the solution for validity."
              << std::endl;
  }

  THETA  = gsl_vector_get(s->x,0);
  thetap = gsl_vector_get(s->x,1);
  phip   = gsl_vector_get(s->x,2);
  gsl_vector_free(x); //frees all memory associated with vector x
  gsl_multiroot_fsolver_free(s); //frees all memory associated with solver

  return;
}

double AKVsolver1D(double THETA, void *params)
{
  gsl_vector * x = gsl_vector_alloc(3);
  gsl_vector_set(x,0,THETA);
  gsl_vector_set(x,1,0.0);
  gsl_vector_set(x,2,0.0);

  gsl_vector * f = gsl_vector_alloc(3);
  gsl_vector_set(f,0,0.0);
  gsl_vector_set(f,1,0.0);
  gsl_vector_set(f,2,0.0);

  AKVsolver(x, params, f);

  return gsl_vector_get(f,0);
}

//function header required for use with gsl multidimensional root finders
int AKVsolver(const gsl_vector * x,
              void *params,
              gsl_vector * f)
{
  const double THETA = gsl_vector_get (x,0);
  const double thetap = gsl_vector_get (x,1);
  const double phip = gsl_vector_get (x,2);

  //for use with gsl root finder
  const DataMesh& theta = static_cast<struct rparams*>(params)->theta;
  const DataMesh& phi = static_cast<struct rparams*>(params)->phi;
  const DataMesh& rp2 = static_cast<struct rparams*>(params)->rp2;
  const SurfaceBasis& sb = static_cast<struct rparams*>(params)->sb;
  const DataMesh& llncf = static_cast<struct rparams*>(params)->llncf;
  const Tensor<DataMesh>& GradRicci = static_cast<struct rparams*>(params)->GradRicci;
  DataMesh& L = static_cast<struct rparams*>(params)->L;
  DataMesh& v = static_cast<struct rparams*>(params)->v;
  const double& L_resid_tol = static_cast<struct rparams*>(params)->L_resid_tol;
  const double& v_resid_tol = static_cast<struct rparams*>(params)->v_resid_tol;
  const bool& printResiduals = static_cast<struct rparams*>(params)->printResiduals;

  SpherePackIterator sit(theta.Extents()[0],theta.Extents()[1]);
  // index for l=0
  const int l00a = sit(0,0,SpherePackIterator::a);
  // indices for l=1
  const int l10a = sit(1,0,SpherePackIterator::a);
  const int l11a = sit(1,1,SpherePackIterator::a);
  const int l11b = sit(1,1,SpherePackIterator::b);

  //eq. 78, 93
  L = cos(thetap)*cos(theta)
        + sin(thetap)*sin(theta)*(cos(phip)*cos(phi) + sin(phip)*sin(phi) );

  //eq. 95
  v = L;

  //eq. 94, make sure only l=1 mode exists by setting all l!=1 modes to zero
  DataMesh L_ha = sb.ComputeCoefficients(L);
  DataMesh v_ha = sb.ComputeCoefficients(v);

  //guarantee only l=1 modes exist in the harmonic analysis
  for(sit.Reset(); sit; ++sit){
    if(sit.l()==1){
      continue;
    } else {
      L_ha[sit()] = 0.0;
      v_ha[sit()] = 0.0;
    }
  }

  //recompute L from the 'fixed' analysis L_ha
  //this is necessary for the RHS (collocation points) in the while loop
  L = sb.Evaluate(L_ha);

  //The main loop
  bool unsolved = true;
  bool refining = false;
  int refine_count = 0; int iter_count = 0; const int iter_max = 50;
  double ic10 = 0.0; double ic1p = 0.0; double ic1m = 0.0;
  DataMesh RHS(DataMesh::Empty);//RHS=right hand side of eq. 97
  DataMesh RHS_ha(sb.CoefficientMesh());
  Tensor<DataMesh> Gradv(2,"1",DataMesh::Empty); //Gradient of v, eq. 97

  while(unsolved){

    //eq. 97, first term: compute Laplacian of L
    RHS = sb.ScalarLaplacian(L_ha);

    //eq. 97, component of third term: compute gradient of v
    Gradv = sb.Gradient(v_ha);

    //compute eq. 97
    RHS = -RHS 
          + (4.0*llncf*L + 0.5*GradRicci(0)*Gradv(0) + 0.5*GradRicci(1)*Gradv(1) -2.0*L)*(1.0-THETA);

    //perform harmonic analysis on RHS
    RHS_ha = sb.ComputeCoefficients(RHS);

    //keep track of l=1 values
    ic10 = RHS_ha[l10a];
    ic1p = RHS_ha[l11a];
    ic1m = RHS_ha[l11b];

    //print residuals
    if(printResiduals){
      std::cout << "ic10 = " << ic10
                << " ic1p = " << ic1p
                << " ic1m = " << ic1m << std::endl;
    }

    //remove the l=1 modes
    RHS_ha[l10a] = 0.0;
    RHS_ha[l11a] = 0.0;
    RHS_ha[l11b] = 0.0;

    //recompute RHS from 'fixed' analysis
    RHS = sb.Evaluate(RHS_ha);

    //compute RHS^2
    RHS *= RHS;

    //evaluate normalization factor
    const double norm_RHS_L = sqrt(sqrt(2.0)*sb.ComputeCoefficients(RHS)[l00a]/4.0);

    //invert (Laplacian + 2)
    for(sit.Reset(); sit; ++sit){
      if(sit.l()!=1){ // l=1 modes are zero
        RHS_ha[sit()] /= 2.0 - sit.l()*(sit.l()+1.0);
      }
    }

    //update harmonic coefficients for L = L_0 + L_1(RHS)
    for(sit.Reset(); sit; ++sit){
      L_ha[sit()] += RHS_ha[sit()];
    }

    //compute L from harmonic coefficients
    L = sb.Evaluate(L_ha);

    //compute laplacian of delta v
    //NOTE: RHS changes definitions here
    RHS = sb.ScalarLaplacian(v_ha);

    RHS = -RHS - 2.0*rp2*rp2*L;

    //remove the l=0 mode from RHS
    RHS_ha = sb.ComputeCoefficients(RHS);
    RHS_ha[l00a] = 0.0;
    RHS = sb.Evaluate(RHS_ha);

    //compute RHS^2
    RHS *= RHS;

    //compute norm
    const double norm_RHS_v = sqrt(sqrt(2.0)*sb.ComputeCoefficients(RHS)[l00a]/4.0);

    //invert (Laplacian + 0)
    for(sit.Reset(); sit; ++sit){
      if(sit.l()!=0){ // l=0 mode is zero
        RHS_ha[sit()] /= -sit.l()*(sit.l()+1.0);
      }
    }

    //update harmonic coefficients for v
    for(sit.Reset(); sit; ++sit){
      v_ha[sit()] += RHS_ha[sit()];
    }

    //end the primary loop
    if(++iter_count > iter_max) unsolved = false;
    if((norm_RHS_L < L_resid_tol) && (norm_RHS_v < v_resid_tol)) {
      refining = true;
      ++refine_count;
    } else if(refining){
      refining = false;
      refine_count = 0;
    }
    if(refine_count > 2) unsolved = false;
  } //end the main while loop

  //return to calling program
  v = sb.Evaluate(v_ha);

  gsl_vector_set (f, 0, ic10);
  gsl_vector_set (f, 1, ic1p);
  gsl_vector_set (f, 2, ic1m);

  return GSL_SUCCESS;
} //end AKVsolver


void print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  std::cout << "iter = " << iter
            << " x = " << std::setprecision(8) << std::setw(10) << gsl_vector_get(s->x,0)
            << " " << std::setprecision(8) << std::setw(10) << gsl_vector_get(s->x,1)
            << " " << std::setprecision(8) << std::setw(10) << gsl_vector_get(s->x,2)
            << " f(x) = " << std::setprecision(8) << std::setw(10) << gsl_vector_get(s->f,0)
            << " " << std::setprecision(8) << std::setw(10) << gsl_vector_get(s->f,1)
            << " " << std::setprecision(8) << std::setw(10) << gsl_vector_get(s->f,2)
            << std::endl << std::flush;
}

//takes a DataMesh on the sphere and rotates it by some amount Theta, Phi
DataMesh RotateOnSphere
         (const DataMesh& collocationvalues,
          const DataMesh& thetaGrid,
          const DataMesh& phiGrid,
          const SurfaceBasis& sb,
          const double Theta,
          const double Phi) 
{
  DataMesh result(sb.CollocationMesh());

  const double Cb = cos(Theta);
  const double Sb = sin(Theta);
  const double Ca = cos(Phi);
  const double Sa = sin(Phi);

  for(int i=0; i<collocationvalues.Size(); i++){
    const double Cpp = cos(phiGrid[i]);
    const double Spp = sin(phiGrid[i]);
    const double Ctp = cos(thetaGrid[i]);
    const double Stp = sin(thetaGrid[i]);

    double CTheta = Cb*Ctp - Sb*Stp*Cpp;

    if(CTheta > 1.0) CTheta = 1.0;  // don't let roundoff error cause problems
    const double STheta = sqrt(1.0 - CTheta*CTheta);
    double newTheta;
    double newPhi;
    if(STheta > 1.e-12) {
      newTheta = acos(CTheta);
      const double SP = (Sb*Sa*Ctp + Ca*Stp*Spp + Cb*Sa*Stp*Cpp)/STheta;
      const double CP = (Sb*Ca*Ctp - Sa*Stp*Spp + Cb*Ca*Stp*Cpp)/STheta;
      newPhi = atan2(SP,CP);
    } else if(CTheta > 0.0) { // evaluate at north pole
      newTheta = 0.0;
      newPhi = 0.0;
    } else { // evaluate at south pole
      newTheta = M_PI;
      newPhi = 0.0;
    }
    result[i] = sb.Evaluate(collocationvalues, newTheta, newPhi);
  }
  return result;
}

double normalizeKVAtAllPoints(const SurfaceBasis& sb,
                              const DataMesh& Psi,
                              const DataMesh& theta,
                              const DataMesh& phi,
                              const DataMesh& v,
                              const double& rad)
{
  //create xi
  Tensor<DataMesh> tmp_xi = sb.Gradient(v);
  Tensor<DataMesh> xi(2,"1",DataMesh::Empty);
  xi(0) = tmp_xi(1);
  xi(1) = -tmp_xi(0);

  //create storage
  double avgScaleFactor = 0.0;
  DataMesh scaleFactor(Mesh(theta.Extents()));

  //for(int i=0; i<theta.Size(); i++)
  for(int i=0; i<scaleFactor.Size(); i++)
      scaleFactor[i] = normalizeKVAtOnePoint(sb,Psi,xi,rad,theta[i],phi[i])-1.;


  //return avgScaleFactor/theta.Size();
  scaleFactor *= scaleFactor;
  const DataMesh r2p4 = rad*rad*Psi*Psi*Psi*Psi;
  const double area = sb.ComputeCoefficients(r2p4)[0];
  //return sqrt(sb.ComputeCoefficients(scaleFactor)[0]/8.);
//std::cout << "Integral (T/(2*pi) - 1) dOmega / 4*pi                     = " 
//          << sqrt(sb.ComputeCoefficients(scaleFactor)[0]/8.) << std::endl;
//std::cout << "Integral (T/(2*pi) - 1) r2p4 dOmega / Integral r2p4 dOmega= " 
//          << sqrt(sb.ComputeCoefficients(scaleFactor*r2p4)[0])/area << std::endl;
  return sqrt(sb.ComputeCoefficients(scaleFactor*r2p4)[0])/area;
}

double normalizeKVAtOnePoint(const SurfaceBasis& sb,
                              const DataMesh& Psi,
                              //const DataMesh& theta,
                              //const DataMesh& phi,
                              const DataMesh& v,
                              const double& rad,
                              const double& thetap, /*=M_PI/2.0*/
                              const double& phip /*=0.0*/)
{
  //create xi
  Tensor<DataMesh> tmp_xi = sb.Gradient(v);
  Tensor<DataMesh> xi(2,"1",DataMesh::Empty);
  xi(0) = tmp_xi(1);
  xi(1) = -tmp_xi(0);

  return normalizeKVAtOnePoint(sb, Psi, xi, rad, thetap, phip);
}

double normalizeKVAtOnePoint(const SurfaceBasis& sb,
                              const DataMesh& Psi,
                              //const DataMesh& theta,
                              //const DataMesh& phi,
                              const Tensor<DataMesh>& xi,
                              const double& rad,
                              const double& thetap, /*=M_PI/2.0*/
                              const double& phip /*=0.0*/)
{

  //rotate Psi
/*
  DataMesh rotated_Psi = RotateOnSphere(Psi,
                         theta,
                         phi,
                         sb,
                         thetap,
                         phip);
*/
  DataMesh rotated_Psi = Psi;

  //Rescale the Killing vector. For a rotational Killing vector,
  //the affine path length should be 2*Pi.  Also, test various
  //paths to ensure we have an actual Killing field
  double t; //affine path length


  bool goodtheta = KillingPath(sb, rotated_Psi, xi, rad, t, thetap, phip);
  REQUIRE(goodtheta, "Killing trajectory did not close " << POSITION);
  const double scale = t/(2.0*M_PI);
/*
  std::cout << "theta = " << std::setprecision(8) << std::setw(10)
            << 0.5*M_PI << " : affine path = "
            << std::setprecision(10) << t //<< std::endl;
            << " Scale factor = " << scale << std::endl;

  KillingPath(sb, rotated_Psi, xi, rad, t, 0.5*M_PI/2.0);
  std::cout << "Theta = " << std::setprecision(8) << std::setw(10)
            << 0.5*M_PI/2.0 << " : affine path = "
            << std::setprecision(10) << t //<< std::endl;
            << "Scale factor = " << t/(2.0*M_PI*scale) << std::endl;

  KillingPath(sb, rotated_Psi, xi, rad, t, 0.25*M_PI);
  std::cout << "Theta = " << std::setprecision(8) << std::setw(10)
            << 0.25*M_PI << " : affine path = "
            << std::setprecision(10) << t //<< std::endl;
            << " Scale factor = " << t/(2.0*M_PI*scale) << std::endl;

  KillingPath(sb, rotated_Psi, xi, rad, t, 0.75*M_PI);
  std::cout << "Theta = " << std::setprecision(8) << std::setw(10)
            << 0.75*M_PI << " : affine path = "
            << std::setprecision(10) << t //<< std::endl;
            << " Scale factor = " << t/(2.0*M_PI*scale) << std::endl;
*/
  return scale;
} //end normalizeKillingVector


//This routine will print the value of normalizeKVAtAllPoints for 
//a range of scale factors in the neighborhood of scaleFactor1, scaleFactor2
void TestScaleFactors(const DataMesh& rotated_v,
                      const DataMesh& rotated_Psi,
                      const double& rad,
                      const SurfaceBasis& sb,
                      const DataMesh& theta,
                      const DataMesh& phi,
                      const double& scaleFactor1,
                      const double& scaleFactor2)
{
  const double difference = fabs(scaleFactor1-scaleFactor2);
  std::cout << "scale factors " << scaleFactor1 << " " << scaleFactor2 << std::endl;
  for(int i=0; i<=400; i++){
    const double deviation = (1.-2.*difference)+(difference/100.)*i;
    double scaleFactor = scaleFactor1*deviation;
    //rotated_v *= scaleFactor;
    const double scaleOverSurface =
                normalizeKVAtAllPoints(sb, rotated_Psi, theta, phi, rotated_v*scaleFactor, rad);
    std::cout << std::setprecision(15) << scaleFactor << " " << scaleOverSurface << std::endl;
  }
}

bool KillingPath(const SurfaceBasis& sb,
                 const DataMesh& Psi,
                 const Tensor<DataMesh>& xi,
                 const double& rad,
                 double& t,
                 const double& theta,
                 const double& phi, /*=0.0*/
                 const bool printSteps /*=false*/)
{
  //bool printSteps = false;
  //perform harmonic analysis on Psi, xi
  //to save from repetitive computation in PathDerivs
  const DataMesh Psi_ha = sb.ComputeCoefficients(Psi);
  const Tensor<DataMesh> xi_ha = sb.ComputeVectorCoefficients(xi);
  struct ODEparams params = {sb, Psi_ha, xi_ha, rad};

  const gsl_odeiv2_step_type *const T = gsl_odeiv2_step_rkck;
  gsl_odeiv2_step *const s = gsl_odeiv2_step_alloc (T, 2);
  gsl_odeiv2_control *const c = gsl_odeiv2_control_y_new (1e-12, 1.e-12);
  gsl_odeiv2_evolve *const e = gsl_odeiv2_evolve_alloc (2);
  gsl_odeiv2_system sys = {PathDerivs, NULL, 2, &params};
     
  t = 0.0;
  double t1 = 1.e10;
  double h = 1e-6;
  double y[2] = {theta, phi-2.*M_PI}; //this makes the stopping criteria easier
  bool limit_h = false;
  double hmax = h;
  const double hmax2 = 2.0*M_PI/100.0;

  if(printSteps){
    std::cout << "START: y = ( " 
	      << std::setprecision(8) << std::setw(10) << y[0] << " , " 
	      << std::setprecision(8) << std::setw(10) << y[1]+2.*M_PI << " ); h = " 
	      << std::setprecision(8) << std::setw(10) << h 
	      << std::endl; 
  }
     
  int iter = 0;//delete this logic later
  while (true && iter<1000) {
    iter++; //delete this logic later
    const double ysave[2] = {y[0],y[1]}; 
    const double tsave = t;
    const int status = gsl_odeiv2_evolve_apply (e, c, s,
						&sys, 
						&t, t1,
						&h, y);
    if(printSteps){
      std::cout << "t = " 
    	        << std::setprecision(8) << std::setw(10) << t 
  	        << " : y = ( " 
  	        << std::setprecision(8) << std::setw(10) << y[0] << " , " 
  	        << std::setprecision(8) << std::setw(10) << y[1]+2.*M_PI << " ); h = " 
  	        << std::setprecision(8) << std::setw(10) << h 
  	        << std::endl;
    }

    ASSERT(status==GSL_SUCCESS,"Path Integration failed");
    if(limit_h && h > hmax) h = hmax;
    if(fabs(y[1] - phi) < 1.e-11) break;
    else if(y[1] > phi) { //if solver went too far...
      if(!limit_h) hmax = h;
      limit_h = true;
      h = hmax *= 0.5; //...reset h, hmax
      y[0] = ysave[0]; y[1] = ysave[1]; t = tsave; //return variables to previous state
      gsl_odeiv2_evolve_reset(e); //return solver to previous state
    } //end ifs
    if(h > hmax2) h = hmax2;
  } //end while

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);

  const bool closedPath = (fabs(y[0] - theta) < 1.e-6);
  if(!closedPath){
    std::cout << "##> Theta diff: " << std::setprecision(6) << y[0] - theta << std::endl;
  }
  return closedPath;
}

int PathDerivs(double t_required_by_solver, const double y[], double f[], void *params)
{
  const SurfaceBasis& sb = static_cast<struct ODEparams*>(params)->sb;
  const DataMesh& Psi_ha = static_cast<struct ODEparams*>(params)->Psi_ha;
  const Tensor<DataMesh>& xi_ha = static_cast<struct ODEparams*>(params)->xi_ha;
  const double& rad = static_cast<struct ODEparams*>(params)->rad;

  double psiAtPoint = sb.EvaluateFromCoefficients(Psi_ha, y[0], y[1]);
  MyVector<double> result = sb.EvaluateVectorFromCoefficients(xi_ha, y[0], y[1]);

  const double norm = 1.0 / (psiAtPoint*psiAtPoint*psiAtPoint*psiAtPoint
                             *rad*rad);
  f[0] = result[0]*norm;
  f[1] = result[1]*norm/sin(y[0]);

  return GSL_SUCCESS;
}

//double* AKVInnerProduct(const DataMesh& v1,
MyVector<double> AKVInnerProduct(const DataMesh& v1,
                       const DataMesh& v2,
                       const DataMesh& Ricci,
                       const DataMesh& rp2,
                       const SurfaceBasis& sb)
{
  const Tensor<DataMesh> Gradv1 = sb.Gradient(v1);
  const Tensor<DataMesh> Gradv2 = sb.Gradient(v2);

  DataMesh integrand = 0.5*Ricci*(Gradv1(0)*Gradv2(0)+Gradv1(1)*Gradv2(1));//(rp2*rp2);
  double area = sb.ComputeCoefficients(rp2*rp2)[0];

  double ip1 = (3.*sqrt(2.)/8.0)*sb.ComputeCoefficients(integrand/(rp2*rp2))[0];
  //std::cout << "ip1 : " 
  //          << 1./ip1 << std::endl;

  double ip2 = sb.ComputeCoefficients(integrand)[0] / (2./3. * area);
  //std::cout << "ip2 : " 
  //          << 1./sqrt(ip2) << std::endl;
  //std::cout << "ip3 : "
  //          << 1./ip2 << std::endl;

  //double ip[3]={1./ip1, 1./sqrt(ip2), 1./ip2};

  return MyVector<double>(MV::fill,1./ip1, 1./sqrt(ip2), 1./ip2) ;
}


void KillingDiagnostics(const SurfaceBasis& sb,
                        const DataMesh& L,
                        const DataMesh& Psi,
                        const Tensor<DataMesh>& xi,
                        const double& rad,
                        const MyVector<bool>& printDiagnostic)
{
  const DataMesh& rp2 = rad*Psi*Psi;
  const double rp22_00 = sb.ComputeCoefficients(rp2*rp2)[0];//(0,0) component of coefficients

  DataMesh div = sb.Divergence(xi)/rp2;

  //-----Divergence-------
  if(printDiagnostic[0]){
    //const DataMesh& div2norm = sb.ComputeCoefficients(div*div);
    std::cout << "L2 Norm of Divergence = "
              << std::setprecision(12)
              //<< sqrt(sqrt(2.)*div2norm[0]/4.) << std::endl;
              << sb.ComputeCoefficients(div*div)[0]/rp22_00 << std::endl;
  }

  if(printDiagnostic[1] || printDiagnostic[2]){
    DataMesh vort = sb.Vorticity(xi)/rp2 - 2.0*rp2*L;

    //-----Vorticity-------
    if(printDiagnostic[1]){
      //const DataMesh vort2norm = sb.ComputeCoefficients(vort*vort);
      std::cout << "L2 Norm of Vorticity = "
                << std::setprecision(12)
                //<< sqrt(sqrt(2.)*vort2norm[0]/4.) << std::endl;
                << sb.ComputeCoefficients(vort*vort)[0]/rp22_00 << std::endl;
    }

    //-----SS-------
    if(printDiagnostic[2]){
      Tensor<DataMesh> Dxtheta = xi;
      Tensor<DataMesh> Dxphi = xi;

      Dxtheta = sb.VectorColatitudeDerivative(xi);

      Tensor<DataMesh> gradlncf = sb.Gradient(log(Psi));


      Dxtheta(0) -= 2.0*( xi(0)*gradlncf(0) - xi(1)*gradlncf(1) );
      Dxtheta(1) -= 2.0*( xi(0)*gradlncf(1) + xi(1)*gradlncf(0) );
      Dxphi(0) = Dxtheta(1) - vort;
      Dxphi(1) = div - Dxtheta(0);

      DataMesh SS =  Dxtheta(0)*Dxtheta(0)
                   + Dxtheta(1)*Dxtheta(1)
                   + Dxphi(0)*Dxphi(0)
                   + Dxphi(1)*Dxphi(1);

      const DataMesh p4r2 = rad*rad*Psi*Psi*Psi*Psi;
      SS /= p4r2;
      SS -= 2.0*p4r2*L*L;

      DataMesh SS_ha = sb.ComputeCoefficients(SS);

      std::cout << "Surface Average S_{ij}S^{ij} = "
                << std::setprecision(12) << sqrt(2.)*SS_ha[0]/4. << std::endl;
    }
  }

  if(printDiagnostic[3] || printDiagnostic[4] || printDiagnostic[5]){
    Tensor<DataMesh> GradL = sb.Gradient(L);
    DataMesh xiGradL = xi(0)*GradL(0) + xi(1)*GradL(1);
    const DataMesh& llncf = sb.ScalarLaplacian(log(Psi));
    //const DataMesh Rm1 = Psi*Psi*Psi*Psi / (1.-2.*llncf);
    const DataMesh Rm1 = (rp2*rp2/2.0)/(1.-2.*llncf);

    //-----Norm of f_L-------
    if(printDiagnostic[4]){
      Tensor<DataMesh> RGradL(2,"1",DataMesh::Empty);
      RGradL(0) = GradL(0)*Rm1;
      RGradL(1) = GradL(1)*Rm1;
      DataMesh f_L = sb.Divergence(RGradL);
      f_L /= rp2;
      f_L += rp2*L;

      std::cout << "L2 Norm of f_L= "
                << std::setprecision(12)
                //<< sqrt(sqrt(2.)*fL_ha[0]/4.) << std::endl;
                << sb.ComputeCoefficients(f_L*f_L)[0]/rp22_00 << std::endl;
    } //end f_L

    //-----Norm of f_Lambda-------
    if(printDiagnostic[5]){
      const Tensor<DataMesh> GradRm1 = sb.Gradient(Rm1);
      DataMesh f_lam = GradRm1(0)*GradL(1) - GradRm1(1)*GradL(0);
      f_lam /= rp2;

      std::cout << "L2 Norm of f_lam= "
                << std::setprecision(12)
                //<< sqrt(sqrt(2.)*flam_ha[0]/4.) << std::endl;
                << sb.ComputeCoefficients(f_lam*f_lam)[0]/rp22_00 << std::endl;
    } // end f_Lambda

    //-----Norm of xi*GradL-------
    if(printDiagnostic[3]){
      xiGradL = xi(0)*GradL(0) + xi(1)*GradL(1);
      xiGradL /= rp2;

      //const DataMesh xiGradL_ha = sb.ComputeCoefficients(xiGradL*xiGradL);

      std::cout << "L2 Norm of xi*Grad(L) = "
                << std::setprecision(12)
                //<< sqrt(sqrt(2.)*xiGradL_ha[0]/4.) << std::endl;
                << sb.ComputeCoefficients(xiGradL*xiGradL)[0]/rp22_00 << std::endl;
    } //end xi*DivL

  } //end print conditionals
}//end KillingDiagnostics

