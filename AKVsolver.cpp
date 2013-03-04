#include "AKVsolver.hpp"
#include "Utils/StringParsing/OptionParser.hpp"
#include "Utils/LowLevelUtils/Position.hpp"
#include "Spectral/BasisFunctions/SpherePackIterator.hpp"
#include "gsl/gsl_odeiv2.h"
#include <iomanip>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <stdlib.h>

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
                           const bool& printSteps /*false*/)
{
  double t;
  double thetap = theta[0];
  double phip = phi[0];

  //create xi
  Tensor<DataMesh> tmp_xi = sb.Gradient(rotated_v);
  Tensor<DataMesh> xi(2,"1",DataMesh::Empty);
  xi(0) = tmp_xi(1);
  xi(1) = -tmp_xi(0);

  KillingPath(sb, rotated_Psi, xi, rad, t, thetap, phip,
              thetaOffAxis, phiOffAxis, printSteps);

  if(fabs(thetaOffAxis*180./M_PI)>1.) return false;

  return true;
}

//prints the RMS deviation from a perfectly scaled surface
void PrintSurfaceNormalization(const SurfaceBasis& sb,
                      const DataMesh& rotated_Psi,
                      const DataMesh& theta,
                      const DataMesh& phi,
                      const DataMesh& rotated_v,
                      const double& scaleFactor,
                      const double& rad,
                      const bool& printSteps)
{
  const double scaleOverSurface =
     NormalizeAKVAtAllPoints(sb, rotated_Psi, theta, phi, rotated_v*scaleFactor, rad, printSteps);
  std::cout << std::setprecision(15) 
            << scaleFactor 
            << " " << scaleOverSurface << std::endl;
}

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
  const bool thetapGuessIsPi = M_PI-fmod(thetap,M_PI) < min_thetap;

  if(thetapGuessIsZero){
    const double phip_saved = phip;
    const double THETA_saved = THETA;
    oneDSolutionFound = FindTHETA(p,THETA,residualSize,verbose);
    if(oneDSolutionFound){
      thetap = 0.0;
      phip = 0.0;
    } else { //thetap initial guess is bad and too close to zero
      thetap = min_thetap;
      phip = phip_saved;
      THETA = THETA_saved;
    }
  }
  if(thetapGuessIsPi){
    const double phip_saved = phip;
    const double THETA_saved = THETA;
    oneDSolutionFound = FindTHETA(p,THETA,residualSize,verbose);
    if(oneDSolutionFound){
      thetap = M_PI;
      phip = 0.0;
      //FindTHETA tests for thetap=0.  If thetap=Pi, v and L can be
      //corrected by multiplying by -1
      p->v *= -1.; p->L *= -1.;
    } else { //thetap initial guess is bad and too close to Pi
      thetap = M_PI - min_thetap;
      phip = phip_saved;
      THETA = THETA_saved;
    }
  }

  //if theta=0 was not a good solution
  //or initial guess was not close to zero
  //try the multidimensional root finder
  if(!oneDSolutionFound){
    FindTtp(p, THETA, thetap, phip, solver,residualSize, verbose);

    //if thetap solution is close to zero,
    //and we didn't already try thetap=0,
    //try thetap=0 now
    if( thetap < min_thetap && !thetapGuessIsZero){
      const double THETA_saved = THETA;
      const DataMesh v_saved = p->v;
      const DataMesh L_saved = p->L;
      oneDSolutionFound = FindTHETA(p,THETA,residualSize,verbose);
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
    if( M_PI - fmod(thetap,M_PI) < min_thetap && !thetapGuessIsPi){
      const double THETA_saved = THETA;
      const DataMesh v_saved = p->v;
      const DataMesh L_saved = p->L;
      oneDSolutionFound = FindTHETA(p,THETA,residualSize,verbose);
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

//performs a 1D root finder for THETA at thetap=0 (phip=0)
//returns true if THETA root is found such that
//(l,m)=(1,m) residuals are within tolerance
bool FindTHETA(struct rparams * p,
               double& THETA_root,
               const double& residual_size,
               //const bool verbose)
               const bool verbose,
               const double thetap /*=0.*/,
               const double phip /*=0.*/)
{
  //std::cout << "Starting the gsl 1D root finder at thetap = 0.0." << std::endl;
  std::cout << "Starting the gsl 1D root finder at (thetap, phip) = (" 
            << thetap << ", " << phip << ")." << std::endl;
  bool goodSolution = false;
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double THETA_lo = -1.0, THETA_hi = 1.0;
  gsl_function F;

  //test structure for thetap, phip != 0
  rparam1D p1D = {*p, thetap, phip};

  F.function = &AKVsolver1D;
  //F.params = p;
  F.params = &p1D;

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
  gsl_vector_set(x,1,thetap);
  gsl_vector_set(x,2,phip);
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

//uses the gsl multidimensional root finder to find
//values for THETA, thetap, phip such that
//(l,m)=(1,m) residuals are within tolerance
void FindTtp(struct rparams * p,
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
  else{
    std::cout << "Solver option '" << solver << "' not valid." << std::endl;
    T = gsl_multiroot_fsolver_dnewton;
  }

  s = gsl_multiroot_fsolver_alloc(T, n);

  gsl_multiroot_fsolver_set(s, &f, x);
  if(verbose) print_state(iter, s); //uncomment this later

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

//function header for use with gsl 1D root finder
//runs AKVsolver at thetap, phip = 0
double AKVsolver1D(double THETA,
                   void *params)
                   //void *params,
                   //const double thetap, /*=0.*/
                   //const double phip /*=0.*/)
{
  const double& thetap = static_cast<struct rparam1D*>(params)->thetap;
  const double& phip = static_cast<struct rparam1D*>(params)->phip;
  rparams& p = static_cast<struct rparam1D*>(params)->p;

  gsl_vector * x = gsl_vector_alloc(3);
  gsl_vector_set(x,0,THETA);
  gsl_vector_set(x,1,thetap);
  gsl_vector_set(x,2,phip);

  gsl_vector * f = gsl_vector_alloc(3);
  gsl_vector_set(f,0,0.0);
  gsl_vector_set(f,1,0.0);
  gsl_vector_set(f,2,0.0);

  //AKVsolver(x, params, f);
  AKVsolver(x, &p, f);

  return gsl_vector_get(f,0);
}

//similar to AKVsolver1D, this function runs AKVsolver once
//at thetap, phip != 0
void EvaluateAKV(double THETA, double thetap, double phip, void *params)
{
  gsl_vector * x = gsl_vector_alloc(3);
  gsl_vector_set(x,0,THETA);
  gsl_vector_set(x,1,thetap);
  gsl_vector_set(x,2,phip);

  gsl_vector * f = gsl_vector_alloc(3);
  gsl_vector_set(f,0,0.0);
  gsl_vector_set(f,1,0.0);
  gsl_vector_set(f,2,0.0);

  AKVsolver(x, params, f);
}

//function header for use with gsl multidimensional root finder
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
  const bool& ricciScaling = static_cast<struct rparams*>(params)->ricciScaling;

  SpherePackIterator sit(theta.Extents()[0],theta.Extents()[1]);
  // index for l=0, 1
  const int l00a = sit(0,0,SpherePackIterator::a);
  const int l10a = sit(1,0,SpherePackIterator::a);
  const int l11a = sit(1,1,SpherePackIterator::a);
  const int l11b = sit(1,1,SpherePackIterator::b);

  //valid solution for both ricciScaling = true, false
  L = cos(thetap)*cos(theta)
        + sin(thetap)*sin(theta)*(cos(phip)*cos(phi) + sin(phip)*sin(phi) );
  v = L;

  //make sure only l=1 mode exists by setting all l!=1 modes to zero
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
  DataMesh RHS(DataMesh::Empty);
  DataMesh RHS_ha(sb.CoefficientMesh());
  Tensor<DataMesh> Gradv(2,"1",DataMesh::Empty); //Gradient of v

  while(unsolved){
    //compute Laplacian of L
    RHS = sb.ScalarLaplacian(L_ha);

    //component of third term: compute gradient of v
    Gradv = sb.Gradient(v_ha);

    if(ricciScaling){
      RHS = -RHS 
          + (4.0*llncf*L + 0.5*GradRicci(0)*Gradv(0) + 0.5*GradRicci(1)*Gradv(1) -2.0*L)*(1.0-THETA);
    } else {
      RHS = -RHS
          + 4.0*llncf*L + 0.5*GradRicci(0)*Gradv(0) + 0.5*GradRicci(1)*Gradv(1) -2.0*L +THETA*rp2*rp2*L;
    }

    //perform harmonic analysis on RHS
    RHS_ha = sb.ComputeCoefficients(RHS);

    //keep track of l=1 values
    ic10 = RHS_ha[l10a];
    ic1p = RHS_ha[l11a];
    ic1m = RHS_ha[l11b];

    //print residuals
    if(printResiduals){
      std::cout << "THETA: " << THETA
                << "thetap: " << thetap
                << "phip: " << phip << std::endl;
      std::cout << "ic10: " << ic10
                << " ic1p: " << ic1p
                << " ic1m: " << ic1m << std::endl;
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

//prints the state of the multidimensional root finder
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
//update variable names Theta->thetap, Phi->phip
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

//provides a simple mapping between (theta, phi) and (thetapp, phipp) coordinates
//through axis rotation of (thetap, phip).  Similar to RotateOnSphere
void CoordinateRotationMapping(const DataMesh& thetaGrid,
                               const DataMesh& phiGrid,
                               const double& thetap,
                               const double& phip)
{
  //copy constructor, which can be fixed later for better coding practice
  DataMesh thetaGridNew = thetaGrid;
  DataMesh thetaGridNewNew = thetaGrid;
  DataMesh phiGridNew = phiGrid;
  DataMesh phiGridNewNew = phiGrid;

  {
  const double Cb = cos(thetap);
  const double Sb = sin(thetap);
  const double Ca = cos(phip);
  const double Sa = sin(phip);
  for(int i=0; i<thetaGrid.Size(); i++){
    const double Cpp = cos(phiGrid[i]);
    const double Spp = sin(phiGrid[i]);
    const double Ctp = cos(thetaGrid[i]);
    const double Stp = sin(thetaGrid[i]);

    double CTheta = Cb*Ctp - Sb*Stp*Cpp;

    if(CTheta > 1.0) CTheta = 1.0;  // don't let roundoff error cause problems
    const double STheta = sqrt(1.0 - CTheta*CTheta);
    //double newTheta;
    //double newPhi;
    if(STheta > 1.e-12) {
      thetaGridNew[i] = acos(CTheta);
      const double SP = (Sb*Sa*Ctp + Ca*Stp*Spp + Cb*Sa*Stp*Cpp)/STheta;
      const double CP = (Sb*Ca*Ctp - Sa*Stp*Spp + Cb*Ca*Stp*Cpp)/STheta;
      phiGridNew[i] = atan2(SP,CP);
    } else if(CTheta > 0.0) { // evaluate at north pole
      thetaGridNew[i] = 0.0;
      phiGridNew[i] = 0.0;
    } else { // evaluate at south pole
      thetaGridNew[i] = M_PI;
      phiGridNew[i] = 0.0;
    }
  }
  }

  //now perform the inverse mapping
  {
  const double Cb = cos(thetap);
  const double Sb = sin(thetap);
  const double Ca = cos(-phip);
  const double Sa = sin(-phip);
  for(int i=0; i<thetaGridNew.Size(); i++){
    const double Cpp = cos(phiGridNew[i]);
    const double Spp = sin(phiGridNew[i]);
    const double Ctp = cos(thetaGridNew[i]);
    const double Stp = sin(thetaGridNew[i]);

    double CTheta = Cb*Ctp - Sb*Stp*Cpp;

    if(CTheta > 1.0) CTheta = 1.0;  // don't let roundoff error cause problems
    const double STheta = sqrt(1.0 - CTheta*CTheta);
    //double newTheta;
    //double newPhi;
    if(STheta > 1.e-12) {
      thetaGridNewNew[i] = acos(CTheta);
      //const double SP = (Sb*Sa*Ctp + Ca*Stp*Spp + Cb*Sa*Stp*Cpp)/STheta;
      //const double CP = (Sb*Ca*Ctp - Sa*Stp*Spp + Cb*Ca*Stp*Cpp)/STheta;
      //const double CP = ( 0.*Ctp + Ca*Stp*Spp + Sa*Stp*Cpp)/STheta;
      //const double SP = ( -Sb*Ctp + Cb*Sa*Stp*Spp + Cb*Ca*Stp*Cpp)/STheta;
      const double SP = (Sb*Sa*Ctp + Ca*Stp*Spp + Cb*Sa*Stp*Cpp)/STheta;
      const double CP = (Sb*Ca*Ctp - Sa*Stp*Spp + Cb*Ca*Stp*Cpp)/STheta;

      phiGridNewNew[i] = atan2(SP,CP);
    } else if(CTheta > 0.0) { // evaluate at north pole
      thetaGridNewNew[i] = 0.0;
      phiGridNewNew[i] = 0.0;
    } else { // evaluate at south pole
      thetaGridNewNew[i] = M_PI;
      phiGridNewNew[i] = 0.0;
    }
  }
  }

  std::cout << "Theta Grid" << std::endl;
  std::cout << thetaGrid << std::endl;
    std::cout << "ThetaNew Grid" << std::endl;
  std::cout << thetaGridNew << std::endl;
    std::cout << "ThetaNewNew Grid" << std::endl;
  std::cout << thetaGridNewNew << std::endl;
  std::cout << "Phi Grid" << std::endl;
  std::cout << phiGrid << std::endl;
  std::cout << "PhiNew Grid" << std::endl;
  std::cout << phiGridNew << std::endl;
  std::cout << "PhiNewNew Grid" << std::endl;
  std::cout << phiGridNewNew << std::endl;
}

//This function computes the approximate Killing vector xi (1-form) given
//the scalar quantity v on the surface
Tensor<DataMesh> ComputeXi(const DataMesh& v, const SurfaceBasis& sb)
{
  Tensor<DataMesh> xi(2,"1",DataMesh::Empty);
  Tensor<DataMesh> tmp_xi = sb.Gradient(v);
  xi(0) = tmp_xi(1);
  xi(1) = -tmp_xi(0);
  return xi;
}

//calls NormalizeAKVAtOnePoint for every point in the mesh,
//then returns the average of all the scale factors
double NormalizeAKVAtAllPoints(const SurfaceBasis& sb,
                              const DataMesh& Psi,
                              const DataMesh& theta,
                              const DataMesh& phi,
                              const DataMesh& v,
                              const double& rad,
                              const bool& printSteps)
{
  //create xi
  Tensor<DataMesh> tmp_xi = sb.Gradient(v);
  Tensor<DataMesh> xi(2,"1",DataMesh::Empty);
  xi(0) = tmp_xi(1);
  xi(1) = -tmp_xi(0);

  //create storage
  DataMesh scaleFactor(Mesh(theta.Extents()));

  for(int i=0; i<scaleFactor.Size(); i++){
    //std::cout << POSITION << " Testing " << i 
    //          << " of " << scaleFactor.Size() << std::endl; //can be deleted
    scaleFactor[i] = NormalizeAKVAtOnePoint(sb,Psi,xi,rad,theta[i],phi[i],printSteps)-1.;
  }

  scaleFactor *= scaleFactor;
  const DataMesh r2p4 = rad*rad*Psi*Psi*Psi*Psi;
  const double area = sb.ComputeCoefficients(r2p4)[0];

  return sqrt(sb.ComputeCoefficients(scaleFactor*r2p4)[0])/area;
}

//determines the path length of following the approximate
//Killing vector around the sphere starting at (theta, phi)
//and returns the ratio of that path to the expected
//value of 2*Pi (the normalization scale factor)
double NormalizeAKVAtOnePoint(const SurfaceBasis& sb,
                              const DataMesh& rotated_Psi,
                              const DataMesh& rotated_v,
                              const double& rad,
                              const double& thetap, /*=M_PI/2.0*/
                              const double& phip, /*=0.0*/
                              const bool& printSteps)
{
  //create xi
  Tensor<DataMesh> tmp_xi = sb.Gradient(rotated_v);
  Tensor<DataMesh> xi(2,"1",DataMesh::Empty);
  xi(0) = tmp_xi(1);
  xi(1) = -tmp_xi(0);

  return NormalizeAKVAtOnePoint(sb, rotated_Psi, xi, rad, thetap, phip, printSteps);
}

double NormalizeAKVAtOnePoint(const SurfaceBasis& sb,
                              const DataMesh& rotated_Psi,
                              const Tensor<DataMesh>& xi,
                              const double& rad,
                              const double& thetap, /*=M_PI/2.0*/
                              const double& phip, /*=0.0*/
                              const bool& printSteps)
{
  //Rescale the Killing vector. For a rotational Killing vector,
  //the affine path length should be 2*Pi.  Also, test various
  //paths to ensure we have an actual Killing field
  double t; //affine path length

  //bool goodtheta = KillingPath(sb, rotated_Psi, xi, rad, t, thetap, phip, printSteps);
  double thetaOffAxis = 0.; double phiOffAxis = 0.;
  bool goodtheta = KillingPath(sb, rotated_Psi, xi, rad, t, thetap, phip,
                               thetaOffAxis, phiOffAxis, printSteps);
  REQUIRE(goodtheta, "Killing trajectory did not close " << POSITION);

  const double scale = t/(2.0*M_PI);

  return scale;
} //end normalizeKillingVector

//This routine will print the value of NormalizeAKVAtAllPoints for 
//a range of scale factors in the neighborhood of scaleFactor1-4.  This function
//assumes that there is an optimal (minimal) scale factor value in this range,
//and uses a modified bisection method to find it within a given tolerance
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
                      const bool& printBisectionResults)
{
  //quick routine to determine the minimum and maximum values of the scale
  //factors that were passed in
  double min = scaleFactor1 < scaleFactor2 ? scaleFactor1 : scaleFactor2;
         min = scaleFactor3 < min ? scaleFactor3 : min;
         min = scaleFactor4 < min ? scaleFactor4 : min;

  double max = scaleFactor1 > scaleFactor2 ? scaleFactor1 : scaleFactor2;
         max = scaleFactor3 > min ? scaleFactor3 : min;
         max = scaleFactor4 > min ? scaleFactor4 : min;
std::cout << POSITION << "temporary edit to get optimization working" << std::endl;
min = scaleFactor1;
max = scaleFactor1;

  //start just outside these values
  min *= 0.9;
  max *= 1.1;

  //evaluate the function at these two scale factors
  double lowSurfaceValue = 
           NormalizeAKVAtAllPoints(sb, rotated_Psi, theta, phi, rotated_v*min, rad, printSteps);
  if(printBisectionResults)
    std::cout << std::setprecision(15) << min << " " << lowSurfaceValue << std::endl;
  double highSurfaceValue = 
           NormalizeAKVAtAllPoints(sb, rotated_Psi, theta, phi, rotated_v*max, rad, printSteps);
  if(printBisectionResults)
    std::cout << std::setprecision(15) << max << " " << highSurfaceValue << std::endl;

  //do a bisection to find the minimum surface value
  //this assumes that the function is continuous and symmetric about the minimum
  //stops when the difference in SurfaceValue between min/max is less than 1e-04%
  while(fabs(lowSurfaceValue-highSurfaceValue)/lowSurfaceValue > 1.e-06){
    const double difference = max-min;
    if(lowSurfaceValue < highSurfaceValue){
      max -= difference/4.;
      highSurfaceValue = NormalizeAKVAtAllPoints(sb,rotated_Psi,theta,phi,rotated_v*max,rad,printSteps);
      if(printBisectionResults)
        std::cout << std::setprecision(15) << max << " " << highSurfaceValue << std::endl;
    } else {
      min += difference/4.;
      lowSurfaceValue = NormalizeAKVAtAllPoints(sb,rotated_Psi,theta,phi,rotated_v*min,rad,printSteps);
      if(printBisectionResults)
        std::cout << std::setprecision(15) << min << " " << lowSurfaceValue << std::endl;
    }
  }

  return (highSurfaceValue < lowSurfaceValue ? max : min);
}

//integrates along a particular Killing path on the surface
bool KillingPath(const SurfaceBasis& sb,
                 const DataMesh& Psi,
                 const Tensor<DataMesh>& xi,
                 const double& rad,
                 double& t,
                 const double& theta,
                 const double& phi,
                 double& thetaOffAxis,
                 double& phiOffAxis,
                 const bool printSteps /*=false*/)
{
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
     
  t = 0.0; //path length
  double t1 = 1.e10;
  double h = 1.e-6; //step size
  double y[2] = {theta, phi-2.*M_PI}; //this makes the stopping criteria easier
  bool limit_h = false;
  double hmax = h; //maximum step size
  //const double hmax2 = 2.0*M_PI/100.0; //hard maximum for h size, original
  const double hmax2 = 0.1; //hard maximum for h size, new

  //useful for calculating the "center point" of the path taken
  double pathTotal[2] = {0.,0.};
  double pathAvg[2] = {0.,0.};


/*
//keep this following segment for final production
  if(printSteps){
    std::cout << "START: y = ( " 
	      << std::setprecision(8) << std::setw(10) << y[0] << " , " 
	      << std::setprecision(8) << std::setw(10) << y[1]+2.*M_PI << " ); h = " 
	      << std::setprecision(8) << std::setw(10) << h 
	      << std::endl; 
  }
*/     
  int iter = 0;
  const int iterMax = 1000;
  while (true && iter<iterMax) {
    iter++;
    const double ysave[2] = {y[0],y[1]};

    const double tsave = t;
    const int status = gsl_odeiv2_evolve_apply (e, c, s,
						&sys, 
						&t, t1,
						&h, y);
/*
    if(iter == 1){
      std::cout << "y = ( "
  	        << std::setprecision(8) << std::setw(10) << y[0] << " , " 
  	        << std::setprecision(8) << std::setw(10) << y[1]+2.*M_PI << " )" 
                << std::endl;
    }
*/
/*
    //if(printSteps){
    //if(printSteps && y[1]+2.*M_PI<.0){
    //if(iter > 0.9*iterMax){
    if(h < 1.e-7){
      std::cout << "iter = " << iter << " t = " 
    	        << std::setprecision(8) << std::setw(10) << t 
  	        << " : y = ( " 
  	        << std::setprecision(8) << std::setw(10) << y[0] << " , " 
  	        << std::setprecision(8) << std::setw(10) << y[1]+2.*M_PI << " ); h = " 
  	        << std::setprecision(8) << std::setw(10) << h 
  	        << std::endl;
    }
*/
    if(printSteps){
      //note that this gives output appropriate for gnuplot's splot function
      std::cout << std::setprecision(8) << std::setw(10) <<  y[1]+2.*M_PI
      << " " << std::setprecision(8) << std::setw(10) << -(y[0]-M_PI/2.) << std::endl;
    }
/*
    if(y[1]>0 && iter%100==0){
      std::cout << "iter = " << iter << " t = " 
    	        << std::setprecision(8) << std::setw(10) << t 
  	        << " : y = ( " 
  	        << std::setprecision(8) << std::setw(10) << y[0] << " , " 
  	        << std::setprecision(8) << std::setw(10) << y[1]+2.*M_PI << " ); h = " 
  	        << std::setprecision(8) << std::setw(10) << h 
  	        << std::endl;
    }
*/
    ASSERT(status==GSL_SUCCESS,"Path Integration failed");

    if(limit_h && h > hmax) h = hmax;
    if(fabs(y[1] - phi) < 1.e-12) break; //original condition assuming navigation about the pole
    //if a path brings you back to the original point,
    //you have closure, and that's what we're interested in
    if(   fabs(y[1]+2.*M_PI - phi) < 1.e-6      //close to original phi
       && t > 1.e-2 ){
       //&& printSteps){
      std::cout << POSITION << " Off-axis closure rules have been triggered." << std::endl;
      std::cout << "Path average (theta, phi) = (" << pathTotal[0]
                << ", " << pathTotal[1]/(2.*M_PI)-M_PI << ")" << std::endl;
      break;
    }
    else if(y[1] > phi) { //if solver went too far...
      if(!limit_h) hmax = h;
      limit_h = true;
      h = hmax *= 0.5; //...reset h, hmax
      y[0] = ysave[0]; y[1] = ysave[1]; t = tsave; //return variables to previous state
      gsl_odeiv2_evolve_reset(e); //return solver to previous state
    } //end ifs
    //alternative else if
    else if(ysave[1] < phi-2.*M_PI && y[1]>phi-2.*M_PI){ //Killing path has been followed too far
      //std::cout << POSITION << " Off-axis overstepping rules have been triggered." << std::endl;
      if(!limit_h) hmax = h;
      limit_h = true;
      h = hmax *= 0.5; //...reset h, hmax
      y[0] = ysave[0]; y[1] = ysave[1]; t = tsave; iter--; //return variables to previous state
      gsl_odeiv2_evolve_reset(e); //return solver to previous state
    }
    if(h > hmax2) h = hmax2;

    //diagnostics on average (theta, phi) position.
    //if average theta != input theta, then the axis must be off
    //if average phi != M_PI, then the axis must be off
    pathTotal[0] += (y[0]-theta) * fabs(y[0]-ysave[0]);
    pathTotal[1] += (y[1]+2.*M_PI-phi) * fabs(y[1]-ysave[1]);

  } //end while

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);

  const bool closedPath = (fabs(y[0] - theta) < 1.e-6);
  if(!closedPath){
    std::cout << "##> Theta diff: " << std::setprecision(6) << y[0] - theta << std::endl;
    std::cout << "Iterations: " << iter << std::endl;
  }

    //pathAvg[0] = pathTotal[0]/1.; pathAvg[1] = pathTotal[1]/1.;
    //std::cout << "Path average 2 (th,ph) = (" << pathAvg[0] << ", " << pathAvg[1] << ")" << std::endl;
    pathAvg[0] = pathTotal[0]/1.; pathAvg[1] = pathTotal[1]/(2.*M_PI)-M_PI;
    //std::cout << "Path average 3 (th,ph) = (" << pathAvg[0] << ", " 
    //          << pathAvg[1] << ") Starting at = (" << theta << ", " << phi << ")" << std::endl;
    thetaOffAxis = pathAvg[0];
    phiOffAxis = pathAvg[1];
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

//this function returns the value of
//   (1./(sqrt(2.)*Pi) Integral ( 0.5 * (Ricci) * D^i (v_1) * D_i (v_2) ) dA
// = (1./(sqrt(2.)*Pi) Integral ( 0.5 * (Ricci) * (1/r2p4) * gradient(v_1)*gradient(v_2) ) * r2p4 dOmega
// = (1./(sqrt(2.)*Pi) Integral ( 0.5 * (Ricci) * gradient(v_1)*gradient(v_2) ) dOmega
double AKVInnerProduct(const DataMesh& v1,
                       const DataMesh& v2,
                       const DataMesh& Ricci,
                       const SurfaceBasis& sb,
                       const bool& withRicciScaling)
{
  const Tensor<DataMesh> Gradv1 = sb.Gradient(v1);
  const Tensor<DataMesh> Gradv2 = sb.Gradient(v2);
  const DataMesh integrand = 0.5*(Gradv1(0)*Gradv2(0)+Gradv1(1)*Gradv2(1));
  if(withRicciScaling){
    return sb.ComputeCoefficients(integrand*Ricci)[0];
  } else {
    std::cout << POSITION << " AKVInnerProduct w/o Ricci" << std::endl;
    return sb.ComputeCoefficients(integrand)[0];
  }
}

//returns the scaling factors related to various forms of the AKVInnerProduct
MyVector<double> InnerProductScaleFactors(const DataMesh& v1,
                       const DataMesh& v2,
                       const DataMesh& Ricci,
                       const DataMesh& r2p4,
                       const SurfaceBasis& sb,
                       const bool& withRicciScaling)
{
  //ip1: Integral ( 0.5 * Ricci * gradient(v_1)*gradient(v_2) ) / r2p4 dOmega = (8.*Pi/3.) s
  //in the Ricci argument, r2p4 compensates for the integrand division by r2p4
  const double s_ip1 = (8./(3.*sqrt(2.))) / AKVInnerProduct(v1,v2,Ricci/r2p4,sb,withRicciScaling);

  const double tmp_ip = AKVInnerProduct(v1,v2,Ricci,sb,withRicciScaling);
  const double area = SurfaceArea(r2p4,sb);

  //ip3: Integral ( 0.5 * Ricci * gradient(v_1)*gradient(v_2) ) dOmega = (2/3)A s
  const double s_ip3 = (2./3.)*area / tmp_ip;

  //ip2: Integral ( 0.5 * Ricci * gradient(v_1)*gradient(v_2) ) dOmega = (2/3)A
  //since s=1, it looks like s_ip2 is just the sqrt of the scale factor for ip3
  const double s_ip2 = sqrt(s_ip3);

  return MyVector<double>(MV::fill, s_ip1, s_ip2, s_ip3) ;
}

//this function returns the proper area of the surface
//   (1./(sqrt(2.)*Pi) Integral 1. dA
// = (1./(sqrt(2.)*Pi) Integral r2p4 dOmega
double SurfaceArea(const DataMesh& r2p4,
                   const SurfaceBasis& sb)
{
  return sb.ComputeCoefficients(r2p4)[0];
}

//Gram Schmidt orthogonalization
//Takes the DataMesh that should not change (fixedMesh), the inner product of that
//DataMesh with itself (fixedInnerProduct), the DataMesh that needs to be orthogonalized
//(flexibleMesh), and the inner product of those two meshes (crossTermInnerProduct)
void GramSchmidtOrthogonalization(const DataMesh& fixedMesh,
                                  const double fixedInnerProduct,
                                  DataMesh& flexibleMesh,
                                  const double crossTermInnerProduct)
{
  flexibleMesh -= (crossTermInnerProduct/fixedInnerProduct) * fixedMesh;
}

//this function will create an initial guess for the next axis of
//symmetry based on previous solutions.  It requires theta, phi for
//prior axes solutions, and an index which indicates whether this is
//the first, second, or third guess
void AxisInitialGuess(double theta[], double phi[], const int index)
{
  switch(index){
    case 0:
      std::cout << "symmetry axis guess : " 
                << theta[0]*180./M_PI << " " << phi[0]*180./M_PI << std::endl;
      break;
    case 1:
      theta[1] = M_PI/2.;
      phi[1] = atan2(sin(phi[0]), -cos(phi[0]));
      std::cout << "symmetry axis guess : "
                << theta[1]*180./M_PI << " " << phi[1]*180./M_PI << std::endl;
      break;
    case 2:
      //perform the cross product of the previous two solutions
      const double alpha = sin(theta[0])*sin(phi[0])*cos(theta[1])
                        -sin(theta[1])*sin(phi[1])*cos(theta[0]);
      const double beta = cos(theta[0])*sin(theta[1])*cos(phi[1])
                        -cos(theta[1])*sin(theta[0])*cos(phi[0]);
      const double gamma = sin(theta[0])*cos(phi[0])*sin(theta[1])*sin(phi[1])
                        -sin(theta[1])*cos(phi[1])*sin(theta[0])*sin(phi[0]);
      theta[2] = atan2(sqrt(alpha*alpha+beta*beta),gamma);
      phi[2] = atan2(beta, alpha);
      std::cout << "symmetry axis guess : "
                << theta[2]*180./M_PI << " " << phi[2]*180./M_PI << std::endl;
      break;
  }
}

//performs diagnostics on the approximate Killing vector solution
void KillingDiagnostics(const SurfaceBasis& sb,
                        const DataMesh& L,
                        const DataMesh& Psi,
                        const Tensor<DataMesh>& xi,
                        const double& rad,
                        const MyVector<bool>& printDiagnostic)
{
  const DataMesh& rp2 = rad*Psi*Psi;
  const DataMesh& r2p2 = rad*rp2;
  const DataMesh& r2p4 = rp2*rp2;
  const double r2p4_00 = sb.ComputeCoefficients(r2p4)[0];//(0,0) component of coefficients

  DataMesh div = sb.Divergence(xi);

  //-----Divergence-------
  if(printDiagnostic[0]){
    DataMesh tmp_div = div/r2p2;
    std::cout << "L2 Norm of Divergence = "
              << std::setprecision(12)
              //<< sb.ComputeCoefficients(div*div)[0])/r2p4_00 << std::endl;
              << sqrt(sqrt(2.)*sb.ComputeCoefficients(tmp_div*tmp_div)[0]/4.) << std::endl;
  }

  if(printDiagnostic[1] || printDiagnostic[2]){
    DataMesh vort = sb.Vorticity(xi);
    DataMesh vortM2L = vort/r2p2 - 2.0*Psi*Psi*L;

    //-----Vorticity-------
    if(printDiagnostic[1]){
      std::cout << "L2 Norm of Vorticity = "
                << std::setprecision(12)
                << sqrt(sqrt(2.)*sb.ComputeCoefficients(vortM2L*vortM2L)[0]/4.) << std::endl;
    }

    //-----S_{ij}S^{ij}-------
    if(printDiagnostic[2]){
      Tensor<DataMesh> gradlncf = sb.Gradient(log(Psi));
      Tensor<DataMesh> Dxtheta = sb.VectorColatitudeDerivative(xi);
      Tensor<DataMesh> Dxphi(2,"1",DataMesh::Empty);

      Dxtheta(0) -= 2.0*( xi(0)*gradlncf(0) - xi(1)*gradlncf(1) );
      Dxtheta(1) -= 2.0*( xi(0)*gradlncf(1) + xi(1)*gradlncf(0) );
      Dxphi(0) = Dxtheta(1) - vort;
      Dxphi(1) = div - Dxtheta(0);

      DataMesh SS =  Dxtheta(0)*Dxtheta(0)
                   + Dxtheta(1)*Dxtheta(1)
                   + Dxphi(0)*Dxphi(0)
                   + Dxphi(1)*Dxphi(1);

      //these factors of r2p4 are setup to do a dA integral
      SS /= r2p4;
      SS -= 2.0*r2p4*L*L;

      const double intSSdA = sb.ComputeCoefficients(SS)[0]*M_PI*sqrt(2.);
      std::cout << "Integral S_{ij}S^{ij} dA " << intSSdA << std::endl;
    }
  }

  if(printDiagnostic[3] || printDiagnostic[4] || printDiagnostic[5]){
    Tensor<DataMesh> GradL = sb.Gradient(L);
    DataMesh xiGradL = xi(0)*GradL(0) + xi(1)*GradL(1);
    const DataMesh& llncf = sb.ScalarLaplacian(log(Psi));
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
                << sb.ComputeCoefficients(f_L*f_L)[0]/r2p4_00 << std::endl;
    } //end f_L

    //-----Norm of f_Lambda-------
    if(printDiagnostic[5]){
      const Tensor<DataMesh> GradRm1 = sb.Gradient(Rm1);
      DataMesh f_lam = GradRm1(0)*GradL(1) - GradRm1(1)*GradL(0);
      f_lam /= rp2;

      std::cout << "L2 Norm of f_lam= "
                << std::setprecision(12)
                << sb.ComputeCoefficients(f_lam*f_lam)[0]/r2p4_00 << std::endl;
    } // end f_Lambda

    //-----Norm of xi*GradL-------
    if(printDiagnostic[3]){
      xiGradL = xi(0)*GradL(0) + xi(1)*GradL(1);
      xiGradL /= rp2;

      std::cout << "L2 Norm of xi*Grad(L) = "
                << std::setprecision(12)
                << sb.ComputeCoefficients(xiGradL*xiGradL)[0]/r2p4_00 << std::endl;
    } //end xi*DivL

  } //end print conditionals
}//end KillingDiagnostics

