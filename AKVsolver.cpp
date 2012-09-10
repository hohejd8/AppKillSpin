#include "AKVsolver.hpp"
#include "Utils/StringParsing/OptionParser.hpp"
#include "Utils/LowLevelUtils/Position.hpp"
#include "Spectral/BasisFunctions/SpherePackIterator.hpp"
#include "gsl/gsl_odeiv2.h"
#include <iomanip>

/*
HELPFUL DIAGNOSTIC TOOLS WHEN COMPARING TO AppKillSpin:

  //print v (or other DataMesh)
  std::cout << "v " << POSITION << std::endl;
  for(int i=0; i<mNth; ++i){
      for(int j=0; j<mNph; ++j) {
            std::cout << std::setprecision(10) << v[i*mNph+j] << " " ;
      }
      std::cout << std::endl;
  }
  std::cout << "\n" << std::endl;

*/

/*
TO DO:
eliminate small MyVectors
put in 1D root finding, try testing (1,0), 
*/
//-----------------------------------------

//function header for use with gsl root finders
int AKVsolver(const gsl_vector * x,
              void *params,
              gsl_vector * f){

  const double THETA = gsl_vector_get (x,0);
  const double thetap = gsl_vector_get (x,1);
  const double phip = gsl_vector_get (x,2);

  //for use with gsl root finder
  const DataMesh& theta = static_cast<struct rparams*>(params)->theta;
  const DataMesh& phi = static_cast<struct rparams*>(params)->phi;
  const double& rad = static_cast<struct rparams*>(params)->rad;
  const SurfaceBasis& sb = static_cast<struct rparams*>(params)->sb;
  const DataMesh& Psi = static_cast<struct rparams*>(params)->Psi;
  DataMesh& L = static_cast<struct rparams*>(params)->L;
  DataMesh& v = static_cast<struct rparams*>(params)->v;
  const double L_resid_tol = static_cast<struct rparams*>(params)->L_resid_tol;
  const double v_resid_tol = static_cast<struct rparams*>(params)->v_resid_tol;
  const bool printResiduals = static_cast<struct rparams*>(params)->printResiduals;

  SpherePackIterator sit(theta.Extents()[0],theta.Extents()[1]);

  //eq. 78, 93
  L = cos(thetap)*cos(theta)
        + sin(thetap)*sin(theta)*(cos(phip)*cos(phi) + sin(phip)*sin(phi) );

  //eq. 95
  v = L*rad*rad;

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

  //^2R, eq. 20
  const DataMesh& llncf = sb.ScalarLaplacian(log(Psi)); //original

  const DataMesh& hR = (1.0-2.0*llncf) / (Psi*Psi*Psi*Psi*rad*rad);

  const Tensor<DataMesh>& hGradR = sb.Gradient(hR);

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
    RHS = -RHS + (4.0*llncf*L + hGradR(0)*Gradv(0) + hGradR(1)*Gradv(1) -2.0*L)*(1.0-THETA);

    //perform harmonic analysis on RHS
    RHS_ha = sb.ComputeCoefficients(RHS);

    //keep track of l=1 values
    ic10 = RHS_ha[sit(1,0,SpherePackIterator::a)];
    ic1p = RHS_ha[sit(1,1,SpherePackIterator::a)];
    ic1m = RHS_ha[sit(1,1,SpherePackIterator::b)];

    //print residuals
    if(printResiduals){
      std::cout << "ic10 = " << ic10
                << " ic1p = " << ic1p
                << " ic1m = " << ic1m << std::endl;
    }

    //remove the l=1 modes
    RHS_ha[sit(1,0,SpherePackIterator::a)] = 0.0;
    RHS_ha[sit(1,1,SpherePackIterator::a)] = 0.0;
    RHS_ha[sit(1,1,SpherePackIterator::b)] = 0.0;

    //recompute RHS from 'fixed' analysis
    RHS = sb.Evaluate(RHS_ha);

    //compute RHS^2
    RHS *= RHS;

    //evaluate normalization factor
    const double norm_RHS_L = sqrt(sqrt(2.0)*sb.ComputeCoefficients(RHS)[0]/4.0);

    //invert (Laplacian + 2)
    for(sit.Reset(); sit; ++sit){
      if(sit.l()==1){ //must zero out n=1 modes
        RHS_ha[sit()] = 0.0;
      } else {
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
    const DataMesh twocf4r2 = 2.0*Psi*Psi*Psi*Psi*rad*rad;
    RHS = -RHS - twocf4r2*L;

    //remove the l=0 mode from RHS
    RHS_ha = sb.ComputeCoefficients(RHS);
    RHS_ha[sit(0,0,SpherePackIterator::a)] = 0.0;
    RHS = sb.Evaluate(RHS_ha);

    //compute RHS^2
    RHS *= RHS;

    //compute norm
    const double norm_RHS_v = sqrt(sqrt(2.0)*sb.ComputeCoefficients(RHS)[0]/4.0);

    //invert (Laplacian + 0)
    for(sit.Reset(); sit; ++sit){
      if(sit.l()==0){
        RHS_ha[sit()] = 0.0;
      } else {
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


//takes a DataMesh on the sphere and rotates it by some amount Theta, Phi
DataMesh RotateOnSphere
         (const DataMesh& collocationvalues,
          const DataMesh& thetaGrid,
          const DataMesh& phiGrid,
          const SurfaceBasis& sb,
          const double Theta,
          const double Phi) 
{
  DataMesh result(DataMesh::Empty);

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

double normalizeKillingVector(const SurfaceBasis& sb,
                              const DataMesh& Psi,
                              DataMesh& v,
                              const double& rad)
{
  //rotate v
//-:  DataMesh rotated_v = RotateOnSphere(v,
//-:                       skwm.Grid().SurfaceCoords()(0),
//-:                       skwm.Grid().SurfaceCoords()(1),
//-:                       sb,
//-:                       thetap,
//-:                       phip);
  DataMesh rotated_v = v;
  //rotate Psi
//-:  DataMesh rotated_Psi = RotateOnSphere(Psi,
//-:                       skwm.Grid().SurfaceCoords()(0),
//-:                       skwm.Grid().SurfaceCoords()(1),
//-:                       sb,
//-:                       thetap,
//-:                       phip);
  DataMesh rotated_Psi = Psi;

  //create xi
  Tensor<DataMesh> tmp_xi = sb.Gradient(rotated_v);//rotated v
  Tensor<DataMesh> xi(2,"1",DataMesh::Empty);
  xi(0) = tmp_xi(1);
  xi(1) = -tmp_xi(0);

  //Rescale the Killing vector. For a rotational Killing vector,
  //the affine path length should be 2*Pi.  Also, test various
  //paths to ensure we have an actual Killing field
  double t; //affine path length

  bool goodtheta = KillingPath(sb, rotated_Psi, xi, rad, t, M_PI/2.0);
  REQUIRE(goodtheta, "Killing trajectory did not close " << POSITION);
  const double scale = t/(2.0*M_PI);
  std::cout << "Theta = " << std::setprecision(8) << std::setw(10)
            << M_PI/2.0 << " : T = "
            << std::setprecision(10) << t
            << " (" << scale << ")" << std::endl;

  KillingPath(sb, rotated_Psi, xi, rad, t, 0.5*M_PI/2.0);
  std::cout << "Theta = " << std::setprecision(8) << std::setw(10)
            << 0.5*M_PI/2.0 << " : T = "
            << std::setprecision(10) << t
            << " (" << t/(2.0*M_PI*scale) << ")" << std::endl;

  KillingPath(sb, rotated_Psi, xi, rad, t, 0.25*M_PI/2.0);
  std::cout << "Theta = " << std::setprecision(8) << std::setw(10)
            << 0.25*M_PI/2.0 << " : T = "
            << std::setprecision(10) << t
            << " (" << t/(2.0*M_PI*scale) << ")" << std::endl;

  KillingPath(sb, rotated_Psi, xi, rad, t, 0.125*M_PI/2.0);
  std::cout << "Theta = " << std::setprecision(8) << std::setw(10)
            << 0.125*M_PI/2.0 << " : T = "
            << std::setprecision(10) << t
            << " (" << t/(2.0*M_PI*scale) << ")" << std::endl;

  return scale;
} //end normalizeKillingVector

bool KillingPath(const SurfaceBasis& sb,
                    const DataMesh& Psi,
                    const Tensor<DataMesh>& xi,
                    const double& rad,
                    double& t,
                    const double theta)
{
  bool printSteps = false;
  //perform harmonic analysis on Psi, xi to save from repetitive computation
  //to save from repetitive computation in PathDerivs
  const DataMesh Psi_ha = sb.Evaluate(Psi);
  //const Tensor<DataMesh> xi_ha = sb.EvaluateVector(xi);
  //struct ODEparams params = {sb, Psi_ha, xi_ha, rad};
  struct ODEparams params = {sb, Psi_ha, xi, rad};

  const gsl_odeiv2_step_type *const T = gsl_odeiv2_step_rkck;
  gsl_odeiv2_step *const s = gsl_odeiv2_step_alloc (T, 2);
  gsl_odeiv2_control *const c = gsl_odeiv2_control_y_new (1e-12, 1.e-12);
  gsl_odeiv2_evolve *const e = gsl_odeiv2_evolve_alloc (2);
  gsl_odeiv2_system sys = {PathDerivs, NULL, 2, &params};
     
  t = 0.0;
  double t1 = 1.e10;
  double h = 1e-6;
  double y[2] = { theta, 0.0 };
  bool limit_h = false;
  double hmax = h;
  const double hmax2 = 2.0*M_PI/100.0;

  if(printSteps){
    std::cout << "START: y = ( " 
	      << std::setprecision(8) << std::setw(10) << y[0] << " , " 
	      << std::setprecision(8) << std::setw(10) << y[1] << " ); h = " 
	      << std::setprecision(8) << std::setw(10) << h 
	      << std::endl; 
  }
     
  while (true) {
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
  	        << std::setprecision(8) << std::setw(10) << y[1] << " ); h = " 
  	        << std::setprecision(8) << std::setw(10) << h 
  	        << std::endl;
    }

    ASSERT(status==GSL_SUCCESS,"Path Integration failed");
    if(limit_h && h > hmax) h = hmax;
    if(fabs(y[1] - 2.*M_PI) < 1.e-10) break;
    else if(y[1] > 2.*M_PI) { //if solver went too far...
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

int PathDerivs(double t, const double y[], double f[], void *params)
{
  const SurfaceBasis& sb = static_cast<struct ODEparams*>(params)->sb;
  const DataMesh& Psi_ha = static_cast<struct ODEparams*>(params)->Psi_ha;
  const Tensor<DataMesh>& xi_ha = static_cast<struct ODEparams*>(params)->xi_ha;
  const double& rad = static_cast<struct ODEparams*>(params)->rad;
//std::cout << POSITION << std::endl;
  double psiAtPoint = sb.EvaluateFromCoefficients(Psi_ha, y[0], y[1]);
  //double psiAtPoint = sb.Evaluate(Psi_ha, y[0], y[1]);
std::cout << psiAtPoint << " " << POSITION << std::endl;
  //MyVector<double> result = sb.EvaluateVectorFromCoefficients(xi_ha, y[0], y[1]);
  MyVector<double> result = sb.EvaluateVector(xi_ha, y[0], y[1]);
std::cout << result << " " << POSITION << std::endl;
  const double norm = 1.0 / (psiAtPoint*psiAtPoint*psiAtPoint*psiAtPoint
                             *rad*rad);
  f[0] = result[0]*norm;
  f[1] = result[1]*norm/sin(y[0]);

  return GSL_SUCCESS;
}

void KillingDiagnostics(const SurfaceBasis& sb,
                        const DataMesh& L,
                        const DataMesh& Psi,
                        const Tensor<DataMesh>& xi,
                        const double& rad,
                        const MyVector<bool>& printDiagnostic)
{
  const DataMesh& p2r2 = rad*rad*Psi*Psi;

  DataMesh div = sb.Divergence(xi) / p2r2;

  //-----Divergence-------
  if(printDiagnostic[0]){
    const DataMesh& div2norm = sb.ComputeCoefficients(div*div);
    std::cout << "L2 Norm of Divergence = "
              << std::setprecision(12) << sqrt(sqrt(2.)*div2norm[0]/4.) << std::endl;
  }

  if(printDiagnostic[1] || printDiagnostic[2]){
    DataMesh vort = sb.Vorticity(xi) / p2r2 - 2.0*Psi*Psi*L;

    //-----Vorticity-------
    if(printDiagnostic[1]){
      const DataMesh vort2norm = sb.ComputeCoefficients(vort*vort);
      std::cout << "L2 Norm of Vorticity = "
                << std::setprecision(12) << sqrt(sqrt(2.)*vort2norm[0]/4.) << std::endl;
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
    const DataMesh Rm1 = Psi*Psi*Psi*Psi / (1.-2.*llncf);
    GradL(0) *= Rm1;
    GradL(1) *= Rm1;

    //-----Norm of f_L-------
    if(printDiagnostic[3]){
      DataMesh f_L = sb.Divergence(GradL);
      f_L *= 0.5/(Psi*Psi);
      f_L += Psi*Psi*L;

      const DataMesh fL_ha = sb.ComputeCoefficients(f_L*f_L);

      std::cout << "L2 Norm of f_L= "
                << std::setprecision(12) << sqrt(sqrt(2.)*fL_ha[0]/4.) << std::endl;
    } //end f_L

    //-----Norm of f_Lambda-------
    if(printDiagnostic[4]){
      const Tensor<DataMesh> GradRm1 = sb.Gradient(Rm1);
      DataMesh f_lam = GradRm1(0)*GradL(1) - GradRm1(1)*GradL(0);
      f_lam *= 0.5/(Psi*Psi);

      const DataMesh flam_ha = sb.ComputeCoefficients(f_lam*f_lam);

      std::cout << "L2 Norm of f_lam= "
                << std::setprecision(12) << sqrt(sqrt(2.)*flam_ha[0]/4.) << std::endl;
    } // end f_Lambda

    //-----Norm of xi*DivL-------
    if(printDiagnostic[5]){
      xiGradL = xi(0)*GradL(0) + xi(1)*GradL(1);
      xiGradL /= (rad*rad*Psi*Psi);

      const DataMesh xiGradL_ha = sb.ComputeCoefficients(xiGradL*xiGradL);

      std::cout << "L2 Norm of xi*Div(L) = "
                << std::setprecision(12) << sqrt(sqrt(2.)*xiGradL_ha[0]/4.) << std::endl;
    } //end xi*DivL

  } //end print conditionals
}//end KillingDiagnostics







