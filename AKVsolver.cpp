#include "AKVsolver.hpp"
#include "Utils/StringParsing/OptionParser.hpp"
#include "Utils/LowLevelUtils/Position.hpp"
#include "Spectral/BasisFunctions/SpherePackIterator.hpp"
#include "gsl/gsl_odeiv2.h"
#include <iomanip>

//#include "Utils/DataMesh/DataMeshNorms.hpp"

//-----------------------------------------
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

  //print xi (or other Tensor<DataMesh>)
  std::cout << "xi(0) " << POSITION << std::endl;
  for(int i=0; i<mNth; ++i){
      for(int j=0; j<mNph; ++j) {
            std::cout << std::setprecision(10) << xi(0)[i*mNph+j] << " " ;
      }
      std::cout << std::endl;
  }
  std::cout << "\n" << std::endl;
  std::cout << "xi(1) " << POSITION << std::endl;
  for(int i=0; i<mNth; ++i){
      for(int j=0; j<mNph; ++j) {
            std::cout << std::setprecision(10) << xi(1)[i*mNph+j] << " " ;
      }
      std::cout << std::endl;
  }
  std::cout << "\n" << std::endl;

*/
//-----------------------------------------

//-----------------------------------------
/*
TO DO:
eliminate small MyVectors
remote unnecessary diagnostic printing
fix tensor<datamesh> constructors
put in 1D root finding, try testing (1,0), 
eliminate skwm
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
  const StrahlkorperWithMesh& skwm = static_cast<struct rparams*>(params)->skwm;
  //const DataMesh& theta = static_cast<struct rparams*>(params)->theta;
  //const DataMesh& phi = static_cast<struct rparams*>(params)->phi;
  //const double& rad = static_cast<struct rparams*>(params)->rad;
  const DataMesh& Psi = static_cast<struct rparams*>(params)->Psi;
  DataMesh& L = static_cast<struct rparams*>(params)->L;
  DataMesh& v = static_cast<struct rparams*>(params)->v;
  const double L_resid_tol = static_cast<struct rparams*>(params)->L_resid_tol;
  const double v_resid_tol = static_cast<struct rparams*>(params)->v_resid_tol;
  const bool verbose = static_cast<struct rparams*>(params)->PrintResiduals;

  const SurfaceBasis sb(skwm.Grid());

  const DataMesh& theta = skwm.Grid().SurfaceCoords()(0);
  const DataMesh& phi   = skwm.Grid().SurfaceCoords()(1);
  SpherePackIterator sit(theta.Extents()[0],theta.Extents()[1]);

  const DataMesh& rad   = skwm.Radius();

  //eq. 78, 93
  L = cos(thetap)*cos(theta)
        + sin(thetap)*sin(theta)*(cos(phip)*cos(phi) + sin(phip)*sin(phi) );

  //eq. 95
  v = L*rad*rad;

  //eq. 94, make sure only l=1 mode exists by setting all l!=1 modes to zero
  DataMesh L_ha = sb.ComputeCoefficients(L);
  DataMesh v_ha = sb.ComputeCoefficients(v);

  //replace with scalar iterator tools
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
  Tensor<DataMesh> Gradv = hGradR; //dummy initialization. Gradient of v, eq. 97

  while(unsolved){
    //eq. 97, first term: compute Laplacian of L
    RHS = sb.ScalarLaplacian(L_ha); //RHS is in collocation terms

    //eq. 97, component of third term: compute gradient of v
    Gradv = sb.Gradient(v_ha); //Gradv is in collocation terms

    //compute eq. 97
    RHS = -RHS + (4.0*llncf*L + hGradR(0)*Gradv(0) + hGradR(1)*Gradv(1) -2.0*L)*(1.0-THETA);

    //perform harmonic analysis on RHS
    RHS_ha = sb.ComputeCoefficients(RHS);

    //keep track of l=1 values
    ic10 = RHS_ha[sit(1,0,SpherePackIterator::a)];
    ic1p = RHS_ha[sit(1,1,SpherePackIterator::a)];
    ic1m = RHS_ha[sit(1,1,SpherePackIterator::b)];

    if(verbose){
      //std::cout << "ic10 = " << ic10 << " "
      //    << "ic1p = " << ic1p << " "
      //    << "ic1m = " << ic1m << std::endl;
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
    RHS = sb.ScalarLaplacian(v_ha); //collocation
    const DataMesh twocf4r2 = 2.0*Psi*Psi*Psi*Psi*rad*rad;
    RHS = -RHS - twocf4r2*L; //collocation

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

double normalizeKillingVector(void *params,
                              const double thetap,
                              const double phip){
/*
double normalizeKillingVector(const DataMesh& theta,
                              const DataMesh& phi,
                              const SurfaceBasis& sb,
                              const DataMesh& Psi,
                              DataMesh& v,
                              const double thetap,
                              const double phip)
{
*/
  const StrahlkorperWithMesh& skwm = static_cast<struct rparams*>(params)->skwm;
  const DataMesh& Psi = static_cast<struct rparams*>(params)->Psi;
  DataMesh& v = static_cast<struct rparams*>(params)->v;
  const SurfaceBasis sb(skwm.Grid());

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
  double ds; //affine path length
  double t; //affine path length

  std::cout << "GSL version" << std::endl;
  bool goodtheta = KillingPathNew(skwm, rotated_Psi, xi, t, M_PI/2.0);
  REQUIRE(goodtheta, "Killing trajectory did not close " << POSITION);
  const double scale1 = t/(2.0*M_PI);
  std::cout << "Theta = " << std::setprecision(8) << std::setw(10)
            << M_PI/2.0 << " : T = "
            << std::setprecision(10) << t
            << " (" << scale1 << ")" << std::endl;


  std::cout << "old version" << std::endl;
  goodtheta = KillingPath(params, rotated_Psi, ds, M_PI/2.0, xi);
  REQUIRE(goodtheta, "Killing trajectory did not close " << POSITION);
  const double scale = ds/(2.0*M_PI);
  std::cout << "Theta = " << std::setprecision(8) << std::setw(10)
            << M_PI/2.0 << " : T = "
            << std::setprecision(10) << ds
            << " (" << scale << ")" << std::endl;


  goodtheta = KillingPath(params, rotated_Psi, ds, 0.5*M_PI/2.0, xi);
  std::cout << "Theta = " << std::setprecision(8) << std::setw(10)
            << 0.5*M_PI/2.0 << " : T = "
            << std::setprecision(10) << ds
            << " (" << ds/(2.0*M_PI*scale) << ")" << std::endl;

  goodtheta = KillingPath(params, rotated_Psi, ds, 0.25*M_PI/2.0, xi);
  std::cout << "Theta = " << std::setprecision(8) << std::setw(10)
            << 0.25*M_PI/2.0 << " : T = "
            << std::setprecision(10) << ds
            << " (" << ds/(2.0*M_PI*scale) << ")" << std::endl;

  goodtheta = KillingPath(params, rotated_Psi, ds, 0.125*M_PI/2.0, xi);
  std::cout << "Theta = " << std::setprecision(8) << std::setw(10)
            << 0.125*M_PI/2.0 << " : T = "
            << std::setprecision(10) << ds
            << " (" << ds/(2.0*M_PI*scale) << ")" << std::endl;

  return scale;

} //end normalizeKillingVector

bool KillingPathNew(const StrahlkorperWithMesh& skwm,
                    const DataMesh& Psi_r,
                    const Tensor<DataMesh>& xi,
                    double& t,
                    const double theta)
{
  struct ODEparams params = {skwm, Psi_r, xi};

  const gsl_odeiv2_step_type *const T = gsl_odeiv2_step_rkck;
  gsl_odeiv2_step *const s = gsl_odeiv2_step_alloc (T, 2);
  gsl_odeiv2_control *const c = gsl_odeiv2_control_y_new (1e-12, 1.e-12);
  gsl_odeiv2_evolve *const e = gsl_odeiv2_evolve_alloc (2);
     
  gsl_odeiv2_system sys = {PathDerivNew, NULL, 2, &params};
     
  t = 0.0;
  double t1 = 1.e10;
  double h = 1e-6;
  double y[2] = { theta, 0.0 };
  bool limit_h = false;
  double hmax = h;
  const double hmax2 = 2.0*M_PI/100.0;

  std::cout << "START: y = ( " 
	    << std::setprecision(8) << std::setw(10) << y[0] << " , " 
	    << std::setprecision(8) << std::setw(10) << y[1] << " ); h = " 
	    << std::setprecision(8) << std::setw(10) << h 
	    << std::endl; 
     
  while (true) {
    const double ysave[2] = {y[0],y[1]}; 
    const double tsave = t;
    const int status = gsl_odeiv2_evolve_apply (e, c, s,
						&sys, 
						&t, t1,
						&h, y);
    std::cout << "t = " 
	      << std::setprecision(8) << std::setw(10) << t 
	      << " : y = ( " 
	      << std::setprecision(8) << std::setw(10) << y[0] << " , " 
	      << std::setprecision(8) << std::setw(10) << y[1] << " ); h = " 
	      << std::setprecision(8) << std::setw(10) << h 
	      << std::endl;

    ASSERT(status==GSL_SUCCESS,"Path Integration failed");
    if(limit_h && h > hmax) h = hmax;
    if(fabs(y[1] - 2.*M_PI) < 1.e-10) break;
    else if(y[1] > 2.*M_PI) {
      if(!limit_h) hmax = h;
      limit_h = true;
      h = hmax *= 0.5;
      y[0] = ysave[0]; y[1] = ysave[1]; t = tsave;
      gsl_odeiv2_evolve_reset(e);
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

//need to rename this function
int PathDerivNew(double t, const double y[], double f[], void *params){
  const StrahlkorperWithMesh& skwm = static_cast<struct ODEparams*>(params)->skwm;
  const DataMesh& Psi_r = static_cast<struct ODEparams*>(params)->Psi_r;
  const Tensor<DataMesh>& xi = static_cast<struct ODEparams*>(params)->xi;

  const SurfaceBasis sb(skwm.Grid());
  const DataMesh& rad = skwm.Radius();

  double psiAtPoint = sb.Evaluate(Psi_r, y[0], y[1]);
  double radAtPoint = sb.Evaluate(rad, y[0], y[1]);

  //reference YlmSpherepack directly?
  MyVector<double> result = sb.EvaluateVector(xi, y[0], y[1]);

  const double norm = 1.0 / (psiAtPoint*psiAtPoint*psiAtPoint*psiAtPoint
                             *radAtPoint*radAtPoint);
  f[0] = result[0]*norm;
  f[1] = result[1]*norm/sin(y[0]);

  return GSL_SUCCESS;
}


bool KillingPath(void *params,
                 const DataMesh& Psi, //rotated Psi
                 double& ds,
                 const double theta,
                 const Tensor<DataMesh>& xi){
  const int n = 100; //number of steps around the circumference (when theta=Pi/2)

  double hmax = 2.0*M_PI/n; //stepsize around circumference
  double h = hmax;
  MyVector<double> Vin(MV::fill, theta, 0.0); //(theta,phi) coords
  MyVector<double> Vout(MV::Size(2));
  ds = 0.0; //affine path length

  int counter = 0;
  while(true){
    counter++;

    Vout = Vin;

    MyVector<double> Vp = PathDerivs(params, Psi, Vout, xi); //vector harmonic synthesis at a point

    MyVector<double> scale(MV::Size(2));
    for(int i=0; i<Vout.Size(); ++i){
      //takes current position, adds some fraction from the scaled Killing vector
      scale[i] = fabs(Vout[i]) + fabs(h*Vp[i]) + 1.e-30;
    }

    double hdid = 0.0;

    //std::cout << "ds = " << std::setprecision(8) << std::setw(10) << ds << std::endl;

    PathRKQC(ds, h, hmax, hdid, Vout, Vp, params, Psi, xi, scale, 1.e-12);

    //std::cout << "ds = " << std::setprecision(8) << std::setw(10) << ds 
              //<< " Vout[1]-2*Pi = " << std::setprecision(8) 
              //<< std::setw(10) << Vout[1] - 2.*M_PI << std::endl;

    if(fabs(Vout[1] - 2.*M_PI) < 1.e-10){
      ds += hdid;
      Vin[0] = Vout[0];
      Vin[1] = Vout[1];
      break;
    } else if(Vout[1] > 2.*M_PI){
      hmax = h *= 0.5;
    } else {
      ds += hdid;
      Vin[0] = Vout[0];
      Vin[1] = Vout[1];
    } //end ifs
  } // end while

  const bool closedPath = (fabs(Vout[0] - theta) < 1.e-6);
  if(!closedPath){
    std::cout << "##> Theta diff: " << std::setprecision(6) << Vout[0] - theta << std::endl;
  }

  return closedPath;
} //end KillingPath


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
              const double epsilon)
{
  //Runge-Kutta Quality Control
  while(true){
    MyVector<double> Vout(MV::Size(2));
    MyVector<double> error(MV::Size(2));
    hdid = h;

    PathRKCK(h,Vin,Vout,Vp,params,Psi,xi,error);
    double errmax = 0.0;

    for(int i=0; i<error.Size(); ++i){
      const double relerr = fabs(error[i]/scale[i]);
      if(relerr > errmax) errmax = relerr;
    }

    errmax /= epsilon;

    if(errmax > 1.0) { //failed step, reduce h
      double hnext = 0.9*h*pow(errmax, -0.25); //decrease h
      if(h/hnext > 10.0) hnext = 0.1*h; //by not by more than 0.1
      const double xtest = ds + (0.1*hnext);
      if(ds == xtest) { //can't step smaller so accept
        for(int i=0; i<2; ++i) Vin[i] = Vout[i];
        break;
      } else { // accept new stepsize and try again
        h = hnext;
        continue;
      }
    } else {
      double hnext = 5.0*h; //increase h, but not by more than 5
      if(errmax > 1.89e-4) hnext = 0.9*h*pow(errmax, -0.2); //increase h
      if(hnext > hmax) hnext = hmax; //and not more than hmax
      //accept and move on
      for(int i=0; i<Vin.Size(); ++i) Vin[i] = Vout[i];

      h = hnext;
      break;
    } //end if
  } //end while
} //end PathRKQC

//returns single position along path
MyVector<double> PathDerivs(void *params,
                const DataMesh& Psi, //rotated Psi
                const MyVector<double>& Vin,
                const Tensor<DataMesh>& xi) {

  const StrahlkorperWithMesh& skwm = static_cast<struct rparams*>(params)->skwm;

  const SurfaceBasis sb(skwm.Grid());
  const DataMesh& rad = skwm.Radius();

  double psiAtPoint = sb.Evaluate(Psi, Vin[0], Vin[1]);
  double radAtPoint = sb.Evaluate(rad, Vin[0], Vin[1]);

  //MyVector<double> result(MV::Size(2),0.0);
  //result = sb.EvaluateVector(xi, Vin[0], Vin[1]);
  MyVector<double> result = sb.EvaluateVector(xi, Vin[0], Vin[1]);

  const double norm = 1.0 / (psiAtPoint*psiAtPoint*psiAtPoint*psiAtPoint
                             *radAtPoint*radAtPoint);
  result[0] = result[0]*norm;
  result[1] = result[1]*norm/sin(Vin[0]);

  return result;
} //end PathDerivs


//changes error, Vout
void PathRKCK(const double& h,
              const MyVector<double>& Vin,
              MyVector<double>& Vout,
              const MyVector<double>& Vp,
              void *params,
              const DataMesh& Psi, //rotated Psi
              const Tensor<DataMesh>& xi,
              MyVector<double>& error)
{
  //This uses the Cash-Karp version of the Runge-Kutta method
  //for solving ordinary differential equations

  //first step
  MyVector<double> t(MV::Size(2));
  for(int j=0; j<t.Size(); ++j) t[j] = Vin[j] + 0.2*Vp[j];

  //second step
  MyVector<double> k2 = PathDerivs(params, Psi, t, xi);
  for(int j=0; j<k2.Size(); j++){
    t[j] = Vin[j] + 0.075*h*(  Vp[j]
                             + 3.0*k2[j]);
  }

  //third step
  MyVector<double> k3 = PathDerivs(params, Psi, t, xi);
  for(int j=0; j<k2.Size(); j++){
    t[j] = Vin[j] + 0.3*h*(  Vp[j]
                           - 3.0*k2[j]
                           + 4.0*k3[j]);
  }

  //fourth step
  MyVector<double> k4 = PathDerivs(params, Psi, t, xi);
  for(int j=0; j<k2.Size(); j++){
    t[j] = Vin[j] + h*( -11.0*Vp[j]
                       + 135.0*k2[j]
                       - 140.0*k3[j]
                       + 70.0*k4[j]) / 54.0;
  }

  //fifth step
  MyVector<double> k5 = PathDerivs(params, Psi, t, xi);
  for(int j=0; j<k2.Size(); j++){
    t[j] = Vin[j] + h*(  3262.0*Vp[j]
                       + 37800.0*k2[j]
                       + 4600.0*k3[j]
                       + 44275.0*k4[j]
                       + 6831.0*k5[j]) / 110592.0;
  }

  //sixth step
  MyVector<double> k6 = PathDerivs(params, Psi, t, xi);
  for(int j=0; j<k2.Size(); j++){
    Vout[j] = Vin[j] + h*(  (37.0/378.0)*Vp[j]
                          + (250.0/621.0)*k3[j]
                          + (125.0/594.0)*k4[j]
                          + (512.0/1771.0)*k6[j] );
  }

  //error estimate
  for(int j=0; j<error.Size(); j++){
    error[j] = h*( -(277.0/64512.0)*Vp[j]
                  + (6925.0/370944.0)*k3[j]
                  - (6925.0/202752.0)*k4[j]
                  - (277.0/14336.0)*k5[j]
                  + (277.0/7084.0)*k6[j] );
  }

}

void KillingDiagnostics(const StrahlkorperWithMesh& skwm,
                        const DataMesh& L,
                        const DataMesh& Psi, //original
                        const Tensor<DataMesh>& xi,
                        const MyVector<bool>& printDiagnostic){

  const SurfaceBasis sb(skwm.Grid());
  const DataMesh& rad = skwm.Radius();
  const DataMesh& p2r2 = rad*rad*Psi*Psi;

//for testing only below here------------------------------
    //const int mNth = skwm.Grid().SurfaceCoords()(0).Extents()[0];
    //const int mNph = skwm.Grid().SurfaceCoords()(0).Extents()[1];
//for testing only above here-----------------------------------------------

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







