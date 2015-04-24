#ifndef _ENFW_H_
#define _ENFW_H_

#include <iostream>
#include <assert.h>
#include <cmath>
#include "../cosmology.h"

using namespace std;

double pow2(double a) { return a*a; }
double pow3(double a) { return a*a*a; }

// conventions:
//   r is always the profile effective radius (ellipsoidally transformed, scaled with scale radius); 
//     it has nothing to do with a physical radius
//   R, phi, theta is the spherical coordinate system with R in units of scale radius
//   rho(R) is a density at R
//   Rho(R) is the mean density inside a sphere of radius R
//   e = sqrt(1-q^2) is the eccentricity if anybody asks

/////////////////////////////////////// 3D DENSITY WITH RHO_0=1 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

double nfw_rho(double r)
// returns 3d NFW density at r=R/R_s, assuming rho_0=1
{
  assert(r>0);
  return 1./(r*pow2(1.+r));
}



double enfw_rho(double R, double theta, double am2, double cm2)
// returns 3d NFW density of prolate ellipsoidal halo with ellipticity am2=a^-2=(1-e^2)^(-1/3), em2=c^-2=(1-e^2)^(2/3)
// where e is the eccentricity
// theta=0 is the major axis (z)
// theta=pi/2 is the minor axis (x,y)
// the result is invariant under rotations around z (phi-invariant)
{
  assert(R>0);
  assert(am2>0);
  assert(cm2>0);
  
  double z=R*cos(theta);
  double xy=R*sin(theta); // (xy)^2=x^2+y^2
  //double y=xy*sin(phi);
  //double x=xy*cos(phi);
  
  // calculate ellipsoidal equivalent radius
  // an ellipsoid of volume 4/3 pi r^3 is described by the equation
  //
  //             ( (1-e^2)^(-1/3)       0              0       )   ( x )
  // (x, y, z) * (       0        (1-e^2)^(-1/3)       0       ) * ( y )  =  r^2
  //             (       0              0        (1-e^2)^(2/3) )   ( z )
  //
  // this r is what we should pass to nfw_rho
  
  
  double r = sqrt((xy*xy)*am2+z*z*cm2);
  
  return nfw_rho(r);  
}


double enfw_rho_trunc(double R, double theta, double am2, double cm2, double tau)
{
  assert(R>0);
  assert(am2>0);
  assert(cm2>0);
  
  double z=R*cos(theta);
  double xy=R*sin(theta); 
  
  double r = sqrt((xy*xy)*am2+z*z*cm2);
  
  return nfw_rho(r)*pow(tau*tau/(r*r+tau*tau),2);  
}

double enfw_rho_trunc3d(double R, double theta, double am2, double cm2, double tau)
{
  assert(R>0);
  assert(am2>0);
  assert(cm2>0);
  
  double z=R*cos(theta);
  double xy=R*sin(theta); 
  
  double r = sqrt((xy*xy)*am2+z*z*cm2);
  
  return nfw_rho(r)*pow(tau*tau/(R*R+tau*tau),2);  
}


double enfw_rho(double R, double theta, double e)
// returns 3d NFW density of prolate ellipsoidal halo with ellipticity e
// theta=0 is the major axis (z)
// theta=pi/2 is the minor axis (x,y)
// the result is invariant under rotations around z (phi-invariant)
{
  double am2=(1.-e*e);
  double cm2=pow(am2,2./3.);
  am2=pow(am2,-1./3.);
  
  return enfw_rho(R, theta, am2, cm2);
}

double enfw_rho_trunc(double R, double theta, double e, double tau)
// returns 3d NFW density of prolate ellipsoidal halo with ellipticity e
// theta=0 is the major axis (z)
// theta=pi/2 is the minor axis (x,y)
// the result is invariant under rotations around z (phi-invariant)
{
  double am2=(1.-e*e);
  double cm2=pow(am2,2./3.);
  am2=pow(am2,-1./3.);
  
  return enfw_rho_trunc(R, theta, am2, cm2, tau);
}

double enfw_rho_trunc3d(double R, double theta, double e, double tau)
// returns 3d NFW density of prolate ellipsoidal halo with ellipticity e
// theta=0 is the major axis (z)
// theta=pi/2 is the minor axis (x,y)
// the result is invariant under rotations around z (phi-invariant)
{
  double am2=(1.-e*e);
  double cm2=pow(am2,2./3.);
  am2=pow(am2,-1./3.);
  
  return enfw_rho_trunc3d(R, theta, am2, cm2, tau);
}


double enfw_rho_ac(double R, double am2, double cm2)
// numerically averaged elliptical NFW density at physical radius R
{
  // integrate one octant phi,theta=[0,pi/2]
  // but luckily phi does not matter at all for a prolate halo
 
  double rho=0.;
 
  // integrating dcostheta
  
  // number of steps 
  const int ncostheta=50;
  // step width
  const double dcostheta=1./double(ncostheta);
 
  
  for(int i=0; i<ncostheta; i++)
  {
   double theta=acos(dcostheta*(double(i)+0.5));
   rho += enfw_rho(R, theta, am2, cm2)*dcostheta;
  }
  
  return rho;
}


double enfw_rho_ac_trunc(double R, double am2, double cm2, double tau)
// numerically averaged elliptical NFW density at physical radius R
{
  // integrate one octant phi,theta=[0,pi/2]
  // but luckily phi does not matter at all for a prolate halo
 
  double rho=0.;
 
  // integrating dcostheta
  
  // number of steps 
  const int ncostheta=50;
  // step width
  const double dcostheta=1./double(ncostheta);
 
  
  for(int i=0; i<ncostheta; i++)
  {
   double theta=acos(dcostheta*(double(i)+0.5));
   rho += enfw_rho_trunc(R, theta, am2, cm2, tau)*dcostheta;
  }
  
  return rho;
}

double enfw_rho_ac_trunc3d(double R, double am2, double cm2, double tau)
// numerically averaged elliptical NFW density at physical radius R
{
  // integrate one octant phi,theta=[0,pi/2]
  // but luckily phi does not matter at all for a prolate halo
 
  double rho=0.;
 
  // integrating dcostheta
  
  // number of steps 
  const int ncostheta=50;
  // step width
  const double dcostheta=1./double(ncostheta);
 
  
  for(int i=0; i<ncostheta; i++)
  {
   double theta=acos(dcostheta*(double(i)+0.5));
   rho += enfw_rho_trunc3d(R, theta, am2, cm2, tau)*dcostheta;
  }
  
  return rho;
}




double enfw_rho(double R, double e)
// numerically averaged elliptical NFW density at physical radius R
{
  // integrate one octant phi,theta=[0,pi/2]
  // but luckily phi does not matter at all for a prolate halo
  
  double am2=(1.-e*e);
  const double cm2=pow(am2,2./3.);
  am2=pow(am2,-1./3.);
  
  return enfw_rho_ac(R, am2, cm2);
}


double enfw_rho_trunc(double R, double e, double tau)
// numerically averaged elliptical NFW density at physical radius R
{
  // integrate one octant phi,theta=[0,pi/2]
  // but luckily phi does not matter at all for a prolate halo
  
  double am2=(1.-e*e);
  const double cm2=pow(am2,2./3.);
  am2=pow(am2,-1./3.);
  
  return enfw_rho_ac_trunc(R, am2, cm2, tau);
}


double enfw_rho_trunc3d(double R, double e, double tau)
// numerically averaged elliptical NFW density at physical radius R
{
  // integrate one octant phi,theta=[0,pi/2]
  // but luckily phi does not matter at all for a prolate halo
  
  double am2=(1.-e*e);
  const double cm2=pow(am2,2./3.);
  am2=pow(am2,-1./3.);
  
  return enfw_rho_ac_trunc3d(R, am2, cm2, tau);
}










////////////////////////////////////////// SPHERICALLY AVERAGED DENSITY WITH RHO_0=1 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


double nfw_Rho(double r)
// returns mean NFW density inside sphere of r=R/R_s, assuming rho_0=1
{
  assert(r>0);
  return 3.*(1./(r+1.)-1.+log(r+1.))/pow3(r);
}

double enfw_Rho(double R, double e)
// numerically averaged density of prolate NFW halo inside sphere of radius R
// this is not the best version because we neglect all matter inside 0.001R_S
{

  const double lR0=-3.;
  const double lR1=log10(R);
  assert(lR1>lR0);
  
  const int nlR=int((lR1-lR0)/0.005); // roughly 1% steps in radius
  double m=0.;
  const double dlR=(lR1-lR0)/double(nlR);
  
  double am2=(1.-e*e);
  const double cm2=pow(am2,2./3.);
  am2=pow(am2,-1./3.);
  
  const double f = 4./3.*M_PI*(pow(10,3.*dlR/2.)-pow(10,-3.*dlR/2.)); // a handy factor for numerical integration below
  
  for(int i=0; i<nlR; i++)
  {
   double Rc=pow(10,lR0+(double(i)+0.5)*dlR); // central log10(radius)
   double rho=enfw_rho_ac(Rc,am2,cm2);
   m += rho*pow3(Rc)*f;
  }
  
  return m/(4./3.*M_PI*pow3(R));
  
}


double enfw_Rho_corecorr(double R, double e)
// numerically averaged density of prolate NFW halo inside sphere of radius R
// this makes a correction for the core (R/R_s<0.001) region, including a correct scaling of the correction with ellipticity
{

  const double lR0=-3.;
  const double lR1=log10(R);
  assert(lR1>=lR0);
  
  double am2=(1.-e*e);
  const double cm2=pow(am2,2./3.);
  am2=pow(am2,-1./3.);
  
  const int nlR=int((lR1-lR0)/0.001); // roughly 0.2% steps in radius TODO: set to 0.005, that's enough
  double m=nfw_Rho(1.e-3)*4./3.*M_PI*pow3(1.e-3);//*enfw_rho_ac(1.e-3,am2,cm2)/nfw_rho(1.e-3);
  const double dlR=(lR1-lR0)/double(nlR);
  
  
  const double f = 4./3.*M_PI*(pow(10,3.*dlR/2.)-pow(10,-3.*dlR/2.)); // a handy factor for numerical integration below
  
  for(int i=0; i<nlR; i++)
  {
   double Rc=pow(10,lR0+(double(i)+0.5)*dlR); // central log10(radius)
   double rho=enfw_rho_ac(Rc,am2,cm2);
   m += rho*pow3(Rc)*f;
  }
  
  return m/(4./3.*M_PI*pow3(R));
}








////////////////////////////////////// RHO_0 for halo of given concentration, redshift and ellipticity \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



double enfw_rho0_200m(double c200m, double z, double e)
// scale density, to be multiplied with the above nfw_rho (depends only on concentration, the mass dependence is in the R/R_s)
// at the chosen value of h
{
 double rhom  = OmegaM*3.*H0*H0/(8.*M_PI*Gs)*pow3(1.+z); // Msol/Mpc^3, at the chosen value of h
 double rho0s = 200.*rhom/enfw_Rho_corecorr(c200m,e);  // Msol/R_s^3
 return rho0s;
}

double nfw_rho0_200m(double c200m, double z)
// scale density, to be multiplied with the above nfw_rho (depends only on concentration, the mass dependence is in the R/R_s)
// at the chosen value of h
{
  double rhom  = OmegaM*3.*H0*H0/(8.*M_PI*Gs)*pow3(1.+z);   // Msol/Mpc^3
  return 200.*rhom*pow3(c200m)/ (3.*(log(1.+c200m)-c200m/(1.+c200m)));
  
  // same result:
  //double r200m = pow(3./(4.*M_PI*200.*rhom),1./3.);       // Mpc for a mass of 1Msol
  //double rs    = r200m/c200m;                             // Mpc
  //double rho0s = 1./(4.*M_PI*rs*rs*rs*(log(1.+c200m)-c200m/(1.+c200m)));   // Msol/R_s^3
  //return rho0s;
}






///////////////////////////////////////// cylindrical projection \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


double nfw_rhop(double v)
// projected density (analytical) of spherical NFW at v[R_s] distance from centre
// this is missing a normalization for R_s/Mpc*rho_0
{
  assert(v>0);
  
  if(v<1.)
  {
	  return 2./(v*v-1.) * ( 1. - 2./sqrt(1.-v*v)*atanh( sqrt((1.-v)/(1.+v)) ) );
  }
  else if(v==1.)
  {
	  return 2./3.;
  }

  return 2./(v*v-1.) * ( 1. - 2./sqrt(v*v-1.)*atan( sqrt((v-1.)/(1.+v)) ) );
 
}

double enfw_rho_cylindrical(double w, double v, double phi, double alpha, double e)
// evaluated in a cylindrical coordinate system 
// with w axis tilted w.r.t major axis of halo by angle alpha and sitting in z,x plane, 
// radius v around that axis and azimuthal angle phi around that axis
{
 double x = w*sin(alpha)+v*cos(alpha)*cos(phi);
 double y = v*sin(phi);
 double z = w*cos(alpha)-v*sin(alpha)*cos(phi);
 double R = sqrt(w*w+v*v);
 double theta = atan(sqrt(x*x+y*y)/z);
 return enfw_rho(R, theta, e);
}

double enfw_rho_cylindrical_trunc(double w, double v, double phi, double alpha, double e, double tau)
// evaluated in a cylindrical coordinate system 
// with w axis tilted w.r.t major axis of halo by angle alpha and sitting in z,x plane, 
// radius v around that axis and azimuthal angle phi around that axis
{
 double x = w*sin(alpha)+v*cos(alpha)*cos(phi);
 double y = v*sin(phi);
 double z = w*cos(alpha)-v*sin(alpha)*cos(phi);
 double R = sqrt(w*w+v*v);
 double theta = atan(sqrt(x*x+y*y)/z);
 return enfw_rho_trunc(R, theta, e, tau);
}

double enfw_rho_cylindrical_trunc3d(double w, double v, double phi, double alpha, double e, double tau)
// evaluated in a cylindrical coordinate system 
// with w axis tilted w.r.t major axis of halo by angle alpha and sitting in z,x plane, 
// radius v around that axis and azimuthal angle phi around that axis
{
 double x = w*sin(alpha)+v*cos(alpha)*cos(phi);
 double y = v*sin(phi);
 double z = w*cos(alpha)-v*sin(alpha)*cos(phi);
 double R = sqrt(w*w+v*v);
 double theta = atan(sqrt(x*x+y*y)/z);
 return enfw_rho_trunc3d(R, theta, e, tau);
}


double enfw_rhop(double v, double phi, double alpha, double e, double lwmin=-4, double lwmax=4, bool estimate_inner=true)
  // this is missing a normalization for R_s/Mpc*rho_0
  // power-law trapezoidal
{
 int    nlw=400;
 double lwstep=(lwmax-lwmin)/double(nlw);
 
 double ln10=log(10.);
 
 double sigma=0;
 if(estimate_inner) sigma=pow(10,lwmin)*(enfw_rho_cylindrical(0, v, phi, alpha, e)+enfw_rho_cylindrical(pow(10,lwmin), v, phi, alpha, e)); 
 // integrated surface density inside innermost part, approximated
 
 double lw0=lwmin;
 double w0=pow(10,lw0);
 double rho0=enfw_rho_cylindrical(w0,v,phi,alpha,e);
 double rho0m=enfw_rho_cylindrical(-w0,v,phi,alpha,e);
 
 double w1w0=pow(10,lwstep);
 
 for(int i=1; i<=nlw; i++)
 {
  double w1=pow(10,lwmin+double(i)*lwstep);
  double rho1=enfw_rho_cylindrical(w1, v, phi, alpha, e);
  double rho1m=enfw_rho_cylindrical(-w1, v, phi, alpha, e);
  
  double alpha=log10(rho1/rho0)/lwstep;    // local power-law slope
  double alpham=log10(rho1m/rho0m)/lwstep; // local power-law slope
  
  sigma += rho0*w0/(alpha+1.)*(pow(w1w0,alpha+1)-1.);
  sigma += rho0m*w0/(alpham+1.)*(pow(w1w0,alpham+1)-1.);
  
  w0=w1;
  rho0=rho1;
  rho0m=rho1m;
 }

 return sigma; // this is missing a normalization for R_s/Mpc*rho_0
}


double enfw_rhop_trunc(double v, double phi, double alpha, double e, double tau, double lwmin=-4, double lwmax=4, bool estimate_inner=true)
  // this is missing a normalization for R_s/Mpc*rho_0
  // power-law trapezoidal
{
 int    nlw=400;
 double lwstep=(lwmax-lwmin)/double(nlw);
 
 double ln10=log(10.);
 
 double sigma=0;
 if(estimate_inner) sigma=pow(10,lwmin)*(enfw_rho_cylindrical_trunc(0, v, phi, alpha, e, tau)+enfw_rho_cylindrical_trunc(pow(10,lwmin), v, phi, alpha, e, tau)); 
 // integrated surface density inside innermost part, approximated
 
 double lw0=lwmin;
 double w0=pow(10,lw0);
 double rho0=enfw_rho_cylindrical_trunc(w0,v,phi,alpha,e,tau);
 double rho0m=enfw_rho_cylindrical_trunc(-w0,v,phi,alpha,e,tau);
 
 double w1w0=pow(10,lwstep);
 
 for(int i=1; i<=nlw; i++)
 {
  double w1=pow(10,lwmin+double(i)*lwstep);
  double rho1=enfw_rho_cylindrical_trunc(w1, v, phi, alpha, e, tau);
  double rho1m=enfw_rho_cylindrical_trunc(-w1, v, phi, alpha, e, tau);
  
  double alpha=log10(rho1/rho0)/lwstep;    // local power-law slope
  double alpham=log10(rho1m/rho0m)/lwstep; // local power-law slope
  
  sigma += rho0*w0/(alpha+1.)*(pow(w1w0,alpha+1)-1.);
  sigma += rho0m*w0/(alpham+1.)*(pow(w1w0,alpham+1)-1.);
  
  w0=w1;
  rho0=rho1;
  rho0m=rho1m;
 }

 return sigma; // this is missing a normalization for R_s/Mpc*rho_0
}


double enfw_rhop_trunc3d(double v, double phi, double alpha, double e, double tau, double lwmin=-4, double lwmax=4, bool estimate_inner=true)
  // this is missing a normalization for R_s/Mpc*rho_0
  // power-law trapezoidal
{
 int    nlw=400;
 double lwstep=(lwmax-lwmin)/double(nlw);
 
 double ln10=log(10.);
 
 double sigma=0;
 if(estimate_inner) sigma=pow(10,lwmin)*(enfw_rho_cylindrical_trunc3d(0, v, phi, alpha, e, tau)+enfw_rho_cylindrical_trunc3d(pow(10,lwmin), v, phi, alpha, e, tau)); 
 // integrated surface density inside innermost part, approximated
 
 double lw0=lwmin;
 double w0=pow(10,lw0);
 double rho0=enfw_rho_cylindrical_trunc3d(w0,v,phi,alpha,e,tau);
 double rho0m=enfw_rho_cylindrical_trunc3d(-w0,v,phi,alpha,e,tau);
 
 double w1w0=pow(10,lwstep);
 
 for(int i=1; i<=nlw; i++)
 {
  double w1=pow(10,lwmin+double(i)*lwstep);
  double rho1=enfw_rho_cylindrical_trunc3d(w1, v, phi, alpha, e, tau);
  double rho1m=enfw_rho_cylindrical_trunc3d(-w1, v, phi, alpha, e, tau);
  
  double alpha=log10(rho1/rho0)/lwstep;    // local power-law slope
  double alpham=log10(rho1m/rho0m)/lwstep; // local power-law slope
  
  sigma += rho0*w0/(alpha+1.)*(pow(w1w0,alpha+1)-1.);
  sigma += rho0m*w0/(alpham+1.)*(pow(w1w0,alpham+1)-1.);
  
  w0=w1;
  rho0=rho1;
  rho0m=rho1m;
 }

 return sigma; // this is missing a normalization for R_s/Mpc*rho_0
}


////////////////////////////////////// radial projected profile \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

double nfw_kappa(double v)
// projected density (analytical) of spherical NFW at v[R_s] distance from centre
// this is missing a normalization for R_s/Mpc*rho_0 and for Sigma_crit
{
 return nfw_rhop(v);
}

double enfw_kappa(double v, double alpha, double e)
// projected density of elliptical (e>0) NFW halo at v[R_s] projected circular distance from centre
// alpha is the tilt angle between projection axis and major axis of halo
// this is missing a normalization for R_s/Mpc*rho_0 and for Sigma_crit
{
  const int nphi = 100;
  const double phistep = M_PI/2./double(nphi);
  
  double sigma=0;
  
  for(int i=0; i<nphi; i++)
  {
   double phi=(double(i)+0.5)*phistep;
   sigma += enfw_rhop(v,phi,alpha,e);
  }
  return sigma/double(nphi);
}

double enfw_kappa_trunc(double v, double alpha, double e, double tau)
// projected density of elliptical (e>0) NFW halo at v[R_s] projected circular distance from centre
// alpha is the tilt angle between projection axis and major axis of halo
// this is missing a normalization for R_s/Mpc*rho_0 and for Sigma_crit
{
  const int nphi = 100;
  const double phistep = M_PI/2./double(nphi);
  
  double sigma=0;
  
  for(int i=0; i<nphi; i++)
  {
   double phi=(double(i)+0.5)*phistep;
   sigma += enfw_rhop_trunc(v,phi,alpha,e,tau);
  }
  return sigma/double(nphi);
}


double enfw_kappa_trunc3d(double v, double alpha, double e, double tau)
// projected density of elliptical (e>0) NFW halo at v[R_s] projected circular distance from centre
// alpha is the tilt angle between projection axis and major axis of halo
// this is missing a normalization for R_s/Mpc*rho_0 and for Sigma_crit
{
  const int nphi = 100;
  const double phistep = M_PI/2./double(nphi);
  
  double sigma=0;
  
  for(int i=0; i<nphi; i++)
  {
   double phi=(double(i)+0.5)*phistep;
   sigma += enfw_rhop_trunc3d(v,phi,alpha,e,tau);
  }
  return sigma/double(nphi);
}


//////////////////// mean projected density inside radius \\\\\\\\\\\\\\\\\


double nfw_Rhop(double v)
{
 assert(v>0);
 
 if(v<1) {
  return 4./(v*v)*(2./sqrt(1.-v*v)*atanh(sqrt((1.-v)/(1.+v))) + log(v/2.)); 
 }
 if(v==1) {
  return 4.*(1.+log(0.5)); 
 }
 
 return 4./(v*v)*(2./sqrt(v*v-1.)*atan(sqrt((v-1.)/(1.+v))) + log(v/2.));
  
}


double nfw_Rhop_trunc(double v, double tau)
{
 assert(v<tau/100); // small scale limit implemented only
 
 return nfw_Rhop(v); // not the greatest approximation here!
}


double enfw_Rhop(double v, double alpha, double e, double v0, double &rhop0)
// numerically determine mean density inside v; 
// use mean density inside v0 (rhop0) if this is known already
{
  double rho=0.; // this will carry the total mass
  double lvmin;
  double lvmax=log10(v);
  
  if(v0>0) {
   rho=M_PI*v0*v0*rhop0; 
   lvmin=log10(v0);
  } else {
   rho=M_PI*1.e-6*nfw_Rhop(1.e-3) / sqrt(pow2(sin(alpha))*pow(1.-e*e,-1./3.)+pow2(cos(alpha))*pow(1.-e*e,2./3.));
   lvmin=-3.;
  }
  
  if(lvmin==lvmax)
  {
     rhop0=rho/(M_PI*v*v);
     return rhop0;
  }
  
  assert(lvmin<lvmax);
 
  const int nlv = (lvmax-lvmin)/0.02+1; // 5% steps in radius
  const double lvstep=(lvmax-lvmin)/double(nlv);
  
  double w0 = pow(10,lvmin);
  double f0 = 2.*M_PI*w0*enfw_kappa(w0,alpha,e);
  
  double w1w0=pow(10,lvstep);
  
  for(int i=1; i<=nlv; i++)
  {
   double lv=lvmin+double(i)*lvstep;
   double w1=pow(10,lv);
   double f1=enfw_kappa(w1,alpha,e)*2.*M_PI*w1;
   
   double alpha=log10(f1/f0)/lvstep;
   
   rho += f0*w0/(alpha+1.)*(pow(w1w0,alpha+1.)-1.);
   
   w0=w1;
   f0=f1;
  }
  
  rhop0=rho/(M_PI*v*v);
    
  return rhop0;
}


double enfw_Rhop_trunc(double v, double alpha, double e, double v0, double tau, double &rhop0)
// numerically determine mean density inside v; 
// use mean density inside v0 (rhop0) if this is known already
{
  double rho=0.; // this will carry the total mass
  double lvmin;
  double lvmax=log10(v);
  
  if(v0>0) {
   rho=M_PI*v0*v0*rhop0; 
   lvmin=log10(v0);
  } else {
   rho=M_PI*1.e-6*nfw_Rhop_trunc(1.e-3,tau) / sqrt(pow2(sin(alpha))*pow(1.-e*e,-1./3.)+pow2(cos(alpha))*pow(1.-e*e,2./3.));
   lvmin=-3.;
  }
  
  if(lvmin==lvmax)
  {
     rhop0=rho/(M_PI*v*v);
     return rhop0;
  }
  
  assert(lvmin<lvmax);
 
  const int nlv = (lvmax-lvmin)/0.02+1; // 5% steps in radius
  const double lvstep=(lvmax-lvmin)/double(nlv);
  
  double w0 = pow(10,lvmin);
  double f0 = 2.*M_PI*w0*enfw_kappa_trunc(w0,alpha,e,tau);
  
  double w1w0=pow(10,lvstep);
  
  for(int i=1; i<=nlv; i++)
  {
   double lv=lvmin+double(i)*lvstep;
   double w1=pow(10,lv);
   double f1=enfw_kappa_trunc(w1,alpha,e,tau)*2.*M_PI*w1;
   
   double alpha=log10(f1/f0)/lvstep;
   
   rho += f0*w0/(alpha+1.)*(pow(w1w0,alpha+1.)-1.);
   
   w0=w1;
   f0=f1;
  }
  
  rhop0=rho/(M_PI*v*v);
    
  return rhop0;
}


double enfw_Rhop_trunc3d(double v, double alpha, double e, double v0, double tau, double &rhop0)
// numerically determine mean density inside v; 
// use mean density inside v0 (rhop0) if this is known already
{
  double rho=0.; // this will carry the total mass
  double lvmin;
  double lvmax=log10(v);
  
  if(v0>0) {
   rho=M_PI*v0*v0*rhop0; 
   lvmin=log10(v0);
  } else {
   rho=M_PI*1.e-6*nfw_Rhop_trunc(1.e-3,tau) / sqrt(pow2(sin(alpha))*pow(1.-e*e,-1./3.)+pow2(cos(alpha))*pow(1.-e*e,2./3.));
   lvmin=-3.;
  }
  
  if(lvmin==lvmax)
  {
     rhop0=rho/(M_PI*v*v);
     return rhop0;
  }
  
  assert(lvmin<lvmax);
 
  const int nlv = (lvmax-lvmin)/0.02+1; // 5% steps in radius
  const double lvstep=(lvmax-lvmin)/double(nlv);
  
  double w0 = pow(10,lvmin);
  double f0 = 2.*M_PI*w0*enfw_kappa_trunc3d(w0,alpha,e,tau);
  
  double w1w0=pow(10,lvstep);
  
  for(int i=1; i<=nlv; i++)
  {
   double lv=lvmin+double(i)*lvstep;
   double w1=pow(10,lv);
   double f1=enfw_kappa_trunc3d(w1,alpha,e,tau)*2.*M_PI*w1;
   
   double alpha=log10(f1/f0)/lvstep;
   
   rho += f0*w0/(alpha+1.)*(pow(w1w0,alpha+1.)-1.);
   
   w0=w1;
   f0=f1;
  }
  
  rhop0=rho/(M_PI*v*v);
    
  return rhop0;
}



//////////////////// delta Sigma = gamma \\\\\\\\\\\\\\\\\\\\\\\\\\\


double nfw_gamma(double v)
{           

        assert(v>0);
  
	if(v<=0.99)
	{
		return 2.*(log((1.+sqrt((1.-v)/(1.+v)))/(1.-sqrt((1.-v)/(1.+v)))))*(2./(v*v*sqrt(1.-v*v))
		      +1./((v*v-1.)*sqrt(1.-v*v))) + 4.*log(v/2.)/(v*v)-2./(v*v-1.);
	}
	
	if(v<1.01) // interpolate rather than calculate, because this might be numerically very unstable
	{
	  const double g099=0.563973487; 
	  const double g101=0.557543610;
	  return ((v-0.99)*g101+(1.01-v)*g099)/0.02;
	}

	return 4.*atan(sqrt((v-1.)/(1.+v))) * (2./(v*v*sqrt(v*v-1.))+1./(pow((v*v-1.),(3./2.)))) + 4.*log(v/2.)/(v*v)-2/(v*v-1); 
  
}

double enfw_gamma(double v, double alpha, double e, double v0, double& rhop0)
{
  // rhop0 reference is set to new value
  return enfw_Rhop(v,alpha,e,v0,rhop0)-enfw_kappa(v,alpha,e);
}


#endif
