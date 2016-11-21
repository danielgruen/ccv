#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include <cmath>
#include <assert.h>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <cstring>

// same cosmology as Becker & Kravtsov, similar to WMAP 7

// physical constants
const double G_si=6.67384E-11; // m^3/(kg*s^2)
const double ckms=299792.458;  // km/s

// fundamental cosmological parameters
const double OmegaM=0.27;
const double OmegaL=1.-OmegaM;
const double OmegaB=0.044;
const double sigma8=0.79;
const double h=0.7;
const double H0=h*100.; // km/s/Mpc
const double ns=0.95;

// derived
const double Gs=G_si*(1.9891E30)/(3.08567758E22)*1.E-6;  // = 4.25e-9;  // Mpc/Msol (km/s)^2
double rho0crit = 3.*H0*H0/(8.*M_PI*Gs); // [ Msol / Mpc^3 ]

// useful
const double arcminperrad=180.*60./M_PI;
const double radperarcmin=1./arcminperrad;

double pow10(double a) { return pow(10,a); }

double angularDiameterDistance(double z1, double z2, int num=1000)
{
   double sum = 0;
   double R2 = 1.0/(1.0+z1);
   double R1 = 1.0/(1.0+z2);
   double dR = (R2-R1)/((double)num);
   double R = R1;
   for (int i = 0 ; i < num ; i++, R += dR)
   {
      double term1 = OmegaL*(R*R*R*R); // constant
      double term2 = OmegaM*R;   // diluted matter

      double val1 = 1.0/sqrt(term1+term2);

      term1 = OmegaL*(R+dR)*(R+dR)*(R+dR)*(R+dR);
      term2 = OmegaM*(R+dR); 

      double val2 = 1.0/sqrt(term1+term2);

      sum += ((val1+val2)/2.0)*dR; // trapezium rule
      // note that da=a^2 dz; we multiply the a^2 into the square root terms
   }

   double result = sum*ckms/100./h/(1.0+z2); // 3000 MPc h^-1
   return result; // Mpc / rad  [ no h^{-1}! ]
}


namespace pklinear {

  // taken from http://background.uchicago.edu/~whu/transfer/power.c

  /* Fitting Formulae for CDM + Baryon + Massive Neutrino (MDM) cosmologies. */
  /* Daniel J. Eisenstein & Wayne Hu, Institute for Advanced Study */

  /* There are two primary routines here, one to set the cosmology, the
  other to construct the transfer function for a single wavenumber k. 
  You should call the former once (per cosmology) and the latter as 
  many times as you want. */

  /* TFmdm_set_cosm() -- User passes all the cosmological parameters as
	  arguments; the routine sets up all of the scalar quantites needed 
	  computation of the fitting formula.  The input parameters are: 
	  1) omega_matter -- Density of CDM, baryons, and massive neutrinos,
				  in units of the critical density. 
	  2) omega_baryon -- Density of baryons, in units of critical. 
	  3) omega_hdm    -- Density of massive neutrinos, in units of critical 
	  4) degen_hdm    -- (Int) Number of degenerate massive neutrino species 
	  5) omega_lambda -- Cosmological constant 
	  6) hubble       -- Hubble constant, in units of 100 km/s/Mpc 
	  7) redshift     -- The redshift at which to evaluate */

  /* TFmdm_onek_mpc() -- User passes a single wavenumber, in units of Mpc^-1.
	  Routine returns the transfer function from the Eisenstein & Hu
	  fitting formula, based on the cosmology currently held in the 
	  internal variables.  The routine returns T_cb (the CDM+Baryon
	  density-weighted transfer function), although T_cbn (the CDM+
	  Baryon+Neutrino density-weighted transfer function) is stored
	  in the global variable tf_cbnu. */

  /* We also supply TFmdm_onek_hmpc(), which is identical to the previous
	  routine, but takes the wavenumber in units of h Mpc^-1. */

  /* We hold the internal scalar quantities in global variables, so that
  the user may access them in an external program, via "extern" declarations. */

  /* Please note that all internal length scales are in Mpc, not h^-1 Mpc! */

  /* -------------------------- Prototypes ----------------------------- */

  int TFmdm_set_cosm(double omega_matter, double omega_baryon, double omega_hdm,
	  int degen_hdm, double omega_lambda, double hubble, double redshift);
  double TFmdm_onek_mpc(double kk);
  double TFmdm_onek_hmpc(double kk);



  /* Convenience from Numerical Recipes in C, 2nd edition */
  static double sqrarg;
  #define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

  /* ------------------------- Global Variables ------------------------ */

  /* The following are set in TFmdm_set_cosm() */
  double   alpha_gamma,	/* sqrt(alpha_nu) */
	  alpha_nu,	/* The small-scale suppression */
	  beta_c,		/* The correction to the log in the small-scale */
	  num_degen_hdm,	/* Number of degenerate massive neutrino species */
	  f_baryon,	/* Baryon fraction */
	  f_bnu,		/* Baryon + Massive Neutrino fraction */
	  f_cb,		/* Baryon + CDM fraction */
	  f_cdm,		/* CDM fraction */
	  f_hdm,		/* Massive Neutrino fraction */
	  growth_k0,	/* D_1(z) -- the growth function as k->0 */
	  growth_to_z0,	/* D_1(z)/D_1(0) -- the growth relative to z=0 */
	  hhubble,	/* Need to pass Hubble constant to TFmdm_onek_hmpc() */
	  k_equality,	/* The comoving wave number of the horizon at equality*/
	  obhh,		/* Omega_baryon * hubble^2 */
	  omega_curv,	/* = 1 - omega_matter - omega_lambda */
	  omega_lambda_z, /* Omega_lambda at the given redshift */
	  omega_matter_z,	/* Omega_matter at the given redshift */
	  omhh,		/* Omega_matter * hubble^2 */
	  onhh,		/* Omega_hdm * hubble^2 */
	  p_c,		/* The correction to the exponent before drag epoch */
	  p_cb,		/* The correction to the exponent after drag epoch */
	  sound_horizon_fit,  /* The sound horizon at the drag epoch */
	  theta_cmb,	/* The temperature of the CMB, in units of 2.7 K */
	  y_drag,		/* Ratio of z_equality to z_drag */
	  z_drag,		/* Redshift of the drag epoch */
	  z_equality;	/* Redshift of matter-radiation equality */

  /* The following are set in TFmdm_onek_mpc() */



  /* By default, these functions return tf_cb */

  /* ------------------------- TFmdm_set_cosm() ------------------------ */
  int TFmdm_set_cosm(double omega_matter, double omega_baryon, double omega_hdm,
	  int degen_hdm, double omega_lambda, double hubble, double redshift)
  /* This routine takes cosmological parameters and a redshift and sets up
  all the internal scalar quantities needed to compute the transfer function. */
  /* INPUT: omega_matter -- Density of CDM, baryons, and massive neutrinos,
				  in units of the critical density. */
  /* 	  omega_baryon -- Density of baryons, in units of critical. */
  /* 	  omega_hdm    -- Density of massive neutrinos, in units of critical */
  /* 	  degen_hdm    -- (Int) Number of degenerate massive neutrino species */
  /*        omega_lambda -- Cosmological constant */
  /* 	  hubble       -- Hubble constant, in units of 100 km/s/Mpc */
  /*        redshift     -- The redshift at which to evaluate */
  /* OUTPUT: Returns 0 if all is well, 1 if a warning was issued.  Otherwise,
	  sets many global variables for use in TFmdm_onek_mpc() */
  {
      double z_drag_b1, z_drag_b2, omega_denom;
      int qwarn;
      qwarn = 0;

      theta_cmb = 2.728/2.7;	/* Assuming T_cmb = 2.728 K */

      /* Look for strange input */
      if (omega_baryon<0.0) {
	  fprintf(stderr,
	    "TFmdm_set_cosm(): Negative omega_baryon set to trace amount.\n");
	  qwarn = 1;
      }
      if (omega_hdm<0.0) {
	  fprintf(stderr,
	    "TFmdm_set_cosm(): Negative omega_hdm set to trace amount.\n");
	  qwarn = 1;
      }
      if (hubble<=0.0) {
	  fprintf(stderr,"TFmdm_set_cosm(): Negative Hubble constant illegal.\n");
	  return -1;
      } else if (hubble>2.0) {
	  fprintf(stderr,"TFmdm_set_cosm(): Hubble constant should be in units of 100 km/s/Mpc.\n");
	  qwarn = 1;
      }
      if (redshift<=-1.0) {
	  fprintf(stderr,"TFmdm_set_cosm(): Redshift < -1 is illegal.\n");
	  return -1;
      } else if (redshift>99.0) {
	  fprintf(stderr,
	    "TFmdm_set_cosm(): Large redshift entered.  TF may be inaccurate.\n");
	  qwarn = 1;
      }
      if (degen_hdm<1) degen_hdm=1;
      num_degen_hdm = (double) degen_hdm;	
	  /* Have to save this for TFmdm_onek_mpc() */
      /* This routine would crash if baryons or neutrinos were zero, 
	  so don't allow that */
      if (omega_baryon<=0) omega_baryon=1e-5;
      if (omega_hdm<=0) omega_hdm=1e-5;

      omega_curv = 1.0-omega_matter-omega_lambda;
      omhh = omega_matter*SQR(hubble);
      obhh = omega_baryon*SQR(hubble);
      onhh = omega_hdm*SQR(hubble);
      f_baryon = omega_baryon/omega_matter;
      f_hdm = omega_hdm/omega_matter;
      f_cdm = 1.0-f_baryon-f_hdm;
      f_cb = f_cdm+f_baryon;
      f_bnu = f_baryon+f_hdm;

      /* Compute the equality scale. */
      z_equality = 25000.0*omhh/SQR(SQR(theta_cmb));	/* Actually 1+z_eq */
      k_equality = 0.0746*omhh/SQR(theta_cmb);

      /* Compute the drag epoch and sound horizon */
      z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
      z_drag_b2 = 0.238*pow(omhh,0.223);
      z_drag = 1291*pow(omhh,0.251)/(1.0+0.659*pow(omhh,0.828))*
		  (1.0+z_drag_b1*pow(obhh,z_drag_b2));
      y_drag = z_equality/(1.0+z_drag);

      sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1.0+10.0*pow(obhh,0.75));

      /* Set up for the free-streaming & infall growth function */
      p_c = 0.25*(5.0-sqrt(1+24.0*f_cdm));
      p_cb = 0.25*(5.0-sqrt(1+24.0*f_cb));

      omega_denom = omega_lambda+SQR(1.0+redshift)*(omega_curv+
			  omega_matter*(1.0+redshift));
      omega_lambda_z = omega_lambda/omega_denom;
      omega_matter_z = omega_matter*SQR(1.0+redshift)*(1.0+redshift)/omega_denom;
      growth_k0 = z_equality/(1.0+redshift)*2.5*omega_matter_z/
	      (pow(omega_matter_z,4.0/7.0)-omega_lambda_z+
	      (1.0+omega_matter_z/2.0)*(1.0+omega_lambda_z/70.0));
      growth_to_z0 = z_equality*2.5*omega_matter/(pow(omega_matter,4.0/7.0)
	      -omega_lambda + (1.0+omega_matter/2.0)*(1.0+omega_lambda/70.0));
      growth_to_z0 = growth_k0/growth_to_z0;	
      
      /* Compute small-scale suppression */
      alpha_nu = f_cdm/f_cb*(5.0-2.*(p_c+p_cb))/(5.-4.*p_cb)*
	  pow(1+y_drag,p_cb-p_c)*
	  (1+f_bnu*(-0.553+0.126*f_bnu*f_bnu))/
	  (1-0.193*sqrt(f_hdm*num_degen_hdm)+0.169*f_hdm*pow(num_degen_hdm,0.2))*
	  (1+(p_c-p_cb)/2*(1+1/(3.-4.*p_c)/(7.-4.*p_cb))/(1+y_drag));
      alpha_gamma = sqrt(alpha_nu);
      beta_c = 1/(1-0.949*f_bnu);
      /* Done setting scalar variables */
      hhubble = hubble;	/* Need to pass Hubble constant to TFmdm_onek_hmpc() */
      return qwarn;
  }

  /* ---------------------------- TFmdm_onek_mpc() ---------------------- */

  double TFmdm_onek_mpc(double kk)
  /* Given a wavenumber in Mpc^-1, return the transfer function for the
  cosmology held in the global variables. */
  /* Input: kk -- Wavenumber in Mpc^-1 */
  /* Output: The following are set as global variables:
	  growth_cb -- the transfer function for density-weighted
			  CDM + Baryon perturbations. 
	  growth_cbnu -- the transfer function for density-weighted
			  CDM + Baryon + Massive Neutrino perturbations. */
  /* The function returns growth_cb */
  {
      double	gamma_eff,	/* Effective \Gamma */
	    growth_cb,	/* Growth factor for CDM+Baryon perturbations */
	    growth_cbnu,	/* Growth factor for CDM+Baryon+Neutrino pert. */
	    max_fs_correction,  /* Correction near maximal free streaming */
	    qq,		/* Wavenumber rescaled by \Gamma */
	    qq_eff,		/* Wavenumber rescaled by effective Gamma */
	    qq_nu,		/* Wavenumber compared to maximal free streaming */
	    tf_master,	/* Master TF */
	    tf_sup,		/* Suppressed TF */
	    y_freestream; 	/* The epoch of free-streaming for a given scale */
    
    
      double tf_sup_L, tf_sup_C;
      double temp1, temp2;

      qq = kk/omhh*SQR(theta_cmb);

      /* Compute the scale-dependent growth functions */
      y_freestream = 17.2*f_hdm*(1+0.488*pow(f_hdm,-7.0/6.0))*
		  SQR(num_degen_hdm*qq/f_hdm);
      temp1 = pow(growth_k0, 1.0-p_cb);
      temp2 = pow(growth_k0/(1+y_freestream),0.7);
      growth_cb = pow(1.0+temp2, p_cb/0.7)*temp1;
      growth_cbnu = pow(pow(f_cb,0.7/p_cb)+temp2, p_cb/0.7)*temp1;

      /* Compute the master function */
      gamma_eff =omhh*(alpha_gamma+(1-alpha_gamma)/
		  (1+SQR(SQR(kk*sound_horizon_fit*0.43))));
      qq_eff = qq*omhh/gamma_eff;

      tf_sup_L = log(2.71828+1.84*beta_c*alpha_gamma*qq_eff);
      tf_sup_C = 14.4+325/(1+60.5*pow(qq_eff,1.11));
      tf_sup = tf_sup_L/(tf_sup_L+tf_sup_C*SQR(qq_eff));

      qq_nu = 3.92*qq*sqrt(num_degen_hdm/f_hdm);
      max_fs_correction = 1+1.2*pow(f_hdm,0.64)*pow(num_degen_hdm,0.3+0.6*f_hdm)/
		  (pow(qq_nu,-1.6)+pow(qq_nu,0.8));
      tf_master = tf_sup*max_fs_correction;

      /* Now compute the CDM+HDM+baryon transfer functions */
      double tf_cb = tf_master*growth_cb/growth_k0;
      //tf_cbnu = tf_master*growth_cbnu/growth_k0;
      return tf_cb;
  }

  /* ---------------------------- TFmdm_onek_hmpc() ---------------------- */

  double TFmdm_onek_hmpc(double kk)
  /* Given a wavenumber in h Mpc^-1, return the transfer function for the
  cosmology held in the global variables. */
  /* Input: kk -- Wavenumber in h Mpc^-1 */
  /* Output: The following are set as global variables:
	  growth_cb -- the transfer function for density-weighted
			  CDM + Baryon perturbations. 
	  growth_cbnu -- the transfer function for density-weighted
			  CDM + Baryon + Massive Neutrino perturbations. */
  /* The function returns growth_cb */
  {
      return TFmdm_onek_mpc(kk*hhubble);
  }

  double Dz0, A;
  bool initialized=false;

  double tophat(double k, double r)
  {
    double x = k*r;
    if (x <=1.e-6) return 1.;
    return 3.*(sin(x)-x*cos(x))/pow(x,3);
  }

  double tk_eisenstein_hu(double k)
  // Transfer function
  // see Eisenstein & Hu (1997)
  {
    return TFmdm_onek_hmpc(k);
  }


  double dz(double z)
  // growth factor w.r.t. z=0 (<1)
  // cf. Lahav & Suto 2003, Eqn. 67-70
  {
    if(!initialized) {
      double hz = OmegaM + OmegaL;
      double oMz0 = OmegaM/hz;
      double oLz0 = OmegaL/hz;
      Dz0 = 2.5*oMz0/(pow(oMz0,4./7.)-oLz0+(1.+oMz0/2.)*(1.+oLz0/70.)); 
    }
    double E2z = OmegaM*pow(1.+z,3)+OmegaL;
    double oMz = OmegaM*pow(1.+z,3)/E2z;
    double oLz = OmegaL/E2z;
    double Dz = 2.5*oMz/(1.+z)/(pow(oMz,4./7.)-oLz+(1.+oMz/2.)*(1.+oLz/70.));
    return Dz/Dz0; 
  }


  double pklin(double k, double z)
  // linear matter power spectrum at z
  {
    assert(initialized);
    //cout << A << " " << pow(k,ns) << " " << dz(z) << " " << tk(k) << endl;
    return A*pow(k,ns)*pow(dz(z)*tk_eisenstein_hu(k),2);
  }

  void initialize()
  {
    
    TFmdm_set_cosm(OmegaM, OmegaB, 0., 0, OmegaL, h, 0.);
    
    initialized=true;
    
    // set growth factor at present epoch
    double hz = OmegaM + OmegaL;
    double oMz = OmegaM/hz;
    double oLz = OmegaL/hz;
    Dz0 = 2.5*oMz/(pow(oMz,4./7.)-oLz+(1.+oMz/2.)*(1.+oLz/70.));
    
    // normalize power spectrum to sigma8
    double rsigma8 = 8.;
    double lnkmin = log(1.E-4);
    double lnkmax = log(40.);
    int nk = 1000;
    double dlnk = (lnkmax-lnkmin)/double(nk);
    double k=0.;
    double sum=0.;
    
    A=1.;
    
    for (int i=0; i<=nk; i++) {
      k = exp(lnkmin + i*dlnk);
      double mk = tophat(k,rsigma8);
      sum += pklin(k,0)*mk*mk*pow(k,3);
      //cout << pklin(k,0) << " " << mk << " " << k << endl;
    }
    sum *= dlnk/(2.*M_PI*M_PI);
    A = pow(sigma8,2)/sum;
    
    //cout << "A=" << A << " " << sum << " " << sigma8 << endl;
    
  }

}

namespace xilinear {
 static double *logr=0;
 static double *logxi=0;
 
 static double alphasmall; // power-law slope at small r
 static double alphalarge; // power-law slope at large r
 
 const double logrstep=0.01;
 const double logrmin=-3;
 const int nsteps=501;
 
 void initialize()
 {
  std::ifstream in("lut/2pc.tab");
  assert(in.is_open());
  
  // check cosmology
  double buf;
  in>>buf;
  assert(OmegaM==buf);
  in>>buf;
  assert(OmegaL==buf);
  in>>buf;
  assert(OmegaB==buf);
  in>>buf;
  assert(sigma8==buf);
  in>>buf;
  assert(ns==buf);
   
  // read
  
  logr = new double[nsteps];
  logxi = new double[nsteps];
  
  for(int i=0; i<nsteps; i++) {
    assert(in.eof()==0);
    in >> logr[i] >> logxi[i];
    assert(logxi[i]>0);
    logxi[i] = log10(logxi[i]);
  }
  in>>buf;
  assert(in.eof());
  
  alphasmall=(logxi[5]-logxi[0])/5.;
  alphalarge=(logxi[nsteps-1]-logxi[nsteps-6])/5.;
  
  assert(alphasmall<0);
  assert(alphalarge<0);
  
  pklinear::initialize();
  
 }
 
 double twopclin(double log10r, double z)
 // log10(r/[h^-1 Mpc] comoving), z
 {
   assert(log10r>=-4);
   assert(log10r<=3);
   
   if(!logr) initialize();
  
   double n=(log10r-logrmin)/logrstep;
   if(n>0 && n<nsteps-1) { // linear interpolation in log r
     return pow10(logxi[int(n)]*(1.-n+int(n))+logxi[int(n)+1]*(n-int(n)))*pow(pklinear::dz(z),2);
   }
   
   else if(n<=0)
     return pow10(logxi[0]+alphasmall*n)*pow(pklinear::dz(z),2);
   //else if(n>=nsteps-1)
     return pow10(logxi[nsteps-1]+alphalarge*(n-nsteps+1))*pow(pklinear::dz(z),2);
     
}
  
}

namespace Wlinear {
  
 static double *logtheta=0;
 static double *logW=0;
 static double zW=-1;
 
 static double alphalarge; // power-law slope at large r
 
 const double logthetastep=log10(1.01);
 const double logthetamin=-2;
 const int nsteps=875;
 
 void initialize(double z)
 {
  zW=z;
  std::ostringstream ss;
  ss << int(z*10000000.+0.5);
  std::ifstream in(("lut/W0."+ss.str()+".tab").c_str());
  assert(in.is_open());
 

  const double eps=1.e-7; 
  // check cosmology
  double buf;
  in>>buf;
  assert(fabs(OmegaM-buf)<eps);
  in>>buf;
  assert(fabs(OmegaL-buf)<eps);
  in>>buf;
  assert(fabs(OmegaB-buf)<eps);
  in>>buf;
  assert(fabs(sigma8-buf)<eps);
  in>>buf;
  assert(fabs(ns-buf)<eps);
  in>>buf;
  assert(fabs(z-buf)<eps);
 
  // read
  
  if(logtheta) delete []logtheta;
  if(logW) delete []logW;
  
  logtheta = new double[nsteps];
  logW = new double[nsteps];
  
  for(int i=0; i<nsteps; i++) {
    assert(in.eof()==0);
    in >> logtheta[i] >> logW[i];
  }
  in>>buf;
  assert(in.eof());
  
  alphalarge=(logW[nsteps-1]-logW[nsteps-6])/5.;
  
  assert(alphalarge<0);  
 }
 
 double Wlog(double log10theta, double z)
 // log10(theta/arcmin), z
 // output: h^{-3} Mpc^3 per arcmin^2
 {
   assert(log10theta<=3);
   
   if(!logtheta) initialize(z);
   if(z!=zW) initialize(z);
  
   double n=(log10theta-logthetamin)/logthetastep;
   if(n>0 && n<nsteps-1) { // linear interpolation in log theta
     return pow10(logW[int(n)]*(1.-n+int(n))+logW[int(n)+1]*(n-int(n)));
   }
   
   else if(n<=0)
     return pow10(logW[0]); // constant at small theta
   //else if(n>=nsteps-1) 
     return pow10(logW[nsteps-1]+alphalarge*(n-nsteps+1.)); // power law at large theta
 }
  
}


namespace halomodel {

 static double *logM_sig=0, *logsig=0;
 const double logMmin=7.;
 const double logMmax=16.;
 const int nsigsteps=181;
 const double logMsigstep=(logMmax-logMmin)/double(nsigsteps-1);
 static bool initialized_sigmam=false;
 
 void initialize_sigmam()
 {
  std::ifstream in("lut/sigmam.tab");
  assert(in.is_open());
  
  // check cosmology
  double buf;
  in>>buf;
  assert(OmegaM==buf);
  in>>buf;
  assert(OmegaL==buf);
  in>>buf;
  assert(OmegaB==buf);
  in>>buf;
  assert(sigma8==buf);
  in>>buf;
  assert(ns==buf);
   
  // read
  
  if(logM_sig) delete logM_sig;
  if(logsig) delete logsig;
  
  logM_sig = new double[nsigsteps];
  logsig = new double[nsigsteps];
  
  for(int i=0; i<nsigsteps; i++) {
    assert(in.eof()==0);
    in >>  buf >> logM_sig[i] >> logsig[i];
    logsig[i] = log10(logsig[i]);
  }
  in>>buf;
  assert(in.eof());  
  
  initialized_sigmam=true;
 }
 
 double sigmam(double M200m, double z) // hinv Msol
 {   
   double logM = log10(M200m);
   //std::cout << M200m << std::endl;
   assert(logM>=logMmin);
   assert(logM<=logMmax);
  
   if(!initialized_sigmam) initialize_sigmam();
   
   double n=(logM-logMmin)/logMsigstep;
   assert(n>=0.);
   assert(n<=nsigsteps-1);
   
   if(n>0 && n<nsigsteps-1) { // linear interpolation in log-log space
     return pow10(logsig[int(n)]*(1.-n+int(n))+logsig[int(n)+1]*(n-int(n)))*pklinear::dz(z);
   }
   
   else if(n==0)
     return pow10(logsig[0])*pklinear::dz(z);
   else 
     return pow10(logsig[nsigsteps-1])*pklinear::dz(z);
 }
  
  // Tinker et al. 2011
  
  static double
   alpha=0.368,
   beta0=0.589,
   gamma0=0.864,
   phi0=-0.729,
   eta0=-0.243,
   deltac=1.686;
   
  static double rhom0=OmegaM*rho0crit;

  double nu(double M200m, double z) // Bias from Tinker's mass function, Delta=200 w.r.t. mean matter density; hinv Msol
  {
    return deltac/sigmam(M200m,z);
  }
    
  double bias_tinker(double M200m, double z) // Bias from Tinker's mass function, Delta=200 w.r.t. mean matter density; hinv Msol
  {       
    double beta=beta0*pow(1.0+z,0.20);
    double phi=phi0*pow(1.0+z,-0.08);
    double eta=eta0*pow(1.0+z,0.27);
    double gamma=gamma0*pow(1.0+z,-0.01);
      
    double NU=nu(M200m, z);
    
    //std::cout << gamma*NU*NU << (1.0+2.0*eta) << deltac << " "  << 2.0*phi/deltac/(1.0+pow(beta*NU,2.0*phi)) << std::endl;
    
    return 1.0+(gamma*NU*NU-(1.0+2.0*eta))/deltac+2.0*phi/deltac/(1.0+pow(beta*NU,2.0*phi));
  }

  static double *logM_b=0,*logM_dndM=0,*logB=0,*logdndM=0;
  static double zh=-1;
  const int nbsteps=901;
  const int ndndmsteps=200;
  const double logMbstep=(logMmax-logMmin)/double(nbsteps-1);
  const double logMdndmstep=(logMmax-logMmin)/double(ndndmsteps-1);
  static bool initialized=false;

  void initialize_bias()
  {
    std::ostringstream ss;
    ss << int(zh*10000000.+0.5);
    std::ifstream in(("lut/bias0."+ss.str()+".tab").c_str());
    assert(in.is_open());
    
    // check cosmology
    double buf;
    in>>buf;
    assert(OmegaM==buf);
    in>>buf;
    assert(OmegaL==buf);
    in>>buf;
    assert(OmegaB==buf);
    in>>buf;
    assert(sigma8==buf);
    in>>buf;
    assert(ns==buf);
    in>>buf;
    assert(zh==buf);
  
    // read
    
    if(logM_b) delete []logM_b;
    if(logB) delete []logB;
    
    logM_b = new double[nbsteps];
    logB = new double[nbsteps];
    
    for(int i=0; i<nbsteps; i++) {
      assert(in.eof()==0);
      in >> logM_b[i] >> buf >> logB[i];
      assert(logB[i]>0);
      logB[i]=log10(logB[i]);
    }
    in>>buf;
    assert(in.eof());
  }

  void initialize_dndM()
  {
    std::ostringstream ss;
    ss << int(zh*10000000.+0.5);
    std::ifstream in(("lut/dndM0."+ss.str()+".tab").c_str());
    assert(in.is_open());
    
    // check cosmology
    double buf;
    in>>buf;
    assert(OmegaM==buf);
    in>>buf;
    assert(OmegaL==buf);
    in>>buf;
    assert(OmegaB==buf);
    in>>buf;
    assert(sigma8==buf);
    in>>buf;
    assert(ns==buf);
    in>>buf;
    assert(zh==buf);
  
    // read
    
    if(logM_dndM) delete []logM_dndM;
    if(logdndM) delete []logdndM;
    
    logM_dndM = new double[ndndmsteps];
    logdndM = new double[ndndmsteps];
    
    for(int i=0; i<ndndmsteps; i++) {
      assert(in.eof()==0);
      in >> logM_dndM[i] >> logdndM[i];
      assert(logM_dndM[i]>0);
      logM_dndM[i]=log10(logM_dndM[i]);
      assert(logdndM[i]>0);
      logdndM[i]=log10(logdndM[i]);
    }
    in>>buf;
    assert(in.eof());
  }

  void initialize(double z)
  {
    zh=z;
    
    initialize_bias();
    initialize_dndM();
    
    initialized=true;
  }

  double dndM(double M200m, double z)
  {  
    if(!initialized || z!=zh) {zh=z; initialize_dndM();}
    
    double logM = log10(M200m);
    assert(logM>=logMmin);
    assert(logM<=logMmax);
    
    double n=(logM-logMmin)/logMdndmstep;
    assert(n>=0.);
    assert(n<=ndndmsteps-1);
    
    if(n>0 && n<ndndmsteps-1) { // linear interpolation in log-log space
      return pow10(logdndM[int(n)]*(1.-n+int(n))+logdndM[int(n)+1]*(n-int(n)));
    }
    
    else if(n==0)
      return pow10(logdndM[0]);
    else 
      return pow10(logdndM[ndndmsteps-1]);
  }

  double bias(double M200m, double z)
  {  
    if(!initialized || z!=zh) initialize(z);
    
    double logM = log10(M200m);
    assert(logM>=logMmin);
    assert(logM<=logMmax);
    
    double n=(logM-logMmin)/logMbstep;
    assert(n>=0.);
    assert(n<=nbsteps-1);
    
    if(n>0 && n<nbsteps-1) { // linear interpolation in log-log space
      return pow10(logB[int(n)]*(1.-n+int(n))+logB[int(n)+1]*(n-int(n)));
    }
    
    else if(n==0)
      return pow10(logB[0]);
    else 
      return pow10(logB[nbsteps-1]);
  }
  
  
  double cmz_200m_duffy(double m200m,double z)
  // returns concentration w.r.t. 200xmean matter density according to Duffy et al. (2008) simulations
  // this is true at WMAP5 cosmology, 
  // m200m should be in h^{-1} Msol units, since the pivot mass 2e12 is
  {
    return 10.14*pow(m200m/2.e12,-0.081)/pow(1.+z,1.01);
  }
  
  
  double r200mRadius(double m200m, double z)
  // returns r_200m of halo [proper Mpc, i.e. physical radius; NO hinv]
  // z: redshift of the halo
  // m200m: mass at radius at which average density is equal to 200*rho_mean,matter, NO hinv
  {
	  return pow(m200m*Gs/(100.*OmegaM*H0*H0),1./3.)/(1.+z); // proper radius [Mpc]!
  }
  
  double rDeltaMRadius(double m, double DeltaM, double z)
  // returns r_[DeltaM]m of halo [proper Mpc, i.e. physical radius; NO hinv]
  // z: redshift of the halo
  // m: mass at radius at which average density is equal to DeltaM*rho_mean,matter, NO hinv
  {  
	  return pow(2.*m*Gs/(DeltaM*OmegaM*H0*H0),1./3.)/(1.+z); // proper radius [Mpc]!
  }
  
  double tau200m(double nu)
  {
    // caveat: fitted at z=0.24533
    return std::max(1.15,3.85-0.73*nu); 
  }
  
  double tau200m(double m200m, double z)
  {
    assert(z==0.24533);
    return tau200m(nu(m200m,z));
  }
}


namespace Pkappa {
  
 static double *logl=0;
 static double *logPkappa=0;
 
 static double alphalarge; // power-law slope at large l
 static double alphasmall; // power-law slope at large l
 
 const double loglmin= -2.;
 const double loglmax=  5.;
 const int nsteps=1002;
 const double loglstep=(loglmax-loglmin)/double(nsteps-1);
 
 static bool initialized=false;
 
 void initialize (std::string filename="lut/Pkappa.tab")
 {

  std::ifstream in(filename.c_str());
  assert(in.is_open());
  
  // read
  
  if(logl) delete []logl;
  if(logPkappa) delete []logPkappa;
  
  logl = new double[nsteps];
  logPkappa = new double[nsteps];
  
  for(int i=0; i<nsteps; i++) {
    assert(in.eof()==0);
    in >> logl[i] >> logPkappa[i];
    assert(logl[i]>0);
    assert(logPkappa[i]>0);
    logl[i] = log10(logl[i]);
    logPkappa[i] = log10(logPkappa[i]);
  }
  std::string buf;
  in>>buf;
  assert(in.eof());
  
  alphalarge=(logPkappa[nsteps-1]-logPkappa[nsteps-6])/5.;
  alphasmall=(logPkappa[5]-logPkappa[0])/5.;
  
  assert(alphalarge<0);
  assert(alphasmall>0);
  
  initialized=true;
 }
 
 double logPkappalog(double log10l)
 // argument and return value are log10
 {
   assert(initialized);
     
   double n=(log10l-loglmin)/loglstep;
   if(n>0 && n<nsteps-1) { // linear interpolation in log theta
     return logPkappa[int(n)]*(1.-n+int(n))+logPkappa[int(n)+1]*(n-int(n));
   }
   
   // power-law extrapolate to really large and really small l
   else if(n<=0)
     return logPkappa[0]+alphasmall*n;
   //else if(n>=nsteps-1) 
     return logPkappa[nsteps-1]+alphalarge*(n-nsteps+1.);
 }
  
}


#endif
