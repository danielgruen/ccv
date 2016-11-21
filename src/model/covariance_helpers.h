#ifndef _COVARIANCE_HELPERS_H_
#define _COVARIANCE_HELPERS_H_

const string ccv_basedir="/Users/fabrice/data/ccv/";
const double ckms=c/1000.;
const double Gs=Gnewton*(1.9891E30)/(3.08567758E22)*1.E-6; // = 4.25e-9;  // Mpc/Msol (km/s)^2
const double h=H0/100.;
const double OmegaM=omegaM;
const double OmegaL=omegaL;
const double OmegaB=0.044;
const double sigma8=0.79;
const double ns=0.95;

double rho0crit = 3.*H0*H0/(8.*M_PI*Gs); // [ Msol / Mpc^3 ]

const double arcminperrad=180.*60./M_PI;
const double radperarcmin=1./arcminperrad;

double pow10(const double a) {
  return pow(10.,a);
}

double pow3(const double a) { return a*a*a; }

namespace pklinear {

  // taken from http://background.uchicago.edu/~whu/transfer/power.c

  double dz(double z)
  // growth factor w.r.t. z=0 (<1)
  // cf. Lahav & Suto 2003, Eqn. 67-70
  {
    double hz = OmegaM + OmegaL;
    double oMz0 = OmegaM/hz;
    double oLz0 = OmegaL/hz;
    double Dz0 = 2.5*oMz0/(pow(oMz0,4./7.)-oLz0+(1.+oMz0/2.)*(1.+oLz0/70.)); 
    double E2z = OmegaM*pow(1.+z,3)+OmegaL;
    double oMz = OmegaM*pow(1.+z,3)/E2z;
    double oLz = OmegaL/E2z;
    double Dz = 2.5*oMz/(1.+z)/(pow(oMz,4./7.)-oLz+(1.+oMz/2.)*(1.+oLz/70.));
    return Dz/Dz0; 
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
  std::ifstream in((ccv_basedir+string("lut/sigmam.tab")).c_str());
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
  
  double r200mRadius(double m200m, double z)
  // returns r_200m of halo [proper Mpc, i.e. physical radius; NO hinv]
  // z: redshift of the halo
  // m200m: mass at radius at which average density is equal to 200*rho_mean,matter, NO hinv
  {
	  return pow(m200m*Gs/(100.*OmegaM*H0*H0),1./3.)/(1.+z); // proper radius [Mpc]!
  }
    
  double cmz_200m_duffy(double m200m,double z)
  // returns concentration w.r.t. 200xmean matter density according to Duffy et al. (2008) simulations
  // this is true at WMAP5 cosmology, 
  // m200m should be in h^{-1} Msol units, since the pivot mass 2e12 is
  {
    return 10.14*pow(m200m/2.e12,-0.081)/pow(1.+z,1.01);
  }
  

  static double *logM_b=0,*logM_dndM=0,*logB=0,*logdndM=0;
  static double zh=-1;
  const int nbsteps=901;
  const int ndndmsteps=200;
  const double logMbstep=(logMmax-logMmin)/double(nbsteps-1);
  const double logMdndmstep=(logMmax-logMmin)/double(ndndmsteps-1);
  static bool initialized=false;
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


#endif

