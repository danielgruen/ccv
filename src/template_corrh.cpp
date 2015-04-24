// calculate covariance matrix of kappa due to correlated haloes of fixed mass

#include <CCfits/CCfits>
using namespace CCfits;

#include "corrh/corrh.h"
#include "corrh/template_corrh.h"
#include "profile/profile.h"
#include "enfw/enfw.h"
#include <iostream>
#include <cstdlib>
#include <string.h>
//#include <omp.h>

using namespace std;
using namespace halosigma;

double rmin[rsteps];       // inner radius, Mpc
double rmax[rsteps];       // outer radius, Mpc
double rmin2[rsteps];       // inner radius^2, Mpc^2
double rmax2[rsteps];       // outer radius^2, Mpc^2
double area[rsteps];       // area, Mpc^2

const double zeros[rsteps] = { };

void clear_cov(double cov[rsteps][rsteps])
{
 // clear cov matrix
 for(int i=0; i<rsteps; i++)
   memcpy(cov[i], zeros, rsteps*sizeof(double));
}

void halo_cov(double dens[rsteps], double cov[rsteps][rsteps],  double kahan[rsteps][rsteps], double dn, int minannulus, int maxannulus)
// fills only upper triangle cov[i][j], j>=i
// Kahan summation
{
 
 for(int i=minannulus; i<=maxannulus; i++)
 {
   double densi = dens[i]*dn;
   double *d=dens+i;
   double *c=cov[i]+i;
   double *k=kahan[i]+i;
   
   for(int j=i; j<=maxannulus; j++)
   {
     double y = densi*(*d) - (*k);
     double t = (*c) + y;
     (*k) = ( t-(*c) ) - y; // not zero due to numerical inaccuracy
     //if ((*k) != 0) // happens frequently!
     //{
     // cout << "!!! " << (*k) << endl;
     //}
     (*c) = t;
     
     c++; d++; k++;
   }
 }

}

void halo_cov_fill(double cov[][rsteps])
{
  // unnecessary: copy to lower triangle
  for(int i=1; i<rsteps; i++)
  {
   for(int j=i-1; j>=0; j--)
   {
     cov[i][j]=cov[j][i];  
   }
  }
}



int main(int argc, char **argv)
{
  if(argc!=4) {
   cerr << "syntax: " << argv[0] << " [z] [log M200m in h^-1 Msol] [filename output]" << endl;
   return 1;
  }
  double z=atof(argv[1]);
  Wlinear::initialize(z);
  double M=pow10(atof(argv[2]));
  
  int nthreads=1;  

  assert(M>0.99e8);
  assert(M<1.01e16);

  double Da=angularDiameterDistance(0,z,10000)*radperarcmin; // Mpc / arcmin
    
  const double c200m=halomodel::cmz_200m_duffy(M,z);
  const double rs=halomodel::r200mRadius(M/h,z)/c200m*h;   // proper h^-1 Mpc
  
  const double nu   =  halomodel::nu(M,z);
  double tau  =  c200m*(3.85-0.73*nu);
  if(nu>3.7) tau = c200m*(3.85-0.73*3.7); 
  const double hardt=  100.;
  
  const double rho0=bmo_rho0_200m(c200m, tau, z)/h/h;      // h^2 Msol per proper Mpc^3
  const double A=rho0*rs; // h Msol per proper Mpc^2
  
  cerr << "# halo of mass M200m/Msol=" << M << " c200m=" << c200m << " rs/Mpc=" << rs << " rho0/Msol*Mpc^3=" << rho0 << " A/Msol*Mpc^2=" << A << endl;
  cerr << "# nu=" << nu << " tau=" << tau << endl;
  
  rmin[0]=thetamin*Da;
  rmax[0]=rmin[0]*rfactor;
  area[0]=M_PI*(rmax[0]*rmax[0]-rmin[0]*rmin[0]);
  rmin2[0]=rmin[0]*rmin[0];
  rmax2[0]=rmax[0]*rmax[0];
  
  for(int i=1; i<rsteps; i++)
  {
    rmin[i]=rmax[i-1];
    rmax[i]=rmin[i]*rfactor;
    rmin2[i]=rmin[i]*rmin[i];
    rmax2[i]=rmax[i]*rmax[i];
    area[i]=M_PI*(rmax[i]*rmax[i]-rmin[i]*rmin[i]);
  }
  
  double thetafac  = pow(rfactor,0.01);
  double theta2fac = M_PI*(pow(thetafac,2)-pow(1./thetafac,2))/2.;
  
  double dens[nthreads][rsteps];
  double cov[nthreads][rsteps][rsteps];   // covariance of surface density, [Msol/Mpc^2]^2
  double kahan[nthreads][rsteps][rsteps]; // buffer for Kahan summation
  
  for(int i=0; i<nthreads; i++)
  {
    clear_cov(cov[i]);
    clear_cov(kahan[i]);
  }
  
  
  int minannulus0, maxannulus0;
  
  // contribution of haloes in the (projected) center,
  // inside thetamin/5.
  
  halosigma::halosigma(thetamin/10.*Da, A, rs, rsteps, rfactor, dens[0], rmin, rmax, rmin2, rmax2, area, zeros, minannulus0, maxannulus0, tau, hardt);
  
  halo_cov(dens[0],cov[0],kahan[0],Wlinear::Wlog(-3,z)*M_PI* ( pow(thetamin/thetafac/5.,2) + 0.5*(pow(thetamin*thetafac/5.,2)-pow(thetamin/5.,2))), minannulus0, maxannulus0);
  
  cout << "# integrating from " << thetamin/5. << " to " << thetamax+rs*tau/Da << " arcmin" << endl;
  
  double lthetamin=log10(thetamin/5.);
  double lthetamax=log10(thetamax+rs*tau/Da);
  double lthetafac=log10(thetafac);
  int ithetamax=(lthetamax-lthetamin)/lthetafac+1;
  
//#pragma omp parallel for schedule(dynamic)
  for(int itheta=0; itheta<ithetamax; itheta++)
  { 
    int threadid = 0; // omp_get_thread_num();
    
    double theta=pow10(lthetamin+itheta*lthetafac);
    if(log10(theta)>3) continue;
    double dtheta2=theta2fac*theta*theta;
    
    int minannulus, maxannulus;
    
//   const double r,	// separation of halo from centre [Mpc]
//   const double A, 	// amplitude to be multiplied with BMO profile [Msol/Mpc^2]
//   const double rs, 	// scale radius [Mpc]
//   const int rsteps, 	// number of sky annuli
//   const double rfactor, // radius factor beween subsequent sky annuli
//   double *dens, 	// density in sky annuli; array will be overwritten with result of the calculation [Msol/Mpc^2]
//   const double *rmin, 	// inner radius of each sky annulus [Mpc]
//   const double *rmax, 	// outer radius of each sky annulus
//   const double *rmin2, 	// square of inner radius of each sky annulus [Mpc^2]
//   const double *rmax2,	// square of outer radius of each sky annulus
//   const double *area,	// area of each sky annulus [Mpc^2]
//   const double *zeros,	// array of zeroes, copied into dens before filling it with actual density
//   int &nonzero_minannulus,
//   int &nonzero_maxannulus)
    halosigma::halosigma(theta*Da, A, rs, rsteps, rfactor, dens[threadid], rmin, rmax, rmin2, rmax2, area, zeros, minannulus, maxannulus, tau, hardt); 
      // dens: h Msol per proper Mpc^2
    
    int buf;
    
    halo_cov(dens[threadid], cov[threadid], kahan[threadid], Wlinear::Wlog(log10(theta),z)*dtheta2, minannulus, maxannulus);
      // cov: h^2 Msol^2 per proper Mpc^4 if there were 1 such halo per com. hinv^3 Mpc^3 Lagrangian volume

    cerr << itheta << " / " << ithetamax << "          \r" << flush;
  }
  
  // collect
//#pragma omp parallel for
  for(int i=0; i<rsteps; i++)
  {
   for(int j=i; j<rsteps; j++)
   {
    for(int n=1; n<nthreads; n++)
    {
      cov[0][i][j] += cov[n][i][j];
    }
    cov[0][j][i] = cov[0][i][j];
   }
  }
  
  
  FITS *pFits = 0; 
  long naxis    =   2;      
  long naxes[2] = { rsteps, rsteps };

  try
  {
    const std::string fileName(argv[3]);            
    pFits = new FITS(fileName, FLOAT_IMG , naxis , naxes );
  }
  catch (FITS::CantCreate)
  {
    cerr << "error, can't create fits file" << endl;
    return 1;
  } 

  valarray <float> image(0., rsteps*rsteps); // we're just going to jam our data in here ...
  int idx=0;
  for(int i=0; i<rsteps; i++)
    for(int j=0; j<rsteps; j++)
    {
     image[idx]=float(cov[0][i][j]);
     idx++;
    }
  
  pFits->pHDU().write(1,rsteps*rsteps,image);
  
  return 0;
}

