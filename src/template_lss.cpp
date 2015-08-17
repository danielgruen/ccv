#include <iostream>
using namespace std;
#include <assert.h>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <CCfits/CCfits>
using namespace CCfits;
#include <omp.h>
#include "cosmology.h"
#include "corrh/template_corrh.h"

const double pi4 = M_PI/4.;
const double pi34 = 3.*M_PI/4.;

double pow2(double a) { return a*a; }
double pow3(double a) { return a*a*a; }

double J0(double x)
{ return gsl_sf_bessel_J0 (x); }

double J0hat(double x, double delta) // J0hat: area-weighted average of J0 in annulus x(1-delta)...x(1+delta)
{ return ((1.+delta)*gsl_sf_bessel_J1(x*(1.+delta))-(1.-delta)*gsl_sf_bessel_J1(x*(1.-delta)) )/(2.*delta*x); }

double J1(double x)
{ return gsl_sf_bessel_J1(x); }

double J0_small(double x)
{ return 1.-x*x/4.; }

double J0_notsosmall(double x)
{ double x24=x*x/4.; return 1.-x24*(1.-x24/16.); }

double J1_small(double x)
{ return x/2.; }

double J0_large(double x)
{ return sqrt(2./(M_PI*x))*cos(x-pi4); }

double J1_large(double x)
{ return sqrt(2./(M_PI*x))*cos(x-pi34);}


const double lmin=1.e-5;          // lower limit for l integration
const double lmax=1.e7;           // upper limit for l integration
const double loglmin=log10(lmin);
const int nsinsteps=360;              // this many integrations steps for one period in kr of the sin(kr)*k^alpha or this many integration steps for a k interval, whichever is smaller
////// WAS 2000!!!
const double radmax=2.*M_PI/nsinsteps; // max d(angle [rad]) for integration of trig fn
//const double lfac=1.01;              // multiplicative factor in l in first integral
const double lfac=1.02;              // multiplicative factor in l in first integral
const double lnlfacinv=log(1./lfac);


double cov[rsteps][rsteps] = { }; // covariance of surface density, [Msol/Mpc^2]^2

const double zeros[rsteps] = { };

void clear_cov()
{
 // clear cov matrix
 for(int i=0; i<rsteps; i++)
   memcpy(cov[i], zeros, rsteps*sizeof(double));
}

double plsi(double Pkk0, double k0, double alpha, double r, double k1)
// integrate k0*P(k0)*(k/k0)^(alpha)*sin(k*r) from k0 to k1
// good idea for (k1-k0)*r< 4 PI since we're doing nsinsteps steps
// Pkk0 = k0*P(k0)
{

  double kstep=(k1-k0)/double(nsinsteps);
  
  double integrandprev=Pkk0*sin(k0*r);
  
  Pkk0 /= pow(k0,alpha); // good for later
  
  double s=0.;
  
  for (int ik=1; ik<nsinsteps; ik++)
  {
   double kmax=k0+double(ik+1)*kstep;
   
   double integrandnext=Pkk0*pow(kmax,alpha)*sin(kmax*r);
   
   s += kstep*(integrandprev+integrandnext)/2.;
   
   integrandprev = integrandnext;
  }
  
  return s;
}

double plsi_periods(double Pkk0, double k0, double alpha, double r, int nperiods)
// integrate k0*P(k)*(k/k0)^(alpha)*sin(k*r) from k0 to k1=k0+2 PI n/r
// good idea for (k1-k0)*r>= 4 PI
{
  assert(nperiods>1);
  double s=plsi(Pkk0,k0,alpha,r,k0+2.*M_PI/r); // s is the integral over the first period
  
  // result is
  //      n   /        2 PI i  \ (alpha)
  // s * sum  |  1  + -------- |
  //     i=0  \         k0 r   /
  // which I fear has no simple form, so we just add them all up
  
  double t=0.;
  double A=2.*M_PI/(k0*r);
  for(int n=0; n<nperiods; n++)
  {
    t += pow(1.+A*n,alpha);
  }
  return t*s;
}



double plci(double Pkk0, double k0, double alpha, double r, double k1)
// integrate k0*P(k0)*(k/k0)^(alpha)*cos(k*r) from k0 to k1
// good idea for (k1-k0)*r< 4 PI since we're doing nsinsteps steps
// Pkk0 = k0*P(k0)
{

  double kstep=(k1-k0)/double(nsinsteps);
  
  double integrandprev=Pkk0*cos(k0*r);
  
  Pkk0 /= pow(k0,alpha); // good for later
  
  double s=0.;
  
  for (int ik=1; ik<nsinsteps; ik++)
  {
   double kmax=k0+double(ik+1)*kstep;
   
   double integrandnext=Pkk0*pow(kmax,alpha)*cos(kmax*r);
   
   s += kstep*(integrandprev+integrandnext)/2.;
   
   integrandprev = integrandnext;
  }
  
  return s;
}

double plci_periods(double Pkk0, double k0, double alpha, double r, int nperiods)
// integrate k0*P(k)*(k/k0)^(alpha)*cos(k*r) from k0 to k1=k0+2 PI n/r
// good idea for (k1-k0)*r>= 4 PI
{
  assert(nperiods>1);
  double s=plci(Pkk0,k0,alpha,r,k0+2.*M_PI/r); // s is the integral over the first period
  
  // result is
  //      n   /        2 PI i  \ (alpha)
  // s * sum  |  1  + -------- |
  //     i=0  \         k0 r   /
  // which I fear has no simple form, so we just add them all up
  
  double t=0.;
  double A=2.*M_PI/(k0*r);
  for(int n=0; n<nperiods; n++)
  {
    t += pow(1.+A*n,alpha);
  }
  return t*s;
}



int main(int argc, char **argv)
{
    if(argc!=3)
    {
     cout << "syntax: " << argv[0] << " [Pkappa file] [FITS output]" << endl; 
     return 1; 
    }
    
    clear_cov();
    Pkappa::initialize(argv[1]);

#pragma omp parallel for schedule(dynamic)
    for(int ilta=0; ilta<rsteps; ilta++)
    for(int iltb=ilta; iltb<rsteps; iltb++)
    {
      const double thetaa = thetamin*pow10(log10step*double(ilta));
      const double thetab = thetamin*pow10(log10step*double(iltb));
      
      const double theta1 = min(thetaa,thetab)*radperarcmin;
      const double theta2 = max(thetaa,thetab)*radperarcmin;
      const double t1pt2  = theta1+theta2;
      const double t2mt1  = theta2-theta1;
      const double t12pt22= theta1*theta1+theta2*theta2;
      
      assert(theta2>=theta1);
      assert(theta1>0);
      assert(theta2>0);
    
      double l;
      double sum=0.;
      double lprev=lmin;
      
      // first integration: two terms, one has an extra factor l*l*(theta1^2+theta2^2)/4
      const double lmax1=0.2/theta2;

      double sum_a=0.;
      double sum_b=0.;
      double integrandprev_a = pow10(Pkappa::logPkappalog(loglmin))*lmin; 
      double integrandprev_b = -integrandprev_a*lmin*lmin*t12pt22/4.;
      
      for(l=lmin*lfac; l<lmax1; l*=lfac)
      {
      
	double integrand_a = pow10(Pkappa::logPkappalog(log10(l)))*l;
	double integrand_b = -integrand_a*l*l*t12pt22/4.;

	double alpha_a = log(integrandprev_a/integrand_a)/lnlfacinv;
	double alpha_b = log(integrandprev_b/integrand_b)/lnlfacinv;
	
	sum_a += integrandprev_a*lprev/(alpha_a+1.)*(pow(lfac,alpha_a+1.)-1.);
	sum_b += integrandprev_b*lprev/(alpha_b+1.)*(pow(lfac,alpha_b+1.)-1.);
	
	
	lprev = l;
	integrandprev_a = integrand_a;
	integrandprev_b = integrand_b;
      }
      
      sum += sum_a;
      sum += sum_b;
      
      double integrandprev = integrandprev_a+integrandprev_b;
      
      // second integration: pain-in-the ass trapezoid
      
      const double lmax2 = 50./theta1;
      
      for(; l<lmax2; )
      {
	double lstep=min((lfac-1.)*l,radmax/(t1pt2)); // dl*(theta1+theta2) must be << 2PI; also, dl/l should be small
	l += lstep;
	
	double integrand = l*J0(l*theta1)*J0(l*theta2)*pow10(Pkappa::logPkappalog(log10(l)));
	
	sum += (integrandprev+integrand)/2.*lstep;
	
	integrandprev = integrand;
      }
      
      sum /= (2.*M_PI); // had left out this factor in the first two integrals
      
      // third integration
      // a) sin part

      lprev = l;
      double Plprev = pow10(Pkappa::logPkappalog(log10(l)));
      sum_a=0.;

      for(double la=l*lfac; la<lmax; la*=lfac)
      {
	bool multipleperiods=false;
	int nperiods=(la-lprev)*t1pt2/(2.*M_PI);
	
	if(nperiods>=2) 
	// at least two full periods are inside our current l interval 
	// change integration strategy to do a fine integration of one 2pi period 
	// and then sum up all periods, re-scaled according to the power-law in l
	{
	  multipleperiods=true;
	
	  la=lprev+2.*M_PI*nperiods/t1pt2; // make sure knext is an integer multiple of 2pi periods in l(t1+t2) away from lprev
	}
	
	// P(l) is a power-law between lprev and la:
	// P(l) = Plprev * (la/lprev)^alpha
	// alpha = log(Plprev/Plnext) / log(la/lprev)
	
	double Plnext=pow10(Pkappa::logPkappalog(log10(la)));
	
	double alpha=log10(Plnext/Plprev)/log10(la/lprev);
	
	      
	if(multipleperiods)
	  sum_a += plsi_periods(Plprev,lprev,alpha,t1pt2,nperiods);
	else
	  sum_a += plsi(Plprev,lprev,alpha,t1pt2,la);
	
	lprev  = la;
	Plprev = Plnext;
	
      }
      
      // b) cos part
      
      lprev = l;
      Plprev = pow10(Pkappa::logPkappalog(log10(l)));
      sum_b=0.;

      if(t2mt1>0) {
      
	for(double lb=l*lfac; lb<lmax; lb*=lfac)
	{
	  bool multipleperiods=false;
	  int nperiods=(lb-lprev)*t2mt1/(2.*M_PI);
	  
	  if(nperiods>=2) 
	  // at least two full periods are inside our current l interval 
	  // change integration strategy to do a fine integration of one 2pi period 
	  // and then sum up all periods, re-scaled according to the power-law in l
	  {
	    multipleperiods=true;
	  
	    lb=lprev+2.*M_PI*nperiods/t2mt1; // make sure knext is an integer multiple of 2pi periods in l(t1+t2) away from lprev
	  }
	  
	  // P(l) is a power-law between lprev and la:
	  // P(l) = Plprev * (la/lprev)^alpha
	  // alpha = log(Plprev/Plnext) / log(la/lprev)
	  
	  double Plnext=pow10(Pkappa::logPkappalog(log10(lb)));
	  
	  double alpha=log10(Plnext/Plprev)/log10(lb/lprev);
	  
		
	  if(multipleperiods)
	    sum_b += plci_periods(Plprev,lprev,alpha,t2mt1,nperiods);
	  else
	    sum_b += plci(Plprev,lprev,alpha,t2mt1,lb);
	  
	  lprev  = lb;
	  Plprev = Plnext;
	  
	}
      
      }
      
      sum += (sum_a+sum_b)/(2.*M_PI*M_PI*sqrt(theta1*theta2));
      
      cout << thetaa << " " << thetab << " " << sum << endl;
      cov[ilta][iltb] = cov[iltb][ilta] = sum;
    }
    
            
  FITS *pFits = 0; 
  long naxis    =   2;      
  long naxes[2] = { rsteps, rsteps };

  try
  {
    const std::string fileName(argv[2]);            
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
     image[idx]=cov[i][j];
     idx++;
    }
  
  pFits->pHDU().write(1,rsteps*rsteps,image);
    
    return 0;
}

