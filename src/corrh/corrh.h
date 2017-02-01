#include "../cosmology.h"
#include "../enfw/enfw.h"
#include "../profile/profile.h"

#include <cstdlib>
#include <string.h>

namespace halosigma { // stuff to calculate surface mass density due to off-center haloes 

const double tmindefault=1.e-3;  // default minimum radius around halo; all matter inside is dumped onto an annulus assuming NFW
const double ringresolution=5.;  // how many steps in halo ring radius per sky annulus?
const double rsresolution=5.;    // how many steps in halo ring radius per halo scale radius?
//const double hardtruncation=100.; // hard truncation of surface mass density at this many rs distance
//const double tau=20.;
//const double hardtruncation=20.; // hard truncation of surface mass density at this many rs distance
//const double tau=10.;

void halosigma(
  
  // calculates matter density due to one halo of fixed mass and off-center position in a set of sky annuli
  // assumes BMO halo with squared truncation at 'tau' and hard cut-off at 'hardtruncation'
  
  const double r,	// separation of halo from centre [Mpc]
  const double A, 	// amplitude to be multiplied with BMO profile [Msol/Mpc^2]
  const double rs, 	// scale radius [Mpc]
  const int rsteps, 	// number of sky annuli
  const double rfactor, // radius factor beween subsequent sky annuli
  double *dens, 	// density in sky annuli; array will be overwritten with result of the calculation [Msol/Mpc^2]
  const double *rmin, 	// inner radius of each sky annulus [Mpc]
  const double *rmax, 	// outer radius of each sky annulus
  const double *rmin2, 	// square of inner radius of each sky annulus [Mpc^2]
  const double *rmax2,	// square of outer radius of each sky annulus
  const double *area,	// area of each sky annulus [Mpc^2]
  const double *zeros,	// array of zeroes, copied into dens before filling it with actual density
  int &nonzero_minannulus,
  int &nonzero_maxannulus,
  const double tau=20.,
  const double hardtruncation=100.  
  )
{
  
  // (0) clear density
  memcpy(dens, zeros, rsteps*sizeof(double));
  
  const double tmax=std::min(rmax[rsteps-1]+r,hardtruncation*rs);
  
  nonzero_minannulus=rsteps-1;
  nonzero_maxannulus=0;
  
  if(r-tmax>rmax[rsteps-1]) return;
  
  // (1) find annulus where this halo belongs
  int r0=0;
  while(r0<rsteps-1 && rmax[r0]<r)
    r0++;
  int r1=r0;
  
  // (2) dump central matter to annulus / annuli
  
  // add density due to BMO halo with amplitude A [Msol per proper Mpc^2], scale radius rs [proper Mpc] at position r [proper Mpc]
  // to array dens starting at rmin [Mpc] in rsteps logarithmic steps of width log10step  
  double tmin;
  if(r0>0 && r0<rsteps-1) { // found the halo core somewhere inside the sky annuli region
    if(fabs(rmin[r0]-r)<1.e-4) { // sitting straight on inner sky annulus border
     tmin=tmindefault; // 1 kpc 
     dens[r0-1]=dens[r0]=nfw_Rhop(tmin/rs)*M_PI*tmin*tmin/2.; // un-truncated NFW because at very small radii that's pretty reasonable
     nonzero_maxannulus=r0; nonzero_minannulus=r0-1;
     assert(!isnan(dens[r0]));
    }
    else if(r-rmin[r0]<rmax[r0]-r) { // closer to inner radius of annulus
     tmin=std::min(tmindefault,r-rmin[r0]);
     dens[r0]=nfw_Rhop(tmin/rs)*M_PI*tmin*tmin;
     nonzero_maxannulus=nonzero_minannulus=r0;
     assert(!isnan(dens[r0]));
    }
    else if(fabs(rmax[r0]-r)<1.e-4) // straight on outer sky annulus border
    {
     tmin=tmindefault; // 1 kpc 
     dens[r0+1]=dens[r0]=nfw_Rhop(tmin/rs)*M_PI*tmin*tmin/2.; // un-truncated NFW because at very small radii that's pretty reasonable
     nonzero_maxannulus=r0+1; nonzero_minannulus=r0;
     assert(!isnan(dens[r0]));
    }
    else { // closer to outer radius of annulus
     tmin=std::min(tmindefault,rmax[r0]-r);
     dens[r0]=nfw_Rhop(tmin/rs)*M_PI*tmin*tmin;
     nonzero_maxannulus=nonzero_minannulus=r0;
     assert(!isnan(dens[r0]));
    }
  }
  else if (r0==0) { // inside innermost sky annulus
   tmin=std::max(rmin[0]-r,1.e-4); 
  }
  else { // r0==rsteps-1, in or outside outermost sky annulus
   if(fabs(r-rmax[r0])<1.e-4) // sitting near outer border of outermost sky annulus
   { 
     tmin=tmindefault;
     dens[r0]=nfw_Rhop(tmin/rs)*M_PI*tmin*tmin/2.; 
     nonzero_maxannulus=nonzero_minannulus=r0;
   }
   else if(r-rmin[r0]<rmax[r0]-r) { // closer to inner radius of annulus
     tmin=std::min(tmindefault,r-rmin[r0]);
     dens[r0]=nfw_Rhop(tmin/rs)*M_PI*tmin*tmin;
     nonzero_maxannulus=nonzero_minannulus=r0;
     assert(!isnan(dens[r0]));
   } else if (rmax[r0]>r) { // closer to outer radius of annulus but inside
     tmin=std::min(tmindefault,rmax[r0]-r);
     dens[r0]=nfw_Rhop(tmin/rs)*M_PI*tmin*tmin; 
     nonzero_maxannulus=nonzero_minannulus=r0;
     assert(!isnan(dens[r0]));
   } else { // outside
     tmin=std::max(tmindefault,r-rmax[r0]);
     nonzero_maxannulus=nonzero_minannulus=r0;
   }
  }
    
  // (3) iterate halo rings and dump their mass on sky annuli
  
  double tprev=tmin;
  double tnext;
  
  double sum=0.;
  
  for(double t=tmin; t<=tmax;) // current radius around halo centre, Mpc
  {
    double rmt=fabs(r-t);
    double rpt=r+t;
    double rt2=2.*r*t;
    double r2t2=r*r+t*t;
    
    //cerr << "halo ring of radius " << t << " reaching from sky radius " << rmt << " - " << rpt << endl;

    // find largest affected radius
    while(r1<rsteps-1 && rmin[r1]<rpt)
      r1++;

    if(r<t) r0=r1;
    
    // find smallest affected annulus
    while(r0>0 && rmax[r0]>rmt)
      r0--;
    
    if(rmt>rmax[rsteps-1]) 
    {
      //cerr << "we are far enough outside the fov; breaking." << endl; 
      break; // far enough outside, break
    }

    double tnext = std::min(t + rmin[r0]*(rfactor-1.)/ringresolution,t+rs/rsresolution); // 1/5 of smallest affected sky annulus
    
    double sigma  = ample_kappa_NFW_trunc2(t/rs,tau);
    double sigang = sigma*0.5*(tnext*tnext-tprev*tprev); 
    // 2Msol / rad; down there we always multiply with half the covered angle; factor 0.5 because we cover all area twice
    
    double thetaprev=0.;
    
    if(r0<nonzero_minannulus)
      nonzero_minannulus=r0;
    if(r1>nonzero_maxannulus)
      nonzero_maxannulus=r1;
        
    for(int ir=r0; ir<=r1; ir++) // iterate all relevant sky annuli
    {
      if(rmax[ir]<rmt) {
	continue;
      }
      if(rmin[ir]>rpt) {
	break;
      }
      
      if(rmt>rmin[ir] && rpt<rmax[ir]) { // halo ring is all inside sky annulus
	dens[ir] += M_PI*sigang;
	assert(!isnan(dens[ir]));
	break;
      }
      
      if(rmt<rmin[ir]) {   // part but not all of halo ring is inside inner border of sky annulus
	if(rpt<rmax[ir]) { // halo ring does not reach outside sky annulus -- halo ring cuts annulus twice on inside border
	  double theta = thetaprev;
	  if(!theta) theta=acos((rmin2[ir]-r2t2)/rt2);
	  if(isnan(theta)) theta=0.;
	  thetaprev=0.;
	  
	  dens[ir] += theta*sigang;
	  assert(!isnan(dens[ir]));
	} else {           // halo ring does reach outside sky annulus -- halo ring cuts annulus four times, twice on each border 
	  double theta1 = thetaprev;
	  if(!theta1) theta1=acos((rmin2[ir]-r2t2)/rt2); // if no previous exists, rare case
	  double theta2;
	  if ((rmax2[ir]-r2t2)/rt2>=1) theta2 = thetaprev = 0.;
	  else theta2 = thetaprev = acos((rmax2[ir]-r2t2)/rt2);
	  dens[ir] += sigang*(theta1-theta2);
	  assert(!isnan(dens[ir]));
	}
      } else if(rmt<rmax[ir]) {  // part of halo ring is inside outer border of sky annulus (but it doesn't cut inside border)
	  double theta=acos((rmax2[ir]-r2t2)/rt2);
	  if(isnan(theta)) theta=M_PI;
	  thetaprev=theta;
	  dens[ir] += sigang*(M_PI-theta);
	  assert(!isnan(dens[ir]));
      }
    }
    
    tprev = t;
    t = tnext;
    
    
    sum += sigang*M_PI;
  }
  
  // (4) normalize mass by area and A --> surface density, Msol/Mpc^2
  
  for(int i=nonzero_minannulus; i<=nonzero_maxannulus; i++) dens[i] *= (A/area[i]);
}

}
