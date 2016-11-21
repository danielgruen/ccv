#ifndef _PROFILE_H_
#define _PROFILE_H_

#include "../cosmology.h"
//#include <complex.h>
//#undef complex
#include <complex>

double b_nfw(double c200m)
{
  return log(1.+c200m)-c200m/(1.+c200m);
}

double b_bmo(double c200m, double tau)
{
  double tau2   = tau*tau;
  double c200m2 = c200m*c200m;
  double tau2pc200m2 = tau2+c200m2;
  double c200mp1 = c200m+1.;
  return tau2 / (2.*pow(tau2+1.,3)*(c200mp1)*(tau2pc200m2)) *   // Oguri & Hamana 2011, Eqn. 10
             (
                (tau2+1.)*c200m*(c200m*c200mp1-tau2*(c200m-1.)*(2.+3.*c200m)-2.*tau2*tau2)
              + tau*c200mp1*tau2pc200m2*
                   (
                     2.*(3.*tau2-1.)*atan(c200m/tau)
                   + tau*(tau2-3.)*log(tau2*pow(c200mp1,2)/tau2pc200m2)
                   )
             );
}

double bmo_rho0_200m(double c200m, double tau, double z)
// scale density of bmo profile
{
  double rhom  = OmegaM*3.*H0*H0/(8.*M_PI*Gs)*pow(1.+z,3);   // Msol/Mpc^3; NO h^2
  return 200.*rhom*pow(c200m,3)/ (3.*b_bmo(c200m,tau));
  
}

// taken from http://kipac.stanford.edu/collab/research/lensing/ample/ample.c
double ample_kappa_NFW_trunc2 ( double x, double tau )
{
  // x: radius in units of r_s
  // tau: truncation radius in units of r_s 
  
  double u = x*x; // quadratic form
  
  // the profile assumed is rho0/(4pi)/(x(1+x)^2)*C^4/(C^2+x^2)^2
  
  if ( fabs ( u - 1.0 ) < 0.00008 ) {
    double kappa1 = ample_kappa_NFW_trunc2 ( sqrt(0.9999), tau );
    double kappa2 = ample_kappa_NFW_trunc2 ( sqrt(1.0001), tau );
    double frac = ( u - 0.9999 ) / 0.0002;
    return kappa2*frac + kappa1*(1.-frac);
  }
  
  double lnC = log ( tau );
  //double ln2C = lnC + M_LN2;
  double um1 = u - 1.0;
  double rootu = sqrt ( u );
  //double lnu = log ( u );
  double u2 = u * u;
  //double u3 = u2 * u;

  complex<double> arccosu = std::acos ( 1.0 / rootu );
  complex<double> sqrtum1, arccossqrt;
  
  sqrtum1 = std::sqrt ( um1 );
  if ( u < 1.0 ) sqrtum1 = -1.*sqrtum1;
  arccossqrt = arccosu / sqrtum1;
  
  double C2 = tau * tau;
  double C4 = C2 * C2;
  double C6 = C4 * C2;
  //double C8 = C4 * C4;
  double rootC2u = sqrt ( C2 + u );
  double messylog = log ( ( rootC2u - tau ) / rootu );

  complex<double> psi1 = C4 / ( 2.0 * ( C2 + 1.0 ) * ( C2 + 1.0 ) * ( C2 + 1.0 ) * u ) *
    ( 2.0 * ( C2 + 4.0 * u - 3.0 ) * arccossqrt 
      + ( M_PI * ( 3.0 * C2 - 1.0 ) + 2.0 * tau * ( C2 - 3.0 ) * lnC ) / tau
      + ( - C2 * tau * M_PI * ( 4.0 * u + 3.0 * C2 - 1.0 )
	  + ( 2.0 * C6 - 6.0 * C4
	      + u * ( 3.0 * C4 - 6.0 * C2 - 1.0 ) ) * messylog ) /
      ( C2 * tau * rootC2u ) );

  complex<double> psi2 = C4 / ( 4.0 * ( C2 + 1.0 ) * ( C2 + 1.0 ) * ( C2 + 1.0 ) *
		  u2  ) *
    ( ( 10.0 - 8.0 * u - 6.0 * C2 ) * arccossqrt - 4.0 * ( C2 - 3.0 ) * lnC
      + ( 2.0 * C6 + 3.0 * C4 * ( u - 2 ) - u - 6.0 * C2 * u ) /
      ( C2 * ( C2 + u ) )
      + ( 2.0 * u * ( C2 + u ) * ( 3.0 * C4 - 6.0 * C2 - 1.0 )
	  - ( 3.0 * u + 2.0 * C2 ) *
	  ( 2.0 * C6 + 3.0 * C4 * ( u - 2 ) - u
	    - 6.0 * C2 * u ) ) / ( C2 * tau * ( C2 + u ) * rootC2u ) * messylog
      - 2.0 * M_PI / tau * ( 3.0 * C2 - 1.0 )
      + M_PI / ( ( C2 + u ) * rootC2u ) *
      ( 2.0 * ( C2 + u ) * ( 3.0 * C2 - 1.0 )
	+ u * ( 4.0 * u + 3.0 * C2 - 1.0 ) )
      + 2.0 * ( C2 + 1.0 ) / um1 * ( 1.0 - arccossqrt ) + 8.0 );

  
  return 4.*real(psi1+u*psi2);
}


#endif

