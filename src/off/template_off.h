#include <cmath>

const double pcen=0.78; // probability of correct centering

// v in units of r_s
const int nv=1500;
const double vmin=-3.5;
const double vfac=1.01;
const double vstep=log10(vfac);

double rayleigh_cumulative(double x)
// P(X<x) where X is a Rayleigh-distributed random variable (i.e. X=sqrt(Gauss1^2+Gauss2^2))
{
  if(x<=0) return 0.;
  return 1.-exp(-0.5*x*x);
}

double lambda_of_m200m(double m200m) // expects hinv Msol input
{
  return 30.*m200m/0.7/2.35e14; // Melchior+2016, but without Evrard+2014 correction
}

double sigma_r_of_m200m(double m200m) // expects hinv Msol input
{
  double lambda = lambda_of_m200m(m200m);
  double r_lambda=0.1*pow(lambda/10.,0.2)/h; // [Mpc]
  double sigma_r=exp(-1.13)*r_lambda; // [Mpc] Rykoff+2016
  return sigma_r;
}
