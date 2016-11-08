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

