#include <cmath>

const double log10cstep=0.005;
const double siglog10c=0.18;
const double log10cdcmax=4.*siglog10c; // 4sigma deviations max 

// v in units of r_200m
const int nv=1500;
const double vmin=-3.5;
const double vfac=1.01;
const double vstep=log10(vfac);
