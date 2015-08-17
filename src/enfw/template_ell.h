const double c200mstep=0.01;

// q=c/a
const double qmin=0.1;
const double qmax=1.0;
const double qstep=0.01;
//const int    nq = (qmax-qmin)/qstep+1;
const int    nq = 91;
const int    iqmin=qmin/qstep;
const int    iqmax=qmax/qstep;

// cos(alpha)
const double camin=0.;
const double camax=1.;
const double castep=0.01;
//const int    nca = (camax-camin)/castep+1;
const int    nca = 101;

const int nv=1500;
const double vmin=-3.;          // can't just change these, would have to run profile_kappa_gamma again
const double vfac=1.01;
const double vstep=log10(vfac);


