// template_corrh.h: some settings for calculation of corrh templates
// spanning angular range 0.04arcmin ... 180arcmin

const int rsteps=846;      // logarithmic steps in radius
const double rfactor=1.01; // 1% width annuli
const double lnrfactor=log(rfactor);
const double thetamin=0.04;  // arcmin
const double log10step=log10(rfactor);
const double thetamax=thetamin*pow(rfactor,rsteps);

const int mmin=80;         // 10*log(M/Msol) (no hinv)
const int mmax=159;

