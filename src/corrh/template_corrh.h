// template_corrh.h: some settings for calculation of corrh templates
// spanning angular range 0.1arcmin ... 60arcmin

const int rsteps=645;      // logarithmic steps in radius
const double rfactor=1.01; // 1% width annuli
const double lnrfactor=log(rfactor);
const double thetamin=0.1;  // arcmin
const double log10step=log10(rfactor);
const double thetamax=thetamin*pow(rfactor,rsteps);

const int mmin=80;         // 10*log(M/Msol) (no hinv)
const int mmax=159;

