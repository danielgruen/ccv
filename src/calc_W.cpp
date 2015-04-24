// test3: calculate W(theta), integrating xi out to 200Mpc

#include <CCfits/CCfits>
using namespace CCfits;

#include "corrh/corrh.h"
#include <iostream>
#include <cstdlib>

using namespace std;
using namespace xilinear;

int main(int argc, char **argv)
{ 
  if(argc!=2) {
    cerr << "syntax: " << argv[0] << " [redshift]" << endl;
    return 1;
  }
 

  double z=atof(argv[1]);
  double Da1pz=angularDiameterDistance(0, z, 10000)*h*(1.+z)/arcminperrad; // comoving Mpc hinv per arcmin

  cout << OmegaM << " " << OmegaL << " " << OmegaB << " " << sigma8 << " " << ns << " " << argv[1] << endl;
    
  for(double theta=0.01; theta<60.; theta*=1.01)
  {
    double Da1pztheta2=pow(Da1pz*theta,2); 
    
    double zetamin=0.;
    double zetamax=10.; // integrate DeltaZeta out to 100 hinv Mpc comoving
    int zetasteps=5000;
    double dzeta=(zetamax-zetamin)/double(zetasteps);
    
    double s=0.;
    double integrandprev=twopclin(0.5*log10(Da1pztheta2), z);
    
    for(int izeta=0; izeta<zetasteps; izeta++) 
    {
      double zetanext=zetamin+double(izeta+1)*dzeta;
      double integrandnext=twopclin(0.5*log10(Da1pztheta2+zetanext*zetanext),z);
      s += dzeta*(integrandprev+integrandnext)/2.;
      
      integrandprev=integrandnext;
    }
    
    zetamin=zetamax;
    zetamax=200.; // integrate DeltaZeta out to 200 hinv Mpc comoving
    zetasteps=10000;
    dzeta=(zetamax-zetamin)/double(zetasteps);
    
    for(int izeta=0; izeta<zetasteps; izeta++) 
    {
      double zetanext=zetamin+double(izeta+1)*dzeta;
      double integrandnext=twopclin(0.5*log10(Da1pztheta2+zetanext*zetanext),z);
      s += dzeta*(integrandprev+integrandnext)/2.;
      
      integrandprev=integrandnext;
    }
    
    s *= 2.*Da1pz*Da1pz;
    
    cout << log10(theta) << " " << log10(s) << endl;  
  }
  
  return 0;
}

