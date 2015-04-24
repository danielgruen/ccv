// combine covariance matrices due to haloes of different mass into one

#include <CCfits/CCfits>
using namespace CCfits;

#include "corrh/corrh.h"
#include "corrh/template_corrh.h"
#include "profile/profile.h"
#include "enfw/enfw.h"
#include <iostream>
#include <cstdlib>
#include <string.h>


using namespace std;


int main(int argc, char **argv)
{
  if(argc!=4) {
   cerr << "syntax: " << argv[0] << " [z] [output .fits] [filename prefix all before 100logm, e.g. templates/corrh/cov_0.2453300_]" << endl;
   return 1;
  }
  double z=atof(argv[1]);

  valarray <double> cov(0., rsteps*rsteps);
  
  valarray <double> dcov_prev(0., rsteps*rsteps);
  
  const double logMstep=0.1;
  const double Mfac=pow10(logMstep);
  
  for(int i=mmin; i<=mmax; i++)
  {
    ostringstream ss;
    ss << argv[3] << i << ".fits";
    cerr << "opening " << ss.str() << endl;
    auto_ptr<FITS> pInfile(new FITS(ss.str(),Read/*,0,true*/)); // opens&reads 
    valarray <double> dcov(0., rsteps*rsteps);
    pInfile->pHDU().read(dcov);
    
    const double M = pow10(i/10.);
    const double bias = halomodel::bias_tinker(M, z);
    const double dndM = halomodel::dndM(M, z);
   
    // dcov: h^2 Msol^2 per proper Mpc^4 if there were 1 such halo per com. hinv^3 Mpc^3 Lagrangian volume
    
    for(int j=0; j<rsteps*rsteps; j++) // power-law integration in each (theta_1,theta_2)
    {
	dcov[j] *= dndM*bias;
    }
    
    // dcov now: h^2 Msol^2 per proper Mpc^4 per 1 solar mass
    
    if(i>mmin) {
    
      for(int j=0; j<rsteps*rsteps; j++) // power-law integration in each (theta_1,theta_2)
      {
	
	if(dcov[j]==0 || dcov_prev[j]==0)
	  continue;  
	
	double alpha = log10(dcov[j]/dcov_prev[j])/logMstep; // mass step is 0.1 dex always
	if(j==0)
	  cout << M << " " << bias << " " << dndM << " " << dcov[j] << " " << dcov_prev[j] << " " << alpha << " " << Mfac << " " << dcov_prev[j]*M/Mfac/(alpha+1.)*(pow(Mfac,alpha+1.)-1.) << endl;
	  
	cov[j] += dcov_prev[j]*M/Mfac/(alpha+1.)*(pow(Mfac,alpha+1.)-1.); // power-law integration
      }
    
    }
        
    dcov_prev = dcov;
  }
  
  FITS *pFits = 0; 
  long naxis    =   2;      
  long naxes[2] = { rsteps, rsteps };

  try
  {
    const std::string fileName(argv[2]);            
    pFits = new FITS(fileName, DOUBLE_IMG , naxis , naxes );
  }
  catch (FITS::CantCreate)
  {
    cerr << "error, can't create fits file" << endl;
    return 1;
  } 

  pFits->pHDU().write(1,rsteps*rsteps,cov);
  
  return 0;
}

