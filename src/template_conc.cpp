#include "conc/template_conc.h"

#include "enfw/enfw.h"

#include <CCfits/CCfits>
using namespace CCfits;


int main(int argc, char **argv)
{
  
  if(argc!=3) {
   cout << "syntax: " << argv[0] << " [c200m] [output filename]" << endl;
   return 1;
  }
      
  FITS *pFits = 0; 
  long naxis    =   2;      
  long naxes[2] = { nv, nv };

  try
  {
    const std::string fileName(argv[2]);            
    pFits = new FITS(fileName, FLOAT_IMG , naxis , naxes );
  }
  catch (FITS::CantCreate)
  {
    cerr << "error, can't create fits file" << endl;
    return 1;
  } 
  
  // (1) parameters
  
  double c200m=atof(argv[1]);
  // concentration
  
  double Tmean[nv];	// mean of Sigma[v]
  double Tbuf[nv];	// buffer of Sigma[v]
  double Tprod[nv][nv]; // mean product of Sigma[v_i][v_j]
  
  for(int i=0; i<nv; i++)
  {
    Tmean[i]=0;
    for(int j=i; j<nv; j++)
      Tprod[i][j]=0;    
  }
  
  double s2siglog10c = sqrt(2.)*siglog10c;
  double dn=0;
  
  for(double dlog10c=-log10cdcmax; dlog10c<=log10cdcmax; dlog10c+=log10cstep) // iterate over concentrations
  {    
    double c = pow(10,dlog10c)*c200m;
    
    double w = erf((dlog10c+log10cstep/2.)/s2siglog10c)-erf((dlog10c-log10cstep/2.)/s2siglog10c);
    
    double rho0 = nfw_rho0_200m(c, 0.)/h/h; // h^2 Msol/R_S[Mpc]^3
    
    for(int i=0; i<nv; i++) 
    {
      double ror200m=pow(10,vmin+i*vstep);
      Tbuf[i] = rho0/c*nfw_rhop(ror200m*c)*w; // need to scale this later by r200m*(1+z)^3 and evaluate at radius in units of r200m
    }
    for(int i=0; i<nv; i++) 
    {
      Tmean[i] += Tbuf[i];    
      for(int j=i; j<nv; j++)
      {
       Tprod[i][j] += Tbuf[i]*Tbuf[j]/w;
      }
    }
    cout << 100.*((dlog10c+log10cdcmax)/(2.*log10cdcmax)) << "%                       \r" << flush; 
    dn += w;
  }
  
  for(int i=0; i<nv; i++)
  {
    Tmean[i] /= dn;
  }  
  
  for(int i=0; i<nv; i++)
  {
    for(int j=i; j<nv; j++)
      Tprod[i][j]=(Tprod[i][j]/dn-Tmean[i]*Tmean[j]); 
  }
  
  for(int i=0; i<nv; i++)
    for(int j=i-1; j>=0; j--) // i>j
      Tprod[i][j] = Tprod[j][i];
  

  valarray <float> image(0., nv*nv); // we're just going to jam our data in here ...
  int idx=0;
  for(int i=0; i<nv; i++)
    for(int j=0; j<nv; j++)
    {
     image[idx]=Tprod[i][j];
     idx++;
    }
  
  pFits->pHDU().write(1,nv*nv,image);
  return 0;
}

