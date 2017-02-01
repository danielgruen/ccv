#include "enfw/enfw.h"

#include <CCfits/CCfits>
using namespace CCfits;

#include "off/template_off.h"
  
double Tmean[nv];	// mean of Sigma[v]
double Tbuf[nv];	// buffer of Sigma[v]
double Tprod[nv][nv]; // mean product of Sigma[v_i][v_j]

int main(int argc, char **argv)
{
  
  if(argc!=3) {
   std::cout << "syntax: " << argv[0] << " [sigma_off/r_s] [output filename]" << std::endl;
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
    std::cerr << "error, can't create fits file" << std::endl;
    return 1;
  } 
  
  // (1) parameters
  
  double sigma_off=atof(argv[1]); // [r_s]
  // off-centering sigma in each of x,y in case it's off-centered
  
  
  for(int i=0; i<nv; i++)
  {
    Tmean[i]=0;
    for(int j=i; j<nv; j++)
      Tprod[i][j]=0;    
  }
  
  double dn=0;
  
  // (a) centered component, w=pcen
  double w=pcen;
  for(int i=0; i<nv; i++) 
  {
    double rors=pow(10,vmin+i*vstep);
    Tbuf[i] = nfw_rhop_offcentered(rors,0.)*w; // need to scale this later by rho0/c*r200m*(1+z)^3 and evaluate at radius in units of r_s
  }
  for(int i=0; i<nv; i++) 
  {
    Tmean[i] += Tbuf[i];    
    for(int j=i; j<nv; j++)
    {
      Tprod[i][j] += Tbuf[i]*Tbuf[j]/w;
    }
  }
  dn += w;
  
  const int psteps=100.;
  const double smax=6.;
  for(int i=0; i<psteps; i++)
  {    
    double obefore=i*smax/psteps;
    double oafter=obefore+smax/psteps;
    w=(1.-pcen)*(rayleigh_cumulative(oafter)-rayleigh_cumulative(obefore));
    double o=(obefore+oafter)/2.*sigma_off;
    
    
    for(int i=0; i<nv; i++) 
    {
      double rors=pow(10,vmin+i*vstep);
      Tbuf[i] = nfw_rhop_offcentered(rors,o)*w; // need to scale this later by rho0/c*r200m*(1+z)^3 and evaluate at radius in units of r_s
    }
    for(int i=0; i<nv; i++) 
    {
      Tmean[i] += Tbuf[i];    
      for(int j=i; j<nv; j++)
      {
       Tprod[i][j] += Tbuf[i]*Tbuf[j]/w;
      }
    }
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

