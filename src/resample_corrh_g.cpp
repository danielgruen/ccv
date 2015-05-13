// compress an angular covariance matrix (corrh, lss) to the system of annuli used in the analysis, tranform to gamma

#include <assert.h>
#include "cosmology.h"
#include "corrh/template_corrh.h"

#include <TMV.h>
#include <iostream>

using namespace std;
using namespace tmv;

#include <CCfits/CCfits>

using namespace CCfits;

#include <fstream>

int main(int argc, char **argv)
{

  if(argc!=5)
  {
   cout << "syntax: " << argv[0] << " [redshift] [annulus definition file] [source] [destination]" << endl;
   return 0;
  }

  // (1) generate annulus structure
  double zlens=atof(argv[1]);
  assert (zlens>0.);
  double DLen = angularDiameterDistance(0,zlens,10000)*h; // hinv Mpc per rad

  ifstream ann(argv[2]);
  int Nb; ann >> Nb;

  double rmin[Nb];
  double rmax[Nb];

  for(int i=0; i<Nb; i++)
  {
    ann >> rmin[i] >> rmax[i];
    assert(ann.eof()==0);
  }
  string buf;
  ann >> buf;
  assert(ann.eof());

  double rmin2[Nb];
  double rmax2[Nb];
  double rmean[Nb]; // mean radius of annulus, area weighted

  for(int i=0; i<Nb; i++)
  {
    rmin2[i]=rmin[i]*rmin[i];
    rmax2[i]=rmax[i]*rmax[i];
    rmean[i]=2.*(pow(rmax[i],3)-pow(rmin[i],3)) / (3.*(rmax2[i]-rmin2[i]));
  }

  // (3) create conversion matrix from one set of annuli to the other
  
  // Nb is set to the final number of bins now, rmin and rmax contain their minimum and maximum radii

  double rhatmin[rsteps];
  double rhatmax[rsteps];
  double rhatmin2[rsteps];
  double rhatmax2[rsteps];
  double rhatmean[rsteps];
  rhatmin[0]=thetamin;
  rhatmax[0]=rhatmin[0]*rfactor;
  for(int j=1; j<rsteps; j++) 
  {
    rhatmin[j]=rhatmax[j-1];
    rhatmax[j]=rhatmin[j]*rfactor;
    rhatmin2[j]=rhatmin[j]*rhatmin[j];
    rhatmax2[j]=rhatmax[j]*rhatmax[j];
    rhatmean[j]=2.*(pow(rhatmax[j],3)-pow(rhatmin[j],3)) / (3.*(pow(rhatmax[j],2)-pow(rhatmin[j],2)));
  }
  
  Matrix<double> Agamma(Nb,rsteps+1);

  cout << "# input covariance matrix: " << rsteps << " bins from " << rhatmin[0] << " to " << rhatmax[rsteps-1] << " arcmin" << endl;
  cout << "# output covariance matrix: " << Nb << " bins from " << rmin[0] << " to " << rmax[Nb-1] << " arcmin" << endl;
  
  // gamma = Agamma kappahat'
  for(int i=0; i<Nb; i++)
  {
    //double sum=0.;
    int j=0;
    while(rhatmean[j]<rmean[i]) {
     Agamma(i,j) = rhatmax2[j]-rhatmin2[j];  // factor of PI consistently left out
     j++;
    }
    assert(j>0);
    for(int J=0; J<j; J++) // normalize
    {
     Agamma(i,J) /= rhatmax2[j-1]; 
     //sum += Agamma(i,J);
    }
    Agamma(i,rsteps) = rhatmin2[0]/rhatmax2[j];
    //sum += Agamma(i,rsteps);
    
    Agamma(i,j) = -1.; // -kappa(theta)
    
    //cout << i << " " << sum << endl;
  }
  
  // (4) read in some covariance matrix
  
  Matrix<double> covhatg(rsteps+1,rsteps+1);
  
  {
    auto_ptr<FITS> pInfile(new FITS(argv[3],Read,true));
    PHDU& image = pInfile->pHDU(); 
    valarray<double>  contents;      
    image.read(contents);
    long ax1(image.axis(0));
    long ax2(image.axis(1));
    
    assert(ax1==rsteps);
    assert(ax2==rsteps);
    
    int idx=0;
    for(int i=0; i<rsteps; i++)
      for(int j=0; j<rsteps; j++)
      {
	covhatg(i,j)=contents[idx];
	idx++;
      }
      
  }
  
  // (5) extrapolate covariance to central circle
  
  // most stupid thing: assume it's the same as innermost annulus
  for(int i=0; i<rsteps; i++)
  {
   covhatg(i,rsteps) = covhatg(rsteps,i) = covhatg(i,0); 
  }
  covhatg(rsteps,rsteps) = covhatg(0,0);
  
  
  Matrix<double> cov(Nb,Nb);
  
  cov = Agamma*covhatg*Agamma.transpose();
  

  {
    FITS *pFits = 0; 
    long naxis    =   2;      
    long naxes[2] = { Nb, Nb };
    try
    {
      const std::string fileName(argv[4]);            
      pFits = new FITS(fileName, DOUBLE_IMG , naxis , naxes );
    }
    catch (FITS::CantCreate)
    {
      cerr << "error, can't create fits file" << endl;
      return 1;
    } 

    valarray <float> image(0., Nb*Nb); // we're just going to jam our data in here ...
    int idx=0;
    for(int i=0; i<Nb; i++)
      for(int j=0; j<Nb; j++)
      {
      image[idx]=cov(i,j);
      idx++;
      }
    
    pFits->pHDU().write(1,Nb*Nb,image);
  }
  
  return 0;
}
