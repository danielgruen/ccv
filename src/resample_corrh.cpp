// compress an angular covariance matrix (corrh, lss) to the system of annuli used in the analysis

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

  // (1) generate annulus structure

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

  for(int i=0; i<Nb; i++)
  {
    rmin2[i]=rmin[i]*rmin[i];
    rmax2[i]=rmax[i]*rmax[i];
  }

  // (3) create conversion matrix from one set of annuli to the other
  
  // Nb is set to the final number of bins now, rmin and rmax contain their minimum and maximum radii

  double rhatmin[rsteps];
  double rhatmax[rsteps];
  rhatmin[0]=thetamin;
  rhatmax[0]=rhatmin[0]*rfactor;
  for(int j=1; j<rsteps; j++) 
  {
    rhatmin[j]=rhatmax[j-1];
    rhatmax[j]=rhatmin[j]*rfactor;
  }
  
  Matrix<double> A(Nb,rsteps);

  cout << "# input covariance matrix: " << rsteps << " bins from " << rhatmin[0] << " to " << rhatmax[rsteps-1] << " arcmin" << endl;
  cout << "# output covariance matrix: " << Nb << " bins from " << rmin[0] << " to " << rmax[Nb-1] << " arcmin" << endl;
  
  // kappa = A kappahat
  for(int i=0; i<Nb; i++)
  {
    double areai = (rmax2[i]-rmin2[i]); // factor of PI consistently left out
    for(int j=0; j<rsteps; j++)
    {
      if(rmax[i]<=rhatmin[j] || rmin[i]>=rhatmax[j]) // no overlap
      {
      A(i,j) = 0.;
      continue;
      }
      double roverlapmin=max(rmin[i],rhatmin[j]);
      double roverlapmax=min(rmax[i],rhatmax[j]);
      A(i,j) = (pow(roverlapmax,2)-pow(roverlapmin,2))/areai;  // factor of PI consistently left out
    }
  }
  
  //cout << A << endl;   
  
  
  // (4) read in some covariance matrix
  
  Matrix<double> covhat(rsteps,rsteps);
  
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
	covhat(i,j)=contents[idx];
	idx++;
      }
      
  }
  
  Matrix<double> cov(Nb,Nb);
  
  cov = A*covhat*A.transpose();
  
  //cout << cov << endl;
  

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
