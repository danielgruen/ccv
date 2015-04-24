// compress an covariance matrix (ell) to the system of annuli used

#include <assert.h>
#include "cosmology.h"
#include "enfw/template_ell.h"
#include <TMV.h>
#include <iostream>

using namespace std;
using namespace tmv;

#include <CCfits/CCfits>

using namespace CCfits;


int main(int argc, char **argv)
{

  if(argc!=5)
  {
   cout << "syntax: " << argv[0] << " [zlens] [annulus definition file] [logM*100] [destination]" << endl;
   return 1;
  }

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
  
  double m200m=pow10(atof(argv[3])/100.);
  double c200m=halomodel::cmz_200m_duffy(m200m,zlens); // assume Duffy here, need to scale covariance to angles correctly
  double r200m=halomodel::r200mRadius(m200m/h, zlens)*h/DLen*arcminperrad; // arcmin
  
  int    nc = int(c200m*10.+0.5);
  assert(nc>=20);
  assert(nc<=100);
  
  double rhatmin[nv];
  double rhatmax[nv];
  rhatmin[0]=pow10(vmin)*r200m/c200m;
  rhatmax[0]=rhatmin[0]*vfac;
  for(int j=1; j<nv; j++)
  {
    rhatmin[j]=rhatmax[j-1];
    rhatmax[j]=rhatmin[j]*vfac;
  }
  
  Matrix<double> A(Nb,nv);

  cout << "# input covariance matrix: " << nv << " bins from " << rhatmin[0] << " to " << rhatmax[nv-1] << " arcmin" << endl;
  cout << "# output covariance matrix: " << Nb << " bins from " << rmin[0] << " to " << rmax[Nb-1] << " arcmin" << endl;
  
  // kappa = A kappahat
  for(int i=0; i<Nb; i++)
  {
    double areai = (rmax2[i]-rmin2[i]); // factor of PI consistently left out
    for(int j=0; j<nv; j++)
    {
      if(rmax[i]<=rhatmin[j] || rmin[i]>=rhatmax[j]) // no overlap
      {
      A(i,j) = 0.;
      continue;
      }
      double roverlapmin=max(rmin[i],rhatmin[j]);
      double roverlapmax=min(rmax[i],rhatmax[j]);
      //cout << i << " " << j << " overlapping between " << roverlapmin << " " << roverlapmax << endl;
      A(i,j) = (pow(roverlapmax,2)-pow(roverlapmin,2))/areai;  // factor of PI consistently left out
    }
  }
    
  // (4) read in some covariance matrix
  
  Matrix<double> covhat(nv,nv);
  
  {
    std::ostringstream ss;
    ss << nc;
    cout << "opening templates/ell/cov_"+ss.str()+".fits" << endl;
    auto_ptr<FITS> pInfile(new FITS("templates/ell/cov_"+ss.str()+".fits",Read,true));
    PHDU& image = pInfile->pHDU(); 
    valarray<double>  contents;      
    image.read(contents);
    long ax1(image.axis(0));
    long ax2(image.axis(1));
    
    assert(ax1==nv);
    assert(ax2==nv);
    
    int idx=0;
    for(int i=0; i<nv; i++)
      for(int j=0; j<nv; j++)
      {
	covhat(i,j)=contents[idx];
	idx++;
      }
      
  }
  
  Matrix<double> cov(Nb,Nb);
  
  cov = A*covhat*A.transpose();

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
