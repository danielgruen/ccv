// compress an covariance matrix (off) to the system of annuli used

#include <assert.h>
#include "cosmology.h"
#include "off/template_off.h"
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
  double rmean[Nb]; // mean radius of annulus, area weighted

  for(int i=0; i<Nb; i++)
  {
    rmin2[i]=rmin[i]*rmin[i];
    rmax2[i]=rmax[i]*rmax[i];
    rmean[i]=2.*(pow(rmax[i],3)-pow(rmin[i],3)) / (3.*(rmax2[i]-rmin2[i]));
  }

  // (3) create conversion matrix from one set of annuli to the other
  
  double m200m=pow10(atof(argv[3])/100.); // hinv Msol
  double c200m=halomodel::cmz_200m_duffy(m200m,zlens); // assume Duffy here, need to scale covariance to angles correctly
  double rs_mpc=halomodel::r200mRadius(m200m/h, zlens)*h/c200m; // [Mpc]
  double rs=rs_mpc/DLen*arcminperrad; // arcmin
  
  
  
  double rhatmin[nv];
  double rhatmax[nv];

  rhatmin[0]=pow10(vmin)*rs;
  rhatmax[0]=rhatmin[0]*vfac;

  for(int j=1; j<nv; j++)
  {
    rhatmin[j]=rhatmax[j-1];
    rhatmax[j]=rhatmin[j]*vfac;
  }
  
  double rhatmin2[nv];
  double rhatmax2[nv];
  double rhatmean[nv];
  for(int j=0; j<nv; j++)
  {
    rhatmin2[j]=rhatmin[j]*rhatmin[j];
    rhatmax2[j]=rhatmax[j]*rhatmax[j];
    rhatmean[j]=2.*(pow(rhatmax[j],3)-pow(rhatmin[j],3)) / (3.*(pow(rhatmax[j],2)-pow(rhatmin[j],2)));
  }
  
  Matrix<double> Agamma(Nb,nv+1);
  Agamma.setZero();

  cout << "# input covariance matrix: " << nv << " bins from " << rhatmin[0] << " to " << rhatmax[nv-1] << " arcmin" << endl;
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
    Agamma(i,nv) = rhatmin2[0]/rhatmax2[j];
    //sum += Agamma(i,nv);
    
    Agamma(i,j) = -1.; // -kappa(theta)
    
    //cout << i << " " << sum << endl;
  }
    
  // (4) find matching version of sigma_R in units of r_s
  double sigma_r=sigma_r_of_m200m(m200m); // [Mpc] Rykoff+2016
  int nc = int(sigma_r/rs_mpc*100.+0.5);
  cout << "using sigma_r=" << sigma_r << " Mpc ~ " << nc/100. << " r_s" << endl;
    
  // (5) read in some covariance matrix
  
  Matrix<double> covhatg(nv+1,nv+1);
  
  {
    std::ostringstream ss;
    ss << nc;
    cout << "opening templates/off/cov_"+ss.str()+".fits" << endl;
    auto_ptr<FITS> pInfile(new FITS("templates/off/cov_"+ss.str()+".fits",Read,true));
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
	covhatg(i,j)=contents[idx];
	idx++;
      }
      
  }
  
  // (5) extrapolate covariance to central circle
  
  // most stupid thing: assume it's the same as innermost annulus
  for(int i=0; i<nv; i++)
  {
   covhatg(i,nv) = covhatg(nv,i) = covhatg(i,0); 
  }
  covhatg(nv,nv) = covhatg(0,0);
  
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
