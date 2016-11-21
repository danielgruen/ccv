#include "enfw/enfw.h"
#include "enfw/template_ell.h"
#include "filter/filter.h"
#include <CCfits/CCfits>

using namespace CCfits;

double Tq[nq];
double Tamp[nq];
double Tmean[nv];	// mean of Sigma[v]
double Tprod[nv][nv]; // mean product of Sigma[v_i][v_j]
double Tk[nv][nq][nca]; // big array to hold the profiles; cos(alpha) striding is the cheapest

double q,c200m; // present value for transformation functions

double weight_c200m_q(double tc200m, double tq)
// bilinear interpolation in e and c200m
{
  if(fabs(tc200m-c200m)>c200mstep) return 0;
  if(fabs(tq-q)>qstep) return 0;
  
  return (c200mstep-fabs(tc200m-c200m))/c200mstep*(qstep-fabs(tq-q))/qstep;
}

int main(int argc, char **argv)
{
  
  if(argc!=3) {
   cout << "syntax: " << argv[0] << " [c200m] [output filename]" << endl;
   return 1;
  }

#ifndef NO_OMP
  omp_set_num_threads(1);
#endif 

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
  
  c200m=atof(argv[1]);
  // concentration
  
  // set up a watchdog: which bin is 1' at z=0.245 and M=1.e14?
  
  double rs_14_1 = halomodel::r200mRadius(1.e14/h, 0.24533)*h  /  (angularDiameterDistance(0,0.24533,10000)*h)  * arcminperrad  /  c200m; // arcmin
  double v_14_1  = 1. / rs_14_1; // v at 1 arcmin
  int    watchdog = (log10(v_14_1)-vmin)/vstep;
  
  cout << "# rs is " << rs_14_1 << " arcmin" << endl;
  

  // (2) library tables
  
  // (2a) density amplitudes as function of c200m and ellipticity
  
  const string cn[] = {"c200m", "q", "rho_0"};
  const Type   ct[] = {doubleType, doubleType, doubleType};
  const int    ck[] = {1, 2, 4};
  
  ifstream i("lut/rho0.tab"); 
  ObjectCollection *amp = new ObjectCollection(i, ck, cn, ct, 3);
    
  // (2b) amplitude for given concentration as a function of ellipticity
 
  
  for(int i=0; i<nq; i++)
  {
   q=qmin+qstep*double(i);
   Tq[i]=q;
   amp->transformColumnNew("w",weight_c200m_q,"c200m","q");
   Tamp[i]=amp->averageWeighted("rho_0","w");
   
   // normalization
   q=1;
   amp->transformColumnNew("w",weight_c200m_q,"c200m","q");
   Tamp[i] /= amp->averageWeighted("rho_0","w");
   
   cout << Tq[i] << " " << Tamp[i] << " " << amp->sum("w") << endl;
  }
  
  delete amp; amp=0;
  
  // (2c) read profiles
  
  double Tv[nv];
  Tv[0]=vmin;
  for(int i=1; i<nv; i++)
  {
   Tv[i]=Tv[i-1]+vstep;
  }
  
  //cout << "watchdog v: " << Tv[watchdog] << " should be " << v_14_1 << endl;
  
  
  {
    auto_ptr<FITS> pInfile(new FITS("lut/profiles.fits",Read,true));
    PHDU& profiles = pInfile->pHDU(); 
    long ax1(profiles.axis(0));
    long ax2(profiles.axis(1));
    long ax3(profiles.axis(2));
      
    assert(ax1==nv);
    assert(ax2==nca);
    assert(ax3==nq);
    
    valarray<double>  contents;      
    profiles.read(contents);
      
    int idx=0;
    for(int iq=iqmin; iq<=iqmax; iq++)
    {
      for(int ica=0; ica<nca; ica++)
      {
	for(int i=0; i<nv; i++)
	{
	  Tk[i][iq-iqmin][ica]=contents[idx]*Tamp[iq-iqmin];
	  /*if(i==watchdog)
          {
	    cout << "Tk[watchdog][" << iq-iqmin << "][" << ica << "]=" << Tk[i][iq-iqmin][ica] << endl;
	  }*/
	  idx++;
	}
      }
    }
  }
  
  
  for(int i=0; i<nv; i++)
  {
    Tmean[i]=0;
    for(int j=0; j<nv; j++)
      Tprod[i][j]=0;    
  }
  
  double  dn=0;  
  double qmean=0.6;
  double s2sigq = sqrt(2.)*0.12;
  
  double qprev = iqmin*qstep;
  
  for(int iq=iqmin; iq<=iqmax; iq++) 
  {    
    double q = (double(iq)+0.5)*qstep; // upper integration limit
    if(q>qmax) q=qmax; // last step
    
    double w = erf((q-qmean)/s2sigq)-erf((qprev-qmean)/s2sigq);
    
    for(int i=0; i<nv; i++) 
    {
     for(int ica=0; ica<nca; ica++)
     {
      double ww=w;
      if(ica==0 || ica==nca-1) ww /= 2.;
      Tmean[i] += Tk[i][iq-iqmin][ica]*ww;
      for(int j=i; j<nv; j++)
      {
       Tprod[i][j] += Tk[i][iq-iqmin][ica]*Tk[j][iq-iqmin][ica]*ww;
      }
     }
    }
    dn += w*(nca-1.);
    qprev = q;
  }
  
  //cout << "watchdog variance " << Tprod[watchdog][watchdog]/dn << "-" << pow(Tmean[watchdog]/dn,2) << "=" << Tprod[watchdog][watchdog]/dn-pow(Tmean[watchdog]/dn,2) << endl;
  
  //cout << "Tmean[0] = " << Tmean[0] << "/" << dn << endl; 
  
  for(int i=0; i<nv; i++)
  {
    Tmean[i] /= dn;
  }  
  
  //cout << "Tprod[0][0] = " << Tprod[0][0]/dn << "-" << Tmean[0]*Tmean[0] << endl;
  
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
     image[idx]=Tprod[i][j]; // in units of  [ R_s*rho_0(q=1) ]^2 because we normalized it to that
     idx++;
    }
  
  pFits->pHDU().write(1,nv*nv,image);
  
  //cout << "Tmean[watchdog]=" << Tmean[watchdog] << endl;
  //cout << "Tvar[watchdog]=" << Tprod[watchdog][watchdog] << endl;
  
  
  return 0;
}

