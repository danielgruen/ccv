#ifndef _COVARIANCE_H_
#define _COVARIANCE_H_

#include <CCfits/CCfits>

using namespace CCfits;

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

#include <TMV.h>
#include <TMV_Sym.h>

#ifndef _COVARIANCE_HELPERS_H_
#include "../cosmology.h"
#include "../enfw/enfw.h"
#endif

using namespace tmv;
using namespace std;

struct CovarianceParameters {
  // stuff one might need to decide what the covariance looks like
  double m200m;
  double nu;
  double bias;
};

CovarianceParameters covarianceParameters(double m200m, double z) {
  CovarianceParameters p = {m200m, halomodel::nu(m200m, z), halomodel::bias_tinker(m200m, z)};
  return p;
}

CovarianceParameters pnull = {0,0,0};

class Covariance { 
//
// a class that represents an abstract covariance
//
  public:
    Covariance() { }
    
    virtual SymMatrix<double> cov(CovarianceParameters &p=pnull) = 0;   
    // pure virtual, must be implemented in derived class
    
    virtual DiagMatrix<double> covDiag(CovarianceParameters &p=pnull) = 0;   
    // pure virtual, must be implemented in derived class
    
    virtual Vector<double> var(CovarianceParameters &p=pnull) = 0;   
    // pure virtual, must be implemented in derived class
        
    virtual double twolnlikelihood(Vector<double> &residual, CovarianceParameters p) { 
    // there could be more efficient implementations depending on your covariance model, so maybe implement in derived classes
      SymMatrix<double> c = cov(p);
      return residual*(residual/c) + c.logDet(); 
    }
    virtual double twolnlikelihooddiag(Vector<double> &residual, CovarianceParameters p) { 
    // there could be more efficient implementations depending on your covariance model, so maybe implement in derived classes
      DiagMatrix<double> d = covDiag(p);
      return residual*(residual/d) + d.logDet();
    } 
    
    static void writeMatrixToFITS(const SymMatrix<double> &m, const string &filename)
    {
	FITS *pFits = 0;
	long naxis    =   2;
	assert(m.nrows()==m.ncols());
	long naxes[2] = { m.nrows(), m.nrows() };

	try
	{
	  pFits = new FITS("!"+filename, DOUBLE_IMG , naxis , naxes );
	}
	catch (FITS::CantCreate)
	{
	  cerr << "error, can't create fits file" << endl;
	  return;
	}

	valarray <double> image(0., m.nrows()*m.nrows()); // we're just going to jam our data in here ...
	int idx=0;
	for(int i=0; i<m.nrows(); i++)
	  for(int j=0; j<m.nrows(); j++)
	  {
	  image[idx]=m(i,j);
	  idx++;
	  }

	pFits->pHDU().write(1,m.nrows()*m.nrows(),image); 
    }
    
    bool isPosDef(CovarianceParameters &p=pnull) // good to know
    {
      try
      {
	SymMatrix<double> m = cov(p);
	m.view().chd();
      }
      catch (tmv::NonPosDef) {
	return false;
      }
      return true;
    }
    
  protected:
    static double max(double a, double b)
    {
     return (a>b)?a:b; 
    }
    
};


class CovarianceMatrix : public Covariance { 
// 
// a covariance as described by a matrix
//
  public:
    CovarianceMatrix() : Covariance() {
      c=0; d=0; vd=0;
    }
    
    ~CovarianceMatrix() {
     delete c;
     delete d;
     delete vd;
    }
 
    SymMatrix<double> cov(CovarianceParameters &p=pnull)
    {
     return *c; 
    }
    DiagMatrix<double> covDiag(CovarianceParameters &p=pnull)
    {
     return *d; 
    }
    Vector<double> var(CovarianceParameters &p=pnull)
    {
     return *vd; 
    }
    
    void rescaleCovariance(double scale) {
     (*c)  *= scale; 
     (*d)  *= scale;
     (*vd) *= scale;
    }
    
    void rescaleCovariance(vector<double> scale) {
     assert(scale.size()==c->colsize());
     
     for(int i=0; i<c->colsize(); i++)
     {
       (*d)(i)    = (*d)(i)*scale[i]*scale[i];
       (*vd)(i)   = (*vd)(i)*scale[i]*scale[i];
       for(int j=i; j<c->colsize(); j++)
       {
         (*c)(i,j) = (*c)(i,j)*scale[i]*scale[j];
       }
     }
    }
    
    void setMatrix(SymMatrix<double> &dc) {
     delete c;
     delete d;
     delete vd;
     c = new SymMatrix<double>(dc);
     d = new DiagMatrix<double>(DiagMatrixViewOf(dc.diag())); 
     vd= new Vector<double>(dc.diag());
    }

    void subtract(CovarianceMatrix &dc) {
     *c  = (*c)-  dc.cov();
     *d  = (*d)-  dc.covDiag();
     *vd = (*vd)- dc.var();
    }
    
    void project(const MatrixView<double> P)
    // project to a new basis described by rows of P
    // C --> P C P^T in that basis
    {
     Matrix<double> buf(P.colsize(),P.colsize());
     buf = P * (*c) * P.transpose();
     
     //cout << "\n\n\nI have rotated a covariance matrix; is it still symmetric? " << endl;
     //cout << buf << endl;
     
     //cout << "lets try make a symmatrix copy of it " << endl;
     SymMatrix<double> bufs(buf);
     //cout << bufs << endl;
     
     setMatrix(bufs);
    }
    
    void addFromWeight(vector<double> weight) {
    // add variance according to a 1/var weight vector to diagonal
     assert(weight.size()==c->colsize());
     for(int i=0; i<weight.size(); i++) {
       (*c)(i,i) = (*c)(i,i)+1./weight[i];
       (*d)(i)   = (*d)(i)+1./weight[i];
       (*vd)(i)  = (*vd)(i)+1./weight[i];
     }
    }

  protected:
    SymMatrix<double>  *c;
    DiagMatrix<double> *d;
    Vector<double>    *vd;
};


class CovarianceSum : public Covariance {

  public:
    CovarianceSum(Covariance *dc1, Covariance *dc2) : Covariance(), c1(dc1), c2(dc2) {
    
    }
    
    SymMatrix<double> cov(CovarianceParameters &p) {
      SymMatrix<double> c(c1->cov(p)+c2->cov(p));
      return c;
    }
    
    DiagMatrix<double> covDiag(CovarianceParameters &p){  
      DiagMatrix<double> c(c1->covDiag(p)+c2->covDiag(p));
      return c;
    }
    
    Vector<double> var(CovarianceParameters &p){
      Vector<double> c(c1->var(p)+c2->var(p));
      return c;
    }

  private:
    Covariance *c1,*c2;
};


class FITSCovarianceMatrix : public CovarianceMatrix {
// 
// a covariance matrix as described in a FITS file
//
  public:
    FITSCovarianceMatrix(string filename) : CovarianceMatrix() {
        
        valarray<double>  contents;  
	long ax1,ax2;
        try {
	  auto_ptr<FITS> pInfile(new FITS(filename.c_str(),Read,true));
	  PHDU& image = pInfile->pHDU(); 
	  image.read(contents);
	  ax1 = image.axis(0);
	  ax2 = image.axis(1);
	  assert(ax1==ax2);
	} catch(...) {
	  cerr << "# cannot open " << filename << "; giving up. You should probably add that redshift to your Makefile and re-run \'make\' -- or maybe you are running the code from somewhere other than the base directory." << endl;
	  throw;
	}
	    

	
	Nb=ax1;
	c = new SymMatrix<double>(Nb);
	d = new DiagMatrix<double>(Nb);
	vd= new Vector<double>(Nb);
	
	int idx=0;
	for(int i=0; i<Nb; i++)
	{
	  for(int j=0; j<Nb; j++)
	  {
	    (*c)(i,j)=contents[idx];
	    if(i==j){
	      (*d)(i)=(*vd)(i)=contents[idx];
	    }
	    idx++;
	  }
	}
    }
    
    ~FITSCovarianceMatrix() {
      delete c; c=0;
      delete d; d=0;
      delete vd; vd=0;
    }
    
    int Nb;
};


class CovarianceMassArray : public Covariance {
//
// covariance as a function of mass, interpolated between discrete values for which we have a FITSCovarianceMatrix
//
  public:
    CovarianceMassArray(string prefix, double z, int immin, int immax, int imstep) : Covariance() {

        cerr << "# initializing CovarianceMassArray with prefix " << prefix << ", z=" << z << ", im=" << immin << "," << immax << "," << imstep << endl;
 
	logmassmin = immin/100.;
	logmassmax = immax/100.;
	logmassstep= imstep/100.;
	
	assert(immax>immin);
	assert(imstep>0);
	assert(z<1);
      
	for(int im=immin; im<=immax; im+=imstep)
	{
	    std::ostringstream ss;
	    ss << im << "_0." << int(z*10000000.+0.5) << ".fits";
	    FITSCovarianceMatrix *tmp = new FITSCovarianceMatrix(prefix+ss.str());
	    
	    matrices.push_back(tmp);
	    logmasses.push_back(im/100.);
	}
	
    }
    
    ~CovarianceMassArray() {
	for(int i=0; i<matrices.size(); i++)
	{
	 delete matrices[i]; 
	}
    }
      
    SymMatrix<double> cov(CovarianceParameters &p)
    {
	double logM = log10(p.m200m);
	
	int ibelow    = int((logM-logmassmin)/logmassstep);
	assert(ibelow>=0);
	assert(ibelow<matrices.size()-1);
	
	double wbelow = (logmasses[ibelow+1]-logM)/logmassstep;
	
	SymMatrix<double> tmpmatrix = (wbelow*(matrices[ibelow]->cov(p))+(1.-wbelow)*(matrices[ibelow+1]->cov(p))); 
	return tmpmatrix;
    }

      
    DiagMatrix<double> covDiag(CovarianceParameters &p)
    {
	double logM = log10(p.m200m);
	
	int ibelow    = int((logM-logmassmin)/logmassstep);
	assert(ibelow>=0);
	assert(ibelow<matrices.size()-1);
	
	double wbelow = (logmasses[ibelow+1]-logM)/logmassstep;
	
	DiagMatrix<double> tmpmatrix = (wbelow*(matrices[ibelow]->covDiag(p))+(1.-wbelow)*(matrices[ibelow+1]->covDiag(p))); 
	return tmpmatrix;
    }
 
    Vector<double> var(CovarianceParameters &p)
    {
	double logM = log10(p.m200m);
	
	int ibelow    = int((logM-logmassmin)/logmassstep);
	assert(ibelow>=0);
	assert(ibelow<matrices.size()-1);
	
	double wbelow = (logmasses[ibelow+1]-logM)/logmassstep;
	
	return (wbelow*(matrices[ibelow]->var(p))+(1.-wbelow)*(matrices[ibelow+1]->var(p))); 
    }

    void rescaleCovariance(double (*f)(double,double,double,double), double p1, double p2, double p3)
    {
	for(int i=0; i<matrices.size(); i++)
	{
	   matrices[i]->rescaleCovariance(f(pow10(logmasses[i]),p1,p2,p3)); 
	}
    }

    void rescaleCovariance(double f)
    {
	for(int i=0; i<matrices.size(); i++)
	{
	   matrices[i]->rescaleCovariance(f); 
	}
    }
    
    void rescaleCovariance(vector<double> f)
    {
	for(int i=0; i<matrices.size(); i++)
	{
	   matrices[i]->rescaleCovariance(f); 
	}    
    }
    
    void project(const MatrixView<double> P)
    // project to a new basis described by rows of P
    // C --> P C P^T in that basis
    {
      
	for(int i=0; i<matrices.size(); i++)
	{
	 matrices[i]->project(P); 
	}
    }
    
    vector<FITSCovarianceMatrix*> matrices;
  private:
    vector<double> logmasses;
    double logmassmin,logmassmax,logmassstep;
};

class CovarianceModel : public Covariance {

  public:
    CovarianceModel(double dzlens, double beta=1.0, 		    	    // lens and source redshift
                    int log10massmin100=1300, int log10massmax100=1600,	    // minimum and maximum log10(M200m/[h^-1Msol])*100
		    double dccorr=5.0, double dcconc=1.0, double dcell=3.7, // semi-analytic re-scaling factors
		    double dcoff=0.0,
                    double dccorr1=0., double dcell1=0.,		    // nu dependence of re-scaling factors
                    string scorr="model/corrh_", 			    // string prefixes of model covariance
                    string sconc="model/conc/conc_m",	
                    string  sell="model/ell/ell_m",	
                    string  soff="model/off/off_m"
                   ) :
                   zlens(dzlens), ccorr(dccorr), cconc(dcconc), cell(dcell), ccorr1(dccorr1), cell1(dcell1), coff(dcoff)
    {
        double DLen=angularDiameterDistance(0,dzlens,1000);
        sigmacrit = ckms*ckms/4./M_PI/Gs/beta/DLen; // h Msol / Mpc^2;

	// (1) correlated haloes; mass dependence comes from bias of cluster halo only
        stringstream ss;
        assert(zlens<1.);
	ss << scorr << "0." << int(zlens*10000000.+0.5) << ".fits";
	Ccorr = new FITSCovarianceMatrix(ss.str()); 
	Ccorr->rescaleCovariance(1./(sigmacrit*sigmacrit));
	
	// (3) concentration scatter
	Cconc = new CovarianceMassArray(sconc, zlens, log10massmin100, log10massmax100, 1);
	Cconc->rescaleCovariance(rescaleCconc, zlens, sigmacrit, DLen);
	
	// (4) halo ellipticity
	Cell  = new CovarianceMassArray(sell, zlens, log10massmin100, log10massmax100, 1);
        Cell->rescaleCovariance(rescaleCell, zlens, sigmacrit, DLen);
	
	// (4) off-centring
	Coff  = new CovarianceMassArray(soff, zlens, log10massmin100, log10massmax100, 1);
        Coff->rescaleCovariance(rescaleCoff, zlens, sigmacrit, DLen);
    }

    ~CovarianceModel() 
    {
      delete Ccorr; delete Cconc; delete Cell; delete Coff;
    }
     
    void project(const MatrixView<double> P)
    // project to a new basis described by rows of P
    // C --> P C P^T in that basis
    {
      if(Ccorr)  Ccorr->project(P);
      if(Cconc)  Cconc->project(P);
      if(Cell)   Cell->project(P);
      if(Coff)   Coff->project(P);
    }
   
    static double rescaleCconc(double m200m, double Szlens, double Ssigmacrit, double SDLen) 
    // re-scaling of FITSCovarianceMatrix for concentration scatter at m200m [hinv Msol]
    {
      double r200m = halomodel::r200mRadius(m200m/h, Szlens)*h/SDLen*arcminperrad; // arcmin
      return pow(pow(1.+Szlens,3)*r200m*SDLen*radperarcmin / Ssigmacrit,2);
    }    
    static double rescaleCell(double m200m, double Szlens, double Ssigmacrit, double SDLen) 
    // re-scaling of FITSCovarianceMatrix for ellipticity scatter at m200m [hinv Msol]
    {
      double r200m = halomodel::r200mRadius(m200m/h, Szlens)*h/SDLen*arcminperrad; // arcmin
      double c200m = halomodel::cmz_200m_duffy(m200m,Szlens);
      return pow(nfw_rho0_200m(c200m, Szlens)/h/h * r200m/c200m *SDLen*radperarcmin / Ssigmacrit,2);
    }
    static double rescaleCoff(double m200m, double Szlens, double Ssigmacrit, double SDLen) 
    // re-scaling of FITSCovarianceMatrix for off-centring at m200m [hinv Msol]
    {
      double r200m = halomodel::r200mRadius(m200m/h, Szlens)*h/SDLen*arcminperrad; // arcmin
      double c200m = halomodel::cmz_200m_duffy(m200m,Szlens);
      double rho0 = nfw_rho0_200m(c200m, 0.)/h/h; // h^2 Msol/R_S[Mpc]^3
      return pow(rho0/c200m*pow(1.+Szlens,3)*r200m*SDLen*radperarcmin / Ssigmacrit,2);
    }
 
    SymMatrix<double> cov(CovarianceParameters &p)
    {
      return  max(0,ccorr+ccorr1*(p.nu-2.5))*p.bias*Ccorr->cov() +
              cconc*Cconc->cov(p) +
              max(0,cell+cell1*(p.nu-2.8))*Cell->cov(p) +
	      coff*Coff->cov(p);
    }

    DiagMatrix<double> covDiag(CovarianceParameters &p)
    {
      return  max(0,ccorr+ccorr1*(p.nu-2.5))*p.bias*Ccorr->covDiag() +
              cconc*Cconc->covDiag(p) +
              max(0,cell+cell1*(p.nu-2.8))*Cell->covDiag(p) +
              coff*Coff->covDiag(p);
    }
    
    Vector<double> var(CovarianceParameters &p)
    {
      return  max(0,ccorr+ccorr1*(p.nu-2.5))*p.bias*Ccorr->var() +
              cconc*Cconc->var(p) +
              max(0,cell+cell1*(p.nu-2.8))*Cell->var(p) +
              coff*Coff->var(p);
    }

    // convenience
    SymMatrix<double> cov(double m200m) // in units of h^{-1} Msol
    {
      CovarianceParameters p = covarianceParameters(m200m, zlens);
      return cov(p);
    }
    DiagMatrix<double> covDiag(double m200m) // in units of h^{-1} Msol
    {
      CovarianceParameters p = covarianceParameters(m200m, zlens);
      return covDiag(p);
    }
    Vector<double> var(double m200m) // in units of h^{-1} Msol
    {
      CovarianceParameters p = covarianceParameters(m200m, zlens);
      return var(p);
    }
    
    void rescaleCovariance(double scale) {
      Ccorr->rescaleCovariance(scale);
      Cconc->rescaleCovariance(scale);
      Cell->rescaleCovariance(scale);
      Coff->rescaleCovariance(scale);
    }
    
    void rescaleCovariance(vector<double> &scale) {
      Ccorr->rescaleCovariance(scale);
      Cconc->rescaleCovariance(scale);
      Cell->rescaleCovariance(scale);
      Coff->rescaleCovariance(scale);
    }

    double zlens, sigmacrit;
    double ccorr, cconc, cell, coff;
    double ccorr1, cell1;

    FITSCovarianceMatrix *Ccorr;
    CovarianceMassArray  *Cconc, *Cell, *Coff;
 
};


static double ranf()
{
 return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}
static double gaussian_random(float m, float s) /* normal random variate generator */
{                                               /* mean m, standard deviation s */
        double x1, x2, w, y1;
        static double y2;
        static int use_last = 0;

        if (use_last)                   /* use value from previous call */
        {
                y1 = y2;
                use_last = 0;
        }
        else
        {
                do {
                        x1 = 2.0 * ranf() - 1.0;
                        x2 = 2.0 * ranf() - 1.0;
                        w = x1 * x1 + x2 * x2;
                } while ( w >= 1.0 );

                w = sqrt( (-2.0 * log( w ) ) / w );
                y1 = x1 * w;
                y2 = x2 * w;
                use_last = 1;
        }

        return( m + y1 * s );
}

vector<double> add_noise(vector<double> &data, Covariance &cov)
// returns the sum of data and a Gaussian noise realization according to cov
{
    Matrix<double> U = cov.cov().chd().getL();

    Vector<double> r(data.size());
    for(int j=0; j<data.size(); j++)
       r(j) = gaussian_random(0,1);
      
    Vector<double> noise = r*U.transpose();
    
    vector<double> out(data);
    for(int i=0; i<data.size(); i++)
       out[i] += noise(i);
       
    return out;
}


#endif

