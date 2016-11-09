// getmodel: write model covariance (in units of kappa) for halo at given mass and redshift to file

#include "model/covariance.h"
#include "nicaea_pz.h"


int main(int argc, char **argv)
{

	if(argc < 6 || argc > 11 ) {
		cerr << "syntax: " << argv[0] << " [output filename] [zlens] [annuli prefix] [zsource or p(z) prefix] [M200m/(h_^{-1}Msol)]"
		                              << " ([cconc=1] [ccorr=5.0] [cell=3.7] [coff=0] [clss=0])" << endl;
		return 1;
	}

	// (1) read command line arguments

	double zlens = atof(argv[2]);
	double m200m = atof(argv[5]);
	int log10mass100 = 100.*log10(m200m);
        string ap(argv[3]);
	
	double beta=0.;
	{
	  double zsource=atof(argv[4]);
	  if(zsource>zlens) { // get beta from fixed zsource
	    beta=angularDiameterDistance(zlens,zsource,1000)/angularDiameterDistance(0,zsource,1000);
	  } else {
	    beta=beta_from_pz(zlens, string(argv[4])+".tab"); 
	  }
	}
	
	double cconc = 1.;  // setting cconc=0 would be the right choice for fitting both c and M
	if(argc>6) cconc = atof(argv[6]);
	double ccorr = 5.;
	if(argc>7) ccorr = atof(argv[7]);
	double cell  = 3.7;
	if(argc>8) cell  = atof(argv[8]);
	double coff  = 0.;
	if(argc>9) coff  = atof(argv[9]);
        double clss  = 0.;
	if(argc>10) clss = atof(argv[10]);
	
	cerr << "# calculating intrinsic covariance for M200m=" << m200m << "h^-1 Msol, z=" << zlens << endl;

	// (2) initialize model

	CovarianceModel model(zlens, beta, log10mass100, log10mass100+1, 
			       ccorr, cconc, cell, coff, 0., 0.,
				"model/corrh_"+ap+"_", 			    // string prefixes of model covariance
				"model/conc/conc_"+ap+"_m",	
				"model/ell/ell_"+ap+"_m",	
				"model/off/off_"+ap+"_m");
       
	SymMatrix<double> cov = model.cov(m200m);

        // (3) write to file

	Covariance::writeMatrixToFITS(cov, argv[1]);

	return 0;
}

