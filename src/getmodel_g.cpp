// getmodel_g: write model covariance (in units of gamma) for halo at given mass and redshift to file

#include "model/covariance.h"

int main(int argc, char **argv)
{

	if(argc < 4 || argc > 8 ) {
		cerr << "syntax: " << argv[0] << " [output filename] [zlens] [M200m/(h_^{-1}Msol)] ([zsource=1] [cconc=1] [ccorr=5.0] [cell=3.7])" << endl;
		return 1;
	}

	// (1) read command line arguments

	double zlens = atof(argv[2]);
	double m200m = atof(argv[3]);
	int log10mass100 = 100.*log10(m200m);

	double zsource = 1.;
	if(argc>4) zsource = atof(argv[4]);
	double cconc = 1.;  // setting cconc=0 would be the right choice for fitting both c and M
	if(argc>5) cconc = atof(argv[5]);
	double ccorr = 5.;
	if(argc>6) ccorr = atof(argv[6]);
	double cell  = 3.7;
	if(argc>7) cell  = atof(argv[7]);

	cerr << "# calculating intrinsic covariance for M200m=" << m200m << "h^-1 Msol, z=" << zlens << endl;

	// (2) initialize model

	CovarianceModel model(zlens, zsource, log10mass100, log10mass100+1, ccorr, cconc, cell, 0., 0., 
			      "model/gamma/corrh_", // string prefixes of model covariance for gamma
			      "model/gamma/conc/conc_m",	
			      "model/gamma/ell/ell_m"); 
		// it's ok to just say CovarianceModel model(zlens, zsource) if you want the default, kappa version
	SymMatrix<double> cov = model.cov(m200m);

        // (3) write to file

	Covariance::writeMatrixToFITS(cov, argv[1]);

	return 0;
}

