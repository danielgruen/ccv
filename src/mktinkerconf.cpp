#include "cosmology.h"
#include <iostream>

using namespace std;

int main(int atoi, char **argv)
{
	if(atoi!=2) {
		cerr << "syntax: " << argv[0] << " [redshift]" << endl;
	}

	cout << "OMEGA_M        " << OmegaM << endl;
	cout << "SIGMA_8        " << sigma8 << endl;
	cout << "RHO_CRIT       " << rho0crit << endl;
	cout << "SPECTRAL_INDX  " << ns << endl;
	cout << "HUBBLE         " << h << endl;
	cout << "OMEGA_B        " << OmegaB << endl;
        cout << "DELTA_CRIT      1.686" << endl;
	cout << "ITRANS          5" << endl;
	cout << "TF_file         transfunc.WMAP3" << endl;
	cout << "DELTA_HALO      200" << endl;
	cout << "REDSHIFT       " << argv[1] << endl;
        cout << "root_filename   tmp" << endl;
        
	return 0;
}

