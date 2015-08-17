#include "cosmology.h"
#include <iostream>

using namespace std;

int main(int atoi, char **argv)
{
	if(atoi!=2) {
		cerr << "syntax: " << argv[0] << " [source redshift distribution file]" << endl;
	}

        ofstream cosmo("cosmo.par");
        ofstream nofz("nofz.par");

        // cosmo.par

	cosmo << "Omega_m        " << OmegaM << endl;
	cosmo << "Omega_de       " << OmegaL << endl;
	cosmo << "w0_de          " << -1.0 << endl;
	cosmo << "w1_de          " <<  0.0 << endl;
	cosmo << "h_100          " << h << endl;
	cosmo << "Omega_b        " << OmegaB << endl;
	cosmo << "Omega_nu_mass  " << 0.0 << endl;
	cosmo << "Neff_nu_mass   " << 0.0 << endl;
	cosmo << "normalization  " << sigma8 << endl;
	cosmo << "n_spec         " << ns << endl;
        cosmo << endl;
        cosmo << "snonlinear     smith03_revised" << endl;
        cosmo << "stransfer      eisenhu_osc" << endl;
        cosmo << "sgrowth        growth_de" << endl;
        cosmo << "sde_param      linder" << endl;
        cosmo << "normmode       0" << endl;
        cosmo << "a_min          0.2" << endl;
          
        // nofz.par

 	nofz << "Nzbin           1" << endl;
        nofz << "snzmode         nz_read_from_files" << endl;
        nofz << "nzfile          " << argv[1] << endl;

	return 0;
}

