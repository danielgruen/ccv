#include <iostream>
#include "cosmology.h"

using namespace std;

double lambda_at_m_z(double m200m, double z)
{
  double lambda = 0.; // this is not an implementation yet :)
  cerr << "I am doing nothing, but at least I know that dn/dM=" << halomodel::dndM(m200m, z) << endl;
  return lambda;
}


int main(int argc, char **argv)
{
    if(argc != 3) {
      cout << "syntax: " << argv[0] << " [M200m at h=0.7, solar masses] [z]" << endl;
      cout << "e.g. ./src/lambda_at_m_z 1e14 0.24533" << endl;
      return 1;
    }
    
    double m200m=atof(argv[1]);
    double z=atof(argv[2]);
    
    cout << "m200m=" << m200m << " z=" << z << " lambda=" << lambda_at_m_z(m200m, z) << endl;

	return 0;
}

