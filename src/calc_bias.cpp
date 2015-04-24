// print nu and bias as a function of mass at fixed z

#include <CCfits/CCfits>
using namespace CCfits;

#include "corrh/corrh.h"
#include <iostream>
#include <cstdlib>

using namespace std;

int main(int argc, char **argv)
{
  double z=atof(argv[1]);
  
  cout << OmegaM << " " << OmegaL << " " << OmegaB << " " << sigma8 << " " << ns << " " << z << endl;
 
  for(double logM=7; logM<=16; logM+=0.01)
  {
    double M = pow(10,logM);
    cout << logM << " " << halomodel::nu(M,z) << " " << halomodel::bias_tinker(M,z) << endl;
  }
  
  return 0;
}

