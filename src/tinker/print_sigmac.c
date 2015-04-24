#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <math.h>
#include "header.h"


int main(int argc, char **argv)
{
  if(argc==1)
    endrun("./print_sigmac bat_file");

  read_parameter_file(argv[1]);

  printf("%e %e %e %e %e\n",OMEGA_M,1.-OMEGA_M,OMEGA_B,SIGMA_8,SPECTRAL_INDX); 
  double norm = SIGMA_8/sigmac(8.);
  double logM;
  for(logM=7; logM<=16.04; logM+=0.05)
    {
      double m = pow(10,logM);
      double x = norm*sigmac(pow(3.0*m/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0));
      printf("%e %e %e\n",pow(3.0*m/(4.0*PI*OMEGA_M*RHO_CRIT),1./3.),logM,x);
    }
    
  return 0;
  
}
