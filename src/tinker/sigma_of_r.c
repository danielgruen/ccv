#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <math.h>
#include "header.h"


int main(int argc, char **argv)
{
  printf("%e\n",SIGMA_8); 
  double norm = SIGMA_8/sigmac(8.);
  double x = norm*sigmac(atof(argv[1]));
  printf("%e %e\n",atof(argv[1]),x);
    
  return 0;
  
}
