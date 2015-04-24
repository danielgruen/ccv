#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "header.h"


void output_halo_mass_function()
{
  int k,nr=200;
  double x,dlogm,m,mvir,cdelta,mmin=1e7,mmax=1e16,delta_vir;
  FILE *fp;
  char aa[100];

  fprintf(stderr,"\n\nCALCULATING HALO MASS FUNCTION.\n");
  fprintf(stderr,    "-------------------------------\n\n");

  sprintf(aa,"%s.dndM",Task.root_filename);
  fp = fopen(aa,"w");

  dlogm = (log(mmax) - log(mmin))/(nr-1);

  fprintf(fp,"%e %e %e %e %e %e\n",OMEGA_M,1.-OMEGA_M,OMEGA_B,SIGMA_8/growthfactor(REDSHIFT),SPECTRAL_INDX,REDSHIFT);

  for(k=0;k<nr;++k)
    {
      m = exp(k*dlogm)*mmin;
      if(m<1.3e14 || m>1.9e14)
       x = halo_mass_function(m);
      else { // interpolate numerically unstable region
       double m0 = 1.2e14;
       double m1 = 2.0e14;
       double x0 = halo_mass_function(m0);
       double x1 = halo_mass_function(m1);
       
       double lm = log(m);
       double lm0 = log(m0);
       double lm1 = log(m1);
       double dlm = lm1-lm0;
       double lx0 = log(x0);
       double lx1 = log(x1);
       
       x = exp(lx0*(lm1-lm)/dlm+lx1*(lm-lm0)/dlm);
      }
      fprintf(fp,"%e %e\n",m,x);
    }
  fclose(fp);
}
