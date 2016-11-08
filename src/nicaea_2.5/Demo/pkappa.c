/* ============================================================ *
 * lensingdemo.c						*
 * Martin Kilbinger 2008-2010					*
 * ============================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

#include "errorlist.h"
#include "io.h"
#include "maths.h"
#include "cosmo.h"
#include "coyote.h"
#include "cmb_bao.h"
#include "lensing.h"
#include "nofz.h"


#define ELL_MIN 0.01
#define ELL_MAX 100000.0
#define NELL    1002
void usage(int ex)
{
   fprintf(stderr, "Usage: pkappa [OPTIONS]\n");
   fprintf(stderr, "OPTIONS\n");
   fprintf(stderr, "  -L 'ELL_MIN ELL_MAX NELL' NELL Fourier modes between ELL_MIN and ELL_MAX (default: %g %g %d)\n",
      ELL_MIN, ELL_MAX, NELL);
   fprintf(stderr, "  -h            This message\n");

   if (ex >= 0) exit(ex);
}


int main(int argc, char** argv)
{
   cosmo_lens* model;
   double ell, pk, pg;
   int i_bin, j_bin, Nzbin, c;
   error *myerr = NULL, **err;
   char *ell_info, *substr;
   FILE *F;
   size_t Nell;
   double ell_min, ell_max, ell_fac;


   printf("# lensingdemo (M. Kilbinger, 2006-2012)\n");


   err = &myerr;


   ell_info = NULL;
   while (1) {

      static struct option long_option[] = {
         {"", required_argument, 0, 'L'},
         {0, 0, 0, 0}
      };

      int option_index = 0;

      c = getopt_long(argc, argv, "L:h", long_option, &option_index);
      switch (c) {
         case 'L':
            ell_info = optarg;
            break;
         case 'h' :
            usage(0);
         case -1 :
            goto end_arg;
         default :
            usage(1);
      }

   }
end_arg:

   /* Setting up cosmology */

   F = fopen_err("cosmo_lens.par", "r", err);
   quitOnError(*err, __LINE__, stderr);
   read_cosmological_parameters_lens(&model, F, err);
   quitOnError(*err, __LINE__, stderr);
   fclose(F);

   Nzbin = model->redshift->Nzbin;

   printf("# Cosmological parameters:\n");
   dump_param_lens(model, stdout, 1, err);
   quitOnError(*err, __LINE__, stderr);

   for (i_bin=0; i_bin<Nzbin; i_bin++) {
      printf("Mean redshift (bin %d) = %.3f\n", i_bin, zmean(model->redshift, i_bin, err));
      quitOnError(*err, __LINE__, stderr);
   }


   /* Density power spectrum */

   if (model->cosmo->nonlinear == coyote10) {
      set_H0_Coyote(model->cosmo, err);
      quitOnError(*err, __LINE__, stderr);
      printf("Coyote10: Hubble parameter set to h = %g\n", model->cosmo->h_100);
   }

   /* Convergence power spectrum */

   if (ell_info != NULL) {
      substr    = strsep(&ell_info, " "); ell_min = atof(substr);
      substr    = strsep(&ell_info, " "); ell_max = atof(substr);
      substr    = strsep(&ell_info, " "); Nell    = atoi(substr);
   } else {
      ell_min   = ELL_MIN;
      ell_max   = ELL_MAX;
      Nell      = NELL;   
   }
   ell_fac   = pow( ell_max / ell_min, 1.0 / ((double)Nell - 1));

   printf("Calculating convergence power spectrum, writing file P_kappa\n");
   F = fopen("P_kappa", "w");
   dump_param_lens(model, F, 1, err);
   quitOnError(*err, __LINE__, stderr);

   fprintf(F, "# Shear cross-power spectra for %d z-bin%s\n", Nzbin, Nzbin==1?"":"s");
   fprintf(F, "# l            ");
   for (i_bin=0; i_bin<Nzbin; i_bin++) {
      for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
	    if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
	    if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
	    fprintf(F, "P_%s^%d%d(l)    ", model->reduced==reduced_K10?"t":"k", i_bin, j_bin);
      }
   }
   fprintf(F, "   ");
   for (i_bin=0; i_bin<Nzbin; i_bin++) {
      for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
         if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
         if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
         if (model->reduced!=reduced_none) {
            fprintf(F, "P_g1^%d%d(l)   ", i_bin, j_bin);
         }
      }
   }
   fprintf(F, "\n");

   for (ell=ell_min; ell<=ell_max*1.001; ell*=ell_fac) {
      fprintf(F, "%e  ", ell);
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            pk = Pshear(model, ell, i_bin, j_bin, err);
            if (getErrorValue(*err)==reduced_fourier_limit) {
               purgeError(err);
            }
            quitOnError(*err, __LINE__, stderr);
            fprintf(F, " %e", pk);
         }
      }
      fprintf(F, "   ");
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            if (model->reduced!=reduced_none) {
               pg = Pg1(model, ell, i_bin, j_bin, err);
               if (getErrorValue(*err)==reduced_fourier_limit) {
                  purgeError(err);
               }
               quitOnError(*err, __LINE__, stderr);
               fprintf(F, " %e", pg);
            } else {
               pg = 0.0;
            }
         }
      }
      fprintf(F, "\n");
   }
   fclose(F);


   printf("Done\n");

   free_parameters_lens(&model);

   return 0;
}
#undef ELL_MIN
#undef ELL_MAX
#undef NELL
