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


#define TH_MIN  0.5
#define TH_FAC  1.1
#define NTHETA  71
#define ELL_MIN 0.01
#define ELL_MAX 10000.0
#define NELL    50
void usage(int ex)
{
   fprintf(stderr, "Usage: lensingdemo [OPTIONS]\n");
   fprintf(stderr, "OPTIONS\n");
   fprintf(stderr, "  -t FNAME               File containing angular scales in arcmin\n");
   fprintf(stderr, "                          (default: %d bins between %g and %g arcmin)\n",
	   NTHETA, TH_MIN, TH_MIN * pow(TH_FAC, NTHETA-1));
   fprintf(stderr, "  -T 'TH_MIN TH_MAX NTH' NTH angular bins between TH_MIN and TH_MAX arcmin\n");
   fprintf(stderr, "  -L 'ELL_MIN ELL_MAX NELL' NELL Fourier modes between ELL_MIN and ELL_MAX (default: %g %g %d)\n",
      ELL_MIN, ELL_MAX, NELL);
   fprintf(stderr, "  -h            This message\n");

   if (ex >= 0) exit(ex);
}


int main(int argc, char** argv)
{
   cosmo_lens* model;
   double k, ell, theta, a, pk, pg, f, z, cumG, Ga, ww;
   int i_bin, j_bin, Nzbin, c;
   error *myerr = NULL, **err;
   char *theta_fname, *theta_info, *ell_info, *substr;
   FILE *F;
   size_t Nrec, Ntheta, Nell;
   double *theta_list, Theta_min, Theta_max, Theta_fac, ell_min, ell_max, ell_fac;


   printf("# lensingdemo (M. Kilbinger, 2006-2012)\n");


   err = &myerr;


   theta_fname = theta_info = ell_info = NULL;
   while (1) {

      static struct option long_option[] = {
         {"", required_argument, 0, 't'},
         {"", required_argument, 0, 'T'},
         {"", required_argument, 0, 'L'},
         {0, 0, 0, 0}
      };

      int option_index = 0;

      c = getopt_long(argc, argv, "t:T:L:h", long_option, &option_index);
      switch (c) {
         case 't' :
            theta_fname = optarg;
            break;
         case 'T' :
            theta_info = optarg;
            break;
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

   /* w(z) */
   F = fopen_err("w_eos", "w", err);
   quitOnError(*err, __LINE__, stderr);
   dump_param(model->cosmo, F);
   for (a=0.01; a<1.0; a+=0.01) {
      fprintf(F, "%f %f\n", a, w_de(model->cosmo, a, err));
      quitOnError(*err, __LINE__, stderr);
   }
   fclose(F);

   /* Growth factor */
   F = fopen("D_plus", "w");
   for (a=0.1; a<1.0; a+=0.01) {
      fprintf(F, "%f %.6f\n", a, D_plus(model->cosmo, a, 1, err));
      quitOnError(*err, __LINE__, stderr);
   }
   fclose(F);


   /* Lensing efficiency */
   F = fopen("G", "w");
   for (a=1.0,cumG=0.0; a>0.1; a-=0.01) {
      ww = w(model->cosmo, a, 0, err);
      quitOnError(*err, __LINE__, stderr);
      Ga = G(model, a, 0, err);
      quitOnError(*err, __LINE__, stderr);
      cumG += Ga * ww;
      fprintf(F, "%f %g %g\n", a, Ga * ww, cumG);
   }
   fclose(F);



   /* Density power spectrum */

   printf("Calculating density power spectrum, writing file P_delta\n");

   if (model->cosmo->nonlinear == coyote10) {
      set_H0_Coyote(model->cosmo, err);
      quitOnError(*err, __LINE__, stderr);
      printf("Coyote10: Hubble parameter set to h = %g\n", model->cosmo->h_100);
   }

   F = fopen("P_delta", "w");
   dump_param_lens(model, F, 1, err);
   quitOnError(*err, __LINE__, stderr);

   //for (k=ksim[0]/model->cosmo->h_100; k<=ksim[nsim-1]/model->cosmo->h_100; k*=1.05) {

   fprintf(F, "# k[Mpc/h] P_NL(k,a)\n");
   fprintf(F, "# k[Mpc/h]    ");
   for (a=1.0; a>=0.5; a-=0.5) {
      fprintf(F, "        a=%g", a);
   }
   fprintf(F, "\n");
   for (k=0.0001; k<=100.0; k*=1.05) {
      f = k*k*k/(2.0*pi_sqr);
      //f = 1.0;
      fprintf(F, "%e ", k);
      for (a=1.0; a>=0.5; a-=0.5) {
	 if (a > model->cosmo->a_min) {
	    fprintf(F, "%e ", P_NL(model->cosmo, a, k, err)*f);
	    quitOnError(*err, __LINE__, stderr);
	 } else {
	    fprintf(F, "undef ");
	 }
      }
      fprintf(F, "\n");
      quitOnError(*err, __LINE__, stderr);
   }
   fclose(F);


   /* Redshift distribution */
   printf("Printing redshift distribution to file nz\n");
   F = fopen("nz", "w");
   dump_param_lens(model, F, 1, err);
   quitOnError(*err, __LINE__, stderr);

   for (z=0; z<=10; z+=0.01) {
      fprintf(F, "%.3f", z);
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
	 fprintf(F, " %g", prob(model->redshift, z, i_bin, err));
	 quitOnError(*err, __LINE__, stderr);
      }
      fprintf(F, "\n");
   }
   fclose(F);


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
      f = ell*(ell+1)/twopi;
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
            fprintf(F, " %e", pk*f);
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
               fprintf(F, " %e", pg*f);
            } else {
               pg = 0.0;
            }
         }
      }
      fprintf(F, "\n");
   }
   fclose(F);


   /* Shear correlation functions */

   if (theta_fname != NULL) {

      Nrec = 0;
      theta_list = (double*)read_any_list_count(theta_fname, &Nrec, "%lg", sizeof(double), &Ntheta, err);
      quitOnError(*err, __LINE__, stderr);

   } else {

      if (theta_info != NULL) {
         substr    = strsep(&theta_info, " "); Theta_min = atof(substr);
         substr    = strsep(&theta_info, " "); Theta_max = atof(substr);
         substr    = strsep(&theta_info, " "); Ntheta    = atoi(substr);
         Theta_fac = pow( Theta_max / Theta_min, 1.0 / ((double)Ntheta - 1.0));
      } else {
         Theta_min  = TH_MIN;
         Theta_fac  = TH_FAC;
         Ntheta     = NTHETA;
      }

      theta_list = malloc_err(sizeof(double) * Ntheta, err);
      quitOnError(*err, __LINE__, stderr);
      for (c=0; c<Ntheta; c++) {
         theta_list[c] = Theta_min * pow(Theta_fac, c);
      }

   }


   printf("Calculating shear correlation function, writing files xi_p, xi_m\n");
 
   F = fopen("xi_p", "w");
   dump_param_lens(model, F, 1, err);
   quitOnError(*err, __LINE__, stderr);

   fprintf(F, "# theta[arcmin] xi+^{mn}(theta)\n# th[']  ");
   for (i_bin=0; i_bin<Nzbin; i_bin++) {
      for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
         if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
         if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
         fprintf(F, "           %d%d", i_bin, j_bin);
      }
   }
   fprintf(F, "\n");

   for (c=0; c<Ntheta; c++) {
      theta = theta_list[c] * arcmin;
      fprintf(F, "%9.3f", theta/arcmin);
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            fprintf(F, " %e", xi(model, +1, theta, i_bin, j_bin, err));
            quitOnError(*err, __LINE__, stderr);
         }
      }
      fprintf(F, "\n");
   }
   fclose(F);

   F = fopen("xi_m", "w");
   dump_param_lens(model, F, 1, err);
   quitOnError(*err, __LINE__, stderr);

   fprintf(F, "# theta[arcmin] xi-^{mn}(theta)\n# th[']  ");
   for (i_bin=0; i_bin<Nzbin; i_bin++) {
      for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
         if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
         if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
         fprintf(F, "           %d%d", i_bin, j_bin);
      }
   }
   fprintf(F, "\n");

   for (c=0; c<Ntheta; c++) {
      theta = theta_list[c] * arcmin;
      fprintf(F, "%9.3f", theta/arcmin);
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            fprintf(F, " %e", xi(model, -1, theta, i_bin, j_bin, err));
            quitOnError(*err, __LINE__, stderr);
         }
      }
      fprintf(F, "\n");
   }
   fclose(F);


   /* Tophat shear variance */

   printf("Calculating tophat shear variance, writing file gammasqr\n");
   F = fopen("gammasqr", "w");
   dump_param_lens(model, F, 1, err);
   quitOnError(*err, __LINE__, stderr);

   fprintf(F, "# theta[arcmin] <|gamma|^2>^{mn}(theta)\n#        ");
   for (i_bin=0; i_bin<Nzbin; i_bin++) {
      for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
         if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
         if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
         fprintf(F, "           %d%d", i_bin, j_bin);
      }
   }
   fprintf(F, "\n");

   for (c=0; c<Ntheta; c++) {
      theta = theta_list[c] * arcmin;
      fprintf(F, "%9.3f", theta/arcmin);
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            fprintf(F, " %e", gamma2(model, theta, i_bin, j_bin, err));
            quitOnError(*err, __LINE__, stderr);
         }
      }
      fprintf(F, "\n");
   }
   fclose(F);


   /* Aperture mass variance */

   printf("Calculating aperture mass variance (polynomial filter), writing files mapsqr, mapsqr_gauss\n");
   F = fopen("mapsqr", "w");
   dump_param_lens(model, F, 1, err);
   quitOnError(*err, __LINE__, stderr);

   fprintf(F, "# theta[arcmin] <M_ap^2>^{mn}(theta){polynomial}\n#        ");
   for (i_bin=0; i_bin<Nzbin; i_bin++) {
      for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
	 if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
	 if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
	 fprintf(F, "           %d%d", i_bin, j_bin);
      }
   }
   fprintf(F, "\n");

   for (c=0; c<Ntheta; c++) {
      theta = theta_list[c] * arcmin;
      fprintf(F, "%9.3f", theta/arcmin);
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
	 for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
	    if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
	    if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
	    fprintf(F, " %e", map2_poly(model, theta, i_bin, j_bin, err));
	    quitOnError(*err, __LINE__, stderr);
	 }
      }
      fprintf(F, "\n");
   }
   fclose(F);

   F = fopen("mapsqr_gauss", "w");
   dump_param_lens(model, F, 1, err);
   quitOnError(*err, __LINE__, stderr);

   fprintf(F, "# theta[arcmin] <M_ap^2>^{mn}(theta){Gaussian}\n#        ");
   for (i_bin=0; i_bin<Nzbin; i_bin++) {
      for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
	 if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
	 if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
	 fprintf(F, "           %d%d", i_bin, j_bin);
      }
   }
   fprintf(F, "\n");

   for (c=0; c<Ntheta; c++) {
      theta = theta_list[c] * arcmin;
      fprintf(F, "%9.3f", theta/arcmin);
      for (i_bin=0; i_bin<Nzbin; i_bin++) {
         for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
            if (model->tomo==tomo_auto_only && i_bin!=j_bin) continue;
            if (model->tomo==tomo_cross_only && i_bin==j_bin) continue;
            fprintf(F, " %e", map2_gauss(model, theta, i_bin, j_bin, err));
            quitOnError(*err, __LINE__, stderr);
         }
      }
      fprintf(F, "\n");
   }
   fclose(F);


   printf("Done\n");

   free_parameters_lens(&model);

   return 0;
}
#undef TH_MIN
#undef TH_FAC
#undef NTHETA
#undef ELL_MIN
#undef ELL_MAX
#undef NELL
