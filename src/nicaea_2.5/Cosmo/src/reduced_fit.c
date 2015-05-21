#include "reduced_fit.h"


/* Parameter limits for (Unused, Omegam, Omegade, w0, Omegab, h100, sigma8, ns) */
const double limits_lower[M_PAR] = {0.0, -0.05, -0.4, -0.6, -0.04, -0.1, -0.15, -0.1};
const double limits_upper[M_PAR] = {0.0,  0.08 , 0.3,  0.4,  0.04,  0.4,  0.13,  0.2};


/* K10 eq. (22) */
double sum_B_a(int alpha, int i, double a)
{
   int j;
   double b_i, pow_a;

   for (j=0,pow_a=1.0,b_i=0.0; j<N_B; j++,pow_a*=a) {
      b_i += B_fit[alpha][i][j]*pow_a;
   }

   return b_i;
}

/* K10 eq. (22) */
double sum_C_a(int alpha, int i, double a)
{
   int j;
   double c_i, pow_a;

   for (j=0,pow_a=1.0,c_i=0.0; j<N_C; j++,pow_a*=a) {
      c_i += C_fit[alpha][i][j]*pow_a;
   }

   return c_i;
}

/* K10 eq. (18) */
double h_piece(double logell, const double b[])
{
   double res, dx, r0, r1, r2, r3, frac;

   if (logell<Y_LOW) {
      res = b[0]*logell + b[1];
   } else if (logell>Y_UP) {
      res = b[2]*logell + b[3];
   } else {
      /* Cubic spline */
      dx = Y_UP - Y_LOW;
      r0 = b[0]*Y_LOW + b[1];
      r1 = dx*b[0];
      r2 = 3.0*(b[2]*Y_UP + b[3]) - dx*b[2] - 3.0*r0 - 2.0*r1;
      r3 = dx*b[2] - 2.0*(b[2]*Y_UP + b[3]) + 2.0*r0 + r1;

      frac = (logell-Y_LOW)/dx;
      res  = ((r3*frac + r2)*frac + r1)*frac + r0;
   }

   return res;
}

/* K10 eq. (21) */
double Q(int alpha, double logell, double a)
{
   double logdQdp, c_i, b[4], res, pow_logell;
   int i;

   for (i=0,logdQdp=0.0,pow_logell=1.0,res=1.0; i<N_POLY; i++,pow_logell*=logell) {
      c_i = sum_C_a(alpha, i, a);
      logdQdp += c_i*pow_logell;
   }
   res *= logdQdp;

   for (i=0; i<N_PL; i++) {
      b[i] = sum_B_a(alpha, i, a);
   }
   res *= exp(h_piece(logell, b));

   return res;
}

/* ============================================================ *
 * Performs the sum over scale factor (K10, eq. 12),  with      *
 * the Taylor-expansion (13,14).				*
 * ============================================================ */
double sum_a_for_Pg1(double logell, double a_min, int Na, double da, const double *fmn, const double **dfmn_dp,
		     const double dpar[M_PAR])
{
   int k, alpha;
   double a, Pg1, Q0, dPg1_dp;

   for (k=0,a=a_min,Pg1=0.0; k<Na; k++,a+=da) {
      Q0   = Q(0, logell, a);
      Pg1 += fmn[k]*Q0;

      for (alpha=1; alpha<M_PAR; alpha++) {
	 /* K10 (14) */
	 dPg1_dp = fmn[k]*Q(alpha, logell, a) + dfmn_dp[alpha][k]*Q0;
	 Pg1    += dpar[alpha]*dPg1_dp;
      }
   }

   Pg1 *= da;

   return Pg1;
}

/* ============================================================ *
 * Returns 0 if difference parameter dpar is within limits      *
 * around 0, and the parameter index (1..M_PAR) otherwise.      *
 * ============================================================ */
int check_limits(const double dpar[M_PAR])
{
   int alpha;

   for (alpha=1; alpha<M_PAR; alpha++) {
      if (dpar[alpha]<limits_lower[alpha]) return alpha;
      if (dpar[alpha]>limits_upper[alpha]) return alpha;
   }

   return 0;
}

/* The fitting matrices, K10 Appendix A */

const double B_fit[M_PAR][N_PL][N_B] = 
  {
    {
      { 1.2157, -1.7061,  3.613, -2.588},
      { 22.054,  14.525, -10.134, -1.2063},
      {-1.7221, -4.2894,  6.9843, -3.7669},
      { 41.719,  34.891, -40.957,  6.3019},
    },
    {
      { 1.2963, -2.2435,  4.7574, -3.4183},
      { 24.802,  15.438, -12.793,  0.30538},
      {-0.42199, -14.769,  28.968, -16.938},
      { 62.035, -22.702, -10.113,  15.12},
    },
    {
      { 1.2314, -1.826,  3.9322, -2.9291},
      { 27.286, -8.0841,  21.999, -18.426},
      {-3.2808,  4.3912, -6.464,  2.8446},
      { 62.447, -73.471,  122.62, -71.766},
    },
    {
      { 1.1772, -1.523,  3.3478, -2.6269},
      { 17.16,  39.066, -56.363,  24.476},
      { 1.7985, -17.673,  25.431, -12.429},
      {-1.2271,  192.08, -259.8 , 113.95},
    },
    {
      { 1.1983, -1.5973,  3.3941, -2.4472},
      { 21.942,  28.174, -32.269,  10.746},
      {-4.4314,  9.5715, -14.092,  6.4373},
      { 75.097, -124.48,  200.67, -109.85},
    },
    {
      { 1.3643, -2.6937,  5.6748, -4.0482},
      { 22.681,  17.078, -13.713,  0.12838},
      {-1.6597, -3.5025,  5.6229, -3.2007},
      { 44.397,  11.923, -3.5699, -10.617},
    },
    {
      { 1.2149, -1.7012,  3.6033, -2.5827},
      { 24.779,  8.4417, -0.27702, -6.3675},
      {-2.5782, -0.084018,  0.53439, -0.56922},
      { 53.073, -11.638,  28.888, -27.294},
    },
    {
      { 1.3469, -3.3556,  6.6862, -4.5164},
      { 25.294,  6.7075,  4.1411, -10.121},
      {-1.0794, -6.1604,  9.2899, -4.7514},
      { 37.736,  49.648, -61.197,  16.44},
    }
  };

const double C_fit[M_PAR][N_POLY][N_C] = 
  {
    {
      { 0.56428,  2.3001, -3.9649,  2.427},
      { 0.12548, -0.94677,  2.59, -2.04},
      { 0.1557, -0.65321,  0.89795, -0.58724},
      {-0.063989,  0.43324, -1.0211,  0.7215},
      { 0.0087569, -0.093102,  0.24571, -0.1693},
      {-0.00055709,  0.008454, -0.022712,  0.015138},
      { 1.4809e-05, -0.00027429,  0.00072501, -0.00046762},
    },
    {
      {-0.43501, -3.1459,  5.6111, -3.4686},
      {-0.087793,  0.80516, -2.761,  2.4422},
      {-0.10227,  0.33267, -0.18646,  0.1432},
      { 0.073105, -0.5126,  1.2891, -0.95999},
      {-0.023766,  0.19417, -0.48231,  0.33371},
      { 0.0030621, -0.024925,  0.059059, -0.039071},
      {-0.0001238,  0.00098714, -0.0022652,  0.0014607},
    },
    {
      {-0.49063, -2.7133,  4.5927, -2.7656},
      { 0.13912, -0.76984,  0.7525,  0.11512},
      {-0.22771,  1.0278, -1.3875,  0.79182},
      {-0.00013136, -0.029954,  0.26883, -0.30505},
      { 0.01479, -0.046285, -0.0031371,  0.04091},
      {-0.0020694,  0.0066385, -0.0028282, -0.0018872},
      { 8.0142e-05, -0.00026039,  0.00016485,  1.1648e-05},
    },
    {
      { 0.25229,  4.0227, -6.6135,  3.7371},
      {-0.52388,  3.3599, -5.5903,  2.5675},
      { 0.35958, -1.7899,  2.4962, -1.27},
      { 0.076671, -0.50746,  0.7332, -0.2437},
      {-0.045615,  0.2431, -0.33587,  0.13382},
      { 0.0052982, -0.026584,  0.03616, -0.014855},
      {-0.00018401,  0.0008918, -0.0012,  0.00049882},
    },
    {
      { 0.12913,  4.7128, -8.0056,  4.5759},
      { 0.13963, -0.11394,  0.61465, -0.97735},
      { 0.26098, -1.3324,  1.9954, -1.1216},
      {-0.11282,  0.48918, -0.99176,  0.70867},
      { 0.016542, -0.075864,  0.19759, -0.1533},
      {-0.0010708,  0.0055907, -0.017359,  0.0139},
      { 2.6106e-05, -0.00015916,  0.00055278, -0.00044773},
    },
    {
      {-0.31944, -3.8608,  6.7541, -4.0259},
      {-0.010529, -0.63367,  1.0269, -0.067878},
      {-0.24657,  1.3721, -2.1149,  1.1977},
      { 0.051632, -0.14509,  0.34221, -0.34002},
      { 0.0010064, -0.030306,  0.0014399,  0.042272},
      {-0.00070106,  0.0052506, -0.0026698, -0.0029094},
      { 3.3959e-05, -0.00020532,  0.00010945,  8.9968e-05},
    },
    {
      { 0.58241,  2.2178, -3.8312,  2.3591},
      { 0.15574, -1.1301,  3.0374, -2.3593},
      { 0.1333, -0.53689,  0.68283, -0.46414},
      {-0.073256,  0.48689, -1.1436,  0.80567},
      { 0.012694, -0.11623,  0.29589, -0.20194},
      {-0.00098447,  0.011061, -0.028339,  0.018735},
      { 2.9427e-05, -0.00036705,  0.00092519, -0.0005942},
    },
    {
      {-0.45031, -3.1685,  5.4862, -3.2435},
      { 0.088613, -0.95773,  1.3824, -0.2107},
      {-0.20975,  1.2013, -1.8198,  1.0069},
      { 0.022754, -0.030319,  0.17275, -0.24886},
      { 0.0067345, -0.052022,  0.030858,  0.028623},
      {-0.0011801,  0.0070072, -0.004861, -0.0020472},
      { 4.8691e-05, -0.00025807,  0.0001707,  6.9607e-05},
    }
  };
