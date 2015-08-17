============================================================================
nicaea
Version 2.5 (10/2014)
NumerIcal Cosmology And lEnsing cAlculations

Web page: http://cosmostat.org/nicaea
(Old web page: http://www2.iap.fr/users/kilbinge/nicaea)

Authors:
Martin Kilbinger
Karim Benabed
Jean Coupon (HOD, halomodel)
Henry J. McCracken (HOD)

Contributors:
Liping Fu (decomp_eb)
Catherine Heymans (intrinsic alignment)
============================================================================


Notes:

============================================================================
0. Documentation
============================================================================

nicaea is the cosmology part of cosmo_pmc, which can be downloaded for free at
www.cosmopmc.info. From that site, the cosmo_pmc manual is available for
further information about nicaea, which are more detailed than covered in this
readme.


============================================================================
1. Download, compile, and run nicaea
============================================================================

---------------------
1.1 Download the code
---------------------

Download the file nicaea_2.5.tgz from http://cosmostat.org/nicaea and un-tar
the archive. The packages fftw3 and gsl are required to compile and run nicaea.
You can install fftw3 from http://www.fftw.org, and gsl from
www.gnu.org/software/gsl.

--------------------------------
1.2 Compile and install the code
--------------------------------

Two options to compile nicaea exist. If nicaea is to be used as a library,
option 1 is recommanded.

Option 1) Using cmake.

	cd build
	cmake ..
	make && make install

	The last command will copy the executable demo programs (e.g. lensingdemo)
	to <BASE>/bin, the library libnicaea.a to <BASE>/lib, and the include
	files to <BASE>/include/nicaea . The default base directory is
	<BASE>=nicaea_2.5 .

	The code can be tested with 'ctest -vv'
	To run the demo programs (see below), go to nicaea_2.5/par_files .

Option 2) Using make.

	cd Demo
	make

	If fftw3 and gsl are not installed in a standard directory (e.g. /usr,
	/usr/local), set the variables 'FFTW' and 'GSL' in the Makefile. The
	header file fftw3.h is looked for in $(FFTW)/include and libfftw3.a in
	$(FFTW)/lib. The gsl header files are looked for in $(GSL)/include, the
	libraries libgsl.a and libgslcblas.a in $(GSL)/lib.

	Various demo programs can be run in ./Demo, see below.

-------------------------
1.3 Run the demo programs
-------------------------

The demo programs need parameter files in the working directory, which can be
found in par_files.

lensingdemo		Weak lensing: density- and lensing power spectrum,
			lensing second-order functions
sn1demo			SNIa: Luminosity distance, distance module
halomodeldemo		Halo model: Power spectrum
cmb_bao_demo 		CMB and BAO: geometrical quantities, e.g. sound
			horizon, angular diameter distance
decomp_eb_demo		Weak lensing: E-/B-mode decomposition (generaalized
			ring statistic)
cosebi_demo		Weak lensing: E-/B-mode decomposition (COSEBIs)
third_order_demo	Weak lensing: Third-order aperture-mass moments


============================================================================
2. Main functions
============================================================================


Second-order shear statistics
-----------------------------

xi(model, pm, theta, i_bin, j_bin, err) 	     Two-point correlation function
	      	     		   		     xi+ (pm=0) and xi- (pm=1) at
	      	     		     		     angular scale theta [rad]

gamma2(model, theta, i_bin, j_bin, err)		     Top-hat shear variance in a circle
	      	     		     		     of radius theta [rad]

map2_poly(model, theta, i_bin, j_bin, err)	     Aperture-mass variance, polynomial filter
map2_gauss(model, theta, i_bin, j_bin, err)	     Aperture-mass variance, Gaussian filter

RR(model, THETA_MIN, THETA_MAX, a, N, poly, pm ,err) 'Ring' statistics, with Chebyshev-filter
	  	     		      	       	     function decomposition, see Fu &Kilbinger
						     (2010)

E_cosebi(model, n, Psimin, Psimax, i_bin, j_bin, path, B_cosebi, err);
						     COSEBIs (Complete Orthogonal E-/B-mode
						     Integrals), Schneider, Eifler & Krause (2010)

map3(model, R, i_bin, i_bin, k_bin, wfilter, err)    Third-order aperture-mass generalized moment
						     Schneider, Kilbinger & Lombardi (2005)


The value of the corresponding two-point function is returned as
double. i_bin and j_bin are the redshift-bin indices.
For the structure cosmo_lens *model and error **err, see Sect. 4.


Power spectra
-------------

P_NL(model, a, k, err)		     3d power spectrum of delta
Pshear(model, s, i_bin, j_bin, err)  2d shear power spectrum: Pkappa or Pkappa+Pg^(1) if reduced-
	      	 	       	      correction is switched on with key "sreduced = K10" in
				      cosmo_lens.par parameter file. In earlier versions, this
				      function was called 'Pkappa'.
Pg1(model, s, i_bin, j_bin, err)     2d reduced-shear correction power spectrum Pg^(1), see Kilbinger
	      	     	    	     (2010). The totel (reduced-shear) power spectrum is
				     Pkappa + Pg1.
arguments:
model				     structure containing the cosmological parameters
theta, THETA_MIN, THETA_MAX          Angular separations/smoothing scales, in rad
a			             Scale factor, max(0.01,1/(1+zmax))<=a<1.0
k	                             3d Fourier wave-mode in h/Mpc
s	                             2d Fourier wave-mode, 1e-2<=ell<=1e6
i_bin, j_bin			     Redshift-bin indices


Ranges
------

The range for k is unlimited except for the coyote10 and coyote13 non-linear emulators.
For k<3.3e-6 h/Mpc and k>333 h/Mpc, the
power spectrum is extrapolated (see below). The limits can be changed
in cosmo.h.

The reduced-shear correction fits are accurate to 2% beetween ell=0.1 and 2*10^5. Outside
that range, Pg^(1) return zero.

The range for theta is very, very large, it is determined
in the routine xi_via_hankel. Although the Hankel transform is
accurate only on a much smaller interval, the range of acceptable
results is still from sub-arcseconds to a couple of degrees.

The limited range of the reduced-shear correction reflects in a smaller valid angular range
of xi+ and xi-. If the reduced-shear is switched on, the ranges within which the second-order
functions are affected to small fractions of a percent are:
xi+          [0.1';1000']
xi-          [0.5';1000']
mapsqr       [0.2';1000']
gammasqr     [0.1':1000']
mapsqr_gauss [0.1';1000']


============================================================================
3. Cosmology
============================================================================

The cosmology is encoded in the structure cosmo. It contains all
relevant cosmological and nuisance parameters, and pre-calculated
tables and constants. If parameters change, these tables are
recomputed once they are needed. All lensing-related variables are
contained in the structure cosmo_lens.

-----------------------------------
3.1. Reading parameters from a file
-----------------------------------

The function

    read_cosmological_parameters_lens(&model, F, err)

reads cosmological and lensing parameters from the file F (type FILE*) and
initialised the structure cosmo_lens *model. The file 'cosmo_lens.par' is an
example file. First, it contains a reference to the basic cosmology file 'cosmo.par',
containing cosmological parameters. Next, redshift information is read from
the file 'nofz.par'. Then, the lensing parameters follow.

-------------------------------
3.2. Initializing the cosmology
-------------------------------

The function

    init_parameters_lens(...)

returns a pointer to the structure cosmo_lens with parameters given by
the arguments and blank tables. If passed to a function (e.g. one
described in Sect.2), the corresponding tables and constants (if
required) are filled and calculated. Successive calls to this function
will be very fast since only a linear interpolation of the tabulated
values is performed.

--------------------------
3.3 Changing the cosmology
--------------------------

If a different cosmology is required, a new cosmo_lens pointer has to be
created, either with

    model_new = init_parameters_lens(...)

as above, or with

    model_new = copy_parameters_lens_only(model, err).
    model_new->param1 = ...
    model_new->param2 = ...
    ...

In both cases, all tables and constants are blanked. A call of

       updateFrom_lens(model_new, model, err)

copies tables from model to model_new if corresponding parameters are
unchanged and leaves those blank which have to be recalculated if
required. This is particularly efficient if only a few or only "fast"
parameters change since a small number of (time-consuming) functions
will be recalculated. E.g., if only the redshift parameters change,
the non-linear power spectra and growth factor need not be
recalculated, only the shear statistics, which is very fast due to the
Hankel transform.

-------------------------
3.4 Parameters and ranges
-------------------------

The following parameters are implemented. Within a given range, the
program should obtain reasonable results or return an error message (see
Sect.4). The program does not check whether a parameter is within its
range. The following ranges have been tested some time ago, probably the code
will work outside of these ranges as well.

Cosmology
---------

Omega_m		total matter density (baryonic + dark)	        [0.1, 1.5]
Omega_de	dark energy density				[0.1, 1.5]
w0_de								[-2.0, -0.5]
w1_de		dark energy eos parametrization (see below)	[-0.6, 0.6]
h_100           Hubble parameter H_0 = 100 h_100 km/s/Mpc	[0.4, 1.0]
Omega_b         baryon density	       	   	 		[0.02, 0.06]
Omega_nu_mass   massive neutrino density			(not tested)
N_eff_mass      Number of massive neutrinos			(not tested)
sigma_8 	power spectrum normalisation sigma_8		[0.1, 1.5]
normalization   = sigma_8 (CMB normalization not supported at the moment)
n_spec		primordial spectral index			[0.7, 1.3]

Redshift parameters
-------------------

The number of redshift bins is Nzbin. For each bin n_bin, the number of
redshift parameters is given by Nnz[n_bin], its type by nofz[n_bin].
The sub-array par_nz[n_bin*Nn_max .. n_bin*Nnz_max+Nnz[n_bin]] contains the
Nnz[n_bin] redshift parameters of bin n_bin. For all types the first two
parameters define the minimum and maximum redshift: par_nz[n_bin*Nn_max]
= zmin par_nz[n_bin*Nn_max+1] = zmax.

The following types exist:

nofz    Nnz    parameters              prob(z) (for zmin<z<zmax)
--------------------------------------------------------------------------------
ludo    5      alpha_p, beta_p, z0     (z/z0)^alpha_p * exp(-(z/z0)^(beta_p))
jonben  5      a, b, c                 z^a/(z^b + c)
ymmk    5      a, b, c                 (z^a + z^(ab))/(z^b + c)
single  2      z0, z0		       delta_Dirac(z - z0)
hist    2n+1			       Histogram with n bins

type=hist assumes a N(z) histogram with n bins. The parameters are
stored in the vector par_nz as follows:

0  1  | 2  3  ... n       | n+1 n+2 ... 2n
-----------------------------------------------
z0 zn | z1 z2 ... z_{n-1} | N0  N1  ... N_{n-1}

The number of parameters is Nnz=2n+1. The redshifts z_i are understood
as the lower bin boundaries with the exception of zn=zmax which is the
limiting redshift. The i-th bin therefore is between z_i and
z_{i+1}, the (unnormalized) number of galaxies is N_i.
zmin=z0 and zmax=zn are in the first two entries, as required.

The function read_par_nz_hist reads the histogram data from a file, sets
Nnz and returns par_nz. The file has to have the following structure:

# hist
  z0  	  N0
  z1	  N1
  ...     ...
  z_{n-1} N_{n-1}
  zn      0.0


All galaxies are at a single redshift z0 can achieved with the following file:
# single
  z0
  z0

(The value z0 has to appear twice. It is both zmin and zmax.)


The normalization for all types, \int_zmin^zmax prob(z) dz = 1, is calculated in
the code.

Flags
-----

nonlinear	linear:            Linear power spectrum (BBKS CDM transfer function)
		pd96:              Peacock & Dodds (1996) fitting formula
		smith03:           Smith et al. (2003) halofit
		smith03_de:	   Smith et al. (2003) halofit + dark-energy correction from icosmo.org
		smith03_revised    Takahashi et al. (2012), revised halofit parameters
		coyote10:	   Coyote emulator v1, Heitmann, Lawrence et al. 2009, 2010
		coyote13:	   Coyote emulator v2, Heitmann et al. 2013

transfer	bbks:              Bardeen et al. (1986) transfer function
		eisenhu:           Eisenstein & Hu (1998) "shape fit"
		camb:              Using camb for T(k) (not yet supported)

growth		heath:             Heath (1977) analytical expression for
			           linear growth factor (valid only for no/a
			           pure cosmological constant, i.e. w0_de=-1,
			           w1_de=0)
		growth_de:         General dark energy model
		camb_vinschter_gr: Growth in camb T(k,z)

de_param	jassal:            w(a) = w0_de + w1_de*a*(1-a)
		linder:            w(a) = w0_de + w1_de*(1-a)
		earlyDE:           w(a) = w0_de/sqrt(1-b_early*log(a))

normmode 	norm_s8:           normalization = sigma_8
		From outside, the pair (normalization, normmode) should be
		used to set the normalization. The function set_norm(mode,
		self) automatically sets sigma_8 correspondingly
		(called from init_parameters and updateFrom). Internally,
		sigma_8 is used throughout the program.

tomo		tomo_all           All redshift-correlations (ij), i<=j
		tomo_auto_only	   Only auto-correlations (ii)
		tomo_cross_only	   Only cross-correlations (i!=j)

reduced		none		   No reduced-shear correction
		K10		   Reduced-shear according to Kilbinger (2010)

q_mag_size	double		   If reduced==K10: q_mag_size = 2*(alpha+beta-1), see
				   K10 eq. 16. Set q_mag_size = 0 if no magnification/size
				   bias correction to be added (reduced-shear only).

sia		none		   No intrinsic alignment (IA)
		HS04		   Hirata & Seljak linear IA model

sia_terms	GI_II		   If sia!=none: IA terms to be added, 'GI_II' (both GI and II),
		only_GI		   'only_GI' (only GI), or 'only_II' (only II).
		only_II

A_ia		double		   If sia!=none: Global amplitude of IA contribution.

The range for w0_de and w1_de correspond to de_param=linder.

The minimum scale factor a_min (used for various integrations) is set using
the function set_amin().


============================================================================
5. Errors and diagnostics
============================================================================

Most of the situations where an error or undefined value occurs are
intercepted by the program. In that case, a variable *err of type error* is
set via the macros

      *err = addError(error_type, "message", *err, __LINE__)
or
      *err = addErrorVA(error_type, "formatted message", *err, __LINE__, VA_LIST)

storing the line in the code, a message and the error type
(ce_xyz). With

      testErrorRet(test, error_type, "message", *err, __LINE__, return_value)
or
      testErrorRetVA(test, error_type, "formatted message", *err, __LINE__, return_value, VA_LIST)


a conditional error is produced if the (Boolean) expression test is
true. The error can be transported up the stack to the calling
function with the macro

      forwardError(*err, __LINE__, return_value)

(omit return_value in case of a void function). This can be used as
diagnostics even for errors deep in the hierarchy of functions. To exit on
an error, use

      exitOnError(*err, FILE).

At the start of the program, or after an error had occurred but one wishes
to continue, maybe with a different cosmology, set *err = NULL.

An error can be caused by undefined values, not initialized parameters,
function arguments outside the valid range. Further, a specific cosmology may
not allow certain functions to be carried out. For example, in a loitering
Universe there is a maximum redshift, and if the redshift distribution extends
over this maximum, the angular diameter distance is undefined and an error is
produced.


============================================================================
6. Extrapolation
============================================================================

In the highly non-linear regime, the power spectrum is
extrapolated. For the linear power spectrum, P(k) \propto
k^{n_spec-4.0} is assumed. In the PD96-case, the stable clustering
result P(k) \propto k^{-2.5} is used. For Smith, the asymptotic form
of the halofit formula is taken, see Rob's paper eq.(61).

In the linear regime at small k, the extrapolation is P(k) \propto k^n_spec.


============================================================================
7. Performance
============================================================================

Time-consuming functions store tabulated values and interpolated when
called after the first time. The tables are recalculated when
cosmological parameters have changed since the previous call. The
correlation functions are calculated using a fast Hankel transform.


============================================================================
8. Known bugs and shortcomings
============================================================================

- Some parameter combinations cause undefined behaviour of the
  program. These are (hopefully) intercepted and an error is created
  (see Sect. 5). E.g., for n_spec<0.7, f_NL (Peacock&Dodds) is not
  defined. For a closed Universe, the probed redshift can be larger
  than the maximum redshift.

- a=1.0 very rarely creates an error, use 0.99999... instead.

- The code is not well suited for Fisher matrix calculations. In particular
  for the inverse Fisher matrix, numerical derivatives have to be very
  accurate, and the interpolations between tabulated values (linear and
  spline) in nicaea introduce numerical noise that can render the Fisher
  matrix numerically singular (Wolz et al. 2012).

- Dark-energy models, in particular with varying w(z), are not recommended
  for the non_linear models smith03, and smith03_de. Instead, use the
  revised halofit model with smith03_revised.

In case of problems please don't hesitate to contact me at
martin.kilbinger@cea.fr . Questions and comments are welcome!


============================================================================
9. Changes compared to the Rob Smith's original halofit
============================================================================

Parts of the program 'cosmo.c' is based on Rob Smiths' halofit (Smith et al.
2003). The code for determining the non-linear power spectrum has been improved
and made more efficient. The main changes are listed below. The code also
includes the non-linear fitting formulae of Peacock & Dodds (1996).

- Tabulation of the linear and non-linear power spectrum, constants
  are calculated only once.
- Integration cutoff for determination of non-linear scale knl
  flexible, as function of smoothing scale rmid; using Romberg
  integration.
- Bisection to find knl is iterative: if the bisection gets stuck at one
  end of the bisecting interval, the interval is shifted accordingly and
  a new bisection is started. If knl is larger than knlstern (I chose
  10^6 h/Mpc), the bisection is canceled and the linear power spectrum
  is used.
- Slope and curvature are calculated only once, after knl is fixed.
- The Eisenstein&Hu (1998) fit for the transfer function is used
  instead of Bond&Efstathiou (1984).
- The exact linear growth factor is used instead of the CPT92 fitting
  formula. Dark energy models are incorporated.


============================================================================
10. Acknowledgements
============================================================================

We thank Alexandre Boucaud, Jan Hartlap, Alina Kiessling, Jasmin Pielorz, Peter
Schneider, Rob E. Smith, Patrick Simon, and Masahiro Takada for helpful
suggestions.


============================================================================
11. Contact
============================================================================

Feel free to email me at martin.kilbinger@cea.fr

Have fun!
   Martin Kilbinger

