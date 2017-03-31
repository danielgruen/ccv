This is an implementation of the model for the intrinsic covariance of projected density profiles of galaxy clusters described in Gruen et al., 2015, MNRAS, 449, 4, 4264, including contributions from (1) uncorrelated large scale structure, (2) subhaloes, (3) halo concentration variations, (4) halo ellipticity and orientation variations and (5) redMaPPer-like offcentring. It provides the intrinsic covariance of kappa, gamma or DeltaSigma for a given set of annuli, a given mass and redshift.


KNOWN LIMITATIONS
-----------------
* only the kappa version of covariance at the snapshot redshifts z=0.24533 and z=0.5 has been validated against simulations
* off-centering covariance cuts some corners in applying the Rykoff+2016 model


WHERE TO BEGIN
--------------
* see Makefile for installation instructions; then type make to build and download pre-computed templates (if available)
* send me <dgruen@stanford.edu> a note or post on https://github.com/danielgruen/ccv in case of any problems
* see src/getmodel* for simple c++ codes that build a model covariance for kappa, gamma and DeltaSigma
