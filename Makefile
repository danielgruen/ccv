######## INSTALLATION INSTRUCTIONS ########
#### (0) Make sure that you have all necessary libraries / software installed, including
####     -- the tmv-cpp library, http://tmv-cpp.googlecode.com/
####     -- the CCfits library, http://heasarc.gsfc.nasa.gov/fitsio/CCfits/
####     -- curl for some automated downloading of templates
#### (1) Edit the following lines to match your system / preferences
#### (2) Run 'make' 
####     -- This should compile all code needed to produce / handle the model.
####     -- This should also download a bunch of pre-computed files from my system.
####        If you are asking for the model at a redshift that I haven't pre-computed, it should 
####        automatically do the calculation for you (but it may take several processor-days)
#### (3) You can get the covariance of kappa in your selected set of annuli for any combination
####     of redshift and mass by running ./src/getmodel, which outputs a fits file.
####     See ./src/getmodel.cpp for how to handle this in c++ directly.
####     The analogous covariance matrix for gamma is output by ./src/getmodel_g
####     The analogous covariance matrix for DeltaSigma is output by ./src/getmodel_ds

###########################################
### edit these to have the right c++ compiler and include/library paths for tmv, blas, CCfits

CPP=g++ -std=c++0x -fopenmp 
# this should work on linux / gcc

#CPP=c++ -DNO_OMP -DUSE_MULTIMAP
# this sould work on MacOS / clang

INCLUDES=-I ~/werc3/include -I ~/include 
# add wherever else tmv and CCfits may be

FFTW_PREFIX=/data/soft/fftw/3.3.6/
GSL_PREFIX=`gsl-config --prefix`

LIBFLAGS=-L /sw/lib -L ~/werc3/lib -L ~/lib -lCCfits -lcfitsio
LIBFLAGS_TMV=-ltmv -ltmv_symband -lblas -lpthread
LIBFLAGS_GSL=`gsl-config --libs --cflags`


### this is how many processors you'd like to use; only worry about this for non-pre-computed redshifts
CORES=1

### cluster definition file
### simple format with one line per cluster with ID z_lens p_z_source annuli_prefix
### where ID is some unique string identifier for the cluster, 
###       z_lens its redshift (choose one between 0.15 and 0.60 with two significant digits to use precomputed templates and save a lot of time), 
###       p_z_source the lensing-weighted source redshift distribution prefix in the cluster field (expected in this directory with suffix .tab)
###                  format is as described in http://www2.iap.fr/users/kilbinge/nicaea/readme-nicaea_2.4.txt
###                  example for a single fixed source redshift is in pz.tab
###                  example for a realistic distribution estimated from data is in stack_pz.tab
###       annuli_prefix is the annuli definition file prefix (expected in this directory with suffix .tab)
###                  simple format with N_annuli in the first line and then one line of theta_min theta_max for each annulus
CLUSTERS=default.tab
#CLUSTERS=codex.tab
#CLUSTERS=zlist.tab

######## END INSTALLATION INSTRUCTIONS ########

VERSION=0.3

REDSHIFTS=$(shell cat $(CLUSTERS) | cut -d \  -f 2 | sort | uniq) # a list of redshifts
PZFILES=$(shell cat $(CLUSTERS) | cut -d \  -f 3 | sort | uniq) # a list of p(z) prefixes
ANNULI=$(shell cat $(CLUSTERS) | cut -d \  -f 4 | sort | uniq) # a list of annuli prefixes




all: software lut templates model info

info:
	@echo z=$(REDSHIFTS)
	@echo pz=$(PZFILES)
	@echo annuli=$(ANNULI)

##### all

software: lut_software template_software model_software tinker nicaea lambda_at_m_z	# programs to calculate covariances

lut: lut_2pc lut_W lut_U lut_sigmam lut_dndM lut_rho0 lut_profiles lut/Pkappa.tab	# look-up tables

templates: templates_corrh templates_conc templates_ell templates_off templates_lss     # templates for intrinsic covariance components

model: model_corrh model_conc model_ell model_off model_lss				# co-added and resampled temples according to some binning scheme


#### software

lambda_at_m_z: src/lambda_at_m_z.cpp src/cosmology.h lut_dndM
	$(CPP) -o src/lambda_at_m_z src/lambda_at_m_z.cpp

template_software: src/template_corrh src/template_corrh_combine src/template_conc src/template_ell src/template_off src/template_lss

model_software: src/resample_ell src/resample_ell_g src/resample_off src/resample_off_g src/resample_conc src/resample_conc_g src/resample_conc_g src/resample_corrh src/resample_corrh_g src/getmodel src/getmodel_g src/getmodel_ds

lut_software: src/calc_W 

template_kappa_gamma_test: src/template_conc_g src/resample_conc_g_test templates_conc_g model_conc_g model_conc

tinker: src/mktinkerconf
	$(MAKE) -C src/tinker

src/mktinkerconf: src/mktinkerconf.cpp src/cosmology.h
	$(CPP) -o src/mktinkerconf src/mktinkerconf.cpp	

nicaea: src/mknicaeaconf
	FFTW_PREFIX=$(FFTW_PREFIX) GSL_PREFIX=$(GSL_PREFIX) LIBFLAGS_GSL=$(LIBFLAGS_GSL) $(MAKE) -C src/nicaea_2.5/Demo


src/mknicaeaconf: src/mknicaeaconf.cpp src/cosmology.h
	$(CPP) -o src/mknicaeaconf src/mknicaeaconf.cpp	

#### lut

lut_2pc:
	@echo "don't know how to make lut_2pc, but will ignore that"

lut_W:
	@echo "========== checking for availability of W lut =========="
	bash helpers/lut_W.sh $(REDSHIFTS)
	@echo "========== finished with W lut =========="

lut_U:
	@echo "don't know how to make lut_U, but will ignore that"

lut_sigmam:
	@echo "don't know how to make lut_sigmam, but will ignore that"

lut_dndM:
	@echo "========== checking for availability of dndM lut =========="
	bash helpers/lut_dndM.sh $(REDSHIFTS)
	@echo "========== finished with dndM lut =========="

lut_rho0:
	@echo "don't know how to make lut_rho0, but will ignore that"

lut_profiles:
	@echo "don't know how to make lut_profiles, but will ignore that"

lut/Pkappa.tab:
	@echo "========== checking availability of Pkappa ==========="
	bash helpers/pkappa.sh $(PZFILES)
	@echo "========== finished with Pkappa ==========="

### lut_software

src/calc_W: src/calc_W.cpp src/corrh/corrh.h src/enfw/enfw.h src/profile/profile.h src/cosmology.h
	$(CPP)  -o src/calc_W $(INCLUDES) $(LIBFLAGS) src/calc_W.cpp

#### templates

templates_corrh:
	@echo "========== checking for availability of corrh templates =========="
	bash helpers/templates_corrh.sh $(REDSHIFTS) $(CORES)
	bash helpers/templates_corrh_combine.sh $(REDSHIFTS)
	@echo "========== finished with corrh templates =========="

templates_conc:
	@echo "========== checking for availability of conc templates =========="
	bash helpers/templates_conc.sh $(CORES)
	@echo "========== finished with conc templates =========="

templates_conc_g:
	@echo "========== checking for availability of conc_g test templates =========="
	bash helpers/templates_conc_g.sh $(CORES)
	@echo "========== finished with conc_g test templates =========="


templates_ell:
	@echo "========== checking for availability of ell templates =========="
	bash helpers/templates_ell.sh $(CORES)
	@echo "========== finished with ell templates =========="

templates_off:
	@echo "========== checking for availability of off templates =========="
	bash helpers/templates_off.sh $(CORES)
	@echo "========== finished with ell templates =========="

templates_lss:
	@echo "========== checking for availability of lss templates =========="
	bash helpers/templates_lss.sh $(PZFILES)                         
	@echo "========== finished with lss templates =========="

### template_software

src/template_corrh: src/template_corrh.cpp src/corrh/template_corrh.h src/corrh/corrh.h src/enfw/enfw.h src/profile/profile.h src/cosmology.h
	$(CPP)  src/template_corrh.cpp -o src/template_corrh $(INCLUDES) $(LIBFLAGS) 

src/template_corrh_combine: src/template_corrh_combine.cpp src/corrh/template_corrh.h src/corrh/corrh.h src/enfw/enfw.h src/profile/profile.h src/cosmology.h
	$(CPP)  src/template_corrh_combine.cpp -o src/template_corrh_combine $(INCLUDES) $(LIBFLAGS)

src/template_conc: src/template_conc.cpp src/conc/template_conc.h src/enfw/enfw.h src/cosmology.h
	$(CPP)  src/template_conc.cpp -o src/template_conc $(INCLUDES) $(LIBFLAGS) 
 
src/template_conc_g: src/template_conc_g.cpp src/conc/template_conc.h src/enfw/enfw.h src/cosmology.h
	$(CPP)  src/template_conc_g.cpp -o src/template_conc_g $(INCLUDES) $(LIBFLAGS) 
 
src/template_ell: src/template_ell.cpp src/enfw/template_ell.h src/enfw/enfw.h src/cosmology.h src/filter/filter.o
	$(CPP)  src/template_ell.cpp src/filter/filter.o -o src/template_ell $(INCLUDES) $(LIBFLAGS) 
 
src/template_off: src/template_off.cpp src/off/template_off.h src/enfw/enfw.h src/cosmology.h src/filter/filter.o
	$(CPP)  src/template_off.cpp src/filter/filter.o -o src/template_off $(INCLUDES) $(LIBFLAGS) 

src/template_lss: src/template_lss.cpp src/corrh/template_corrh.h src/cosmology.h
	$(CPP)  src/template_lss.cpp -o src/template_lss $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_GSL) 

#### model

model_corrh: 
	@echo "========== checking for availability of corrh model =========="
	$(eval ARGS=$(shell cat $(CLUSTERS) | cut -d \  -f 2,4 )) # a list of redshifts and annuli
	bash helpers/model_corrh.sh $(ARGS)
	@echo "========== finished with corrh model =========="

model_conc:
	@echo "========== checking for availability of conc model =========="
	$(eval ARGS=$(shell cat $(CLUSTERS) | cut -d \  -f 2,4 )) # a list of redshifts and annuli
	bash helpers/model_conc.sh $(ARGS) $(CORES)
	@echo "========== finished with conc model =========="

model_conc_g:
	@echo "========== checking for availability of conc_g test model =========="
	$(eval ARGS=$(shell cat $(CLUSTERS) | cut -d \  -f 2,4 )) # a list of redshifts and annuli
	bash helpers/model_conc_g.sh $(ARGS) $(CORES)
	@echo "========== finished with conc model =========="


model_ell:
	@echo "========== checking for availability of ell model =========="
	$(eval ARGS=$(shell cat $(CLUSTERS) | cut -d \  -f 2,4 )) # a list of redshifts and annuli
	bash helpers/model_ell.sh $(ARGS) $(CORES)
	@echo "========== finished with ell model =========="

model_off:
	@echo "========== checking for availability of off model =========="
	$(eval ARGS=$(shell cat $(CLUSTERS) | cut -d \  -f 2,4 )) # a list of redshifts and annuli
	bash helpers/model_off.sh $(ARGS) $(CORES)
	@echo "========== finished with off model =========="

model_lss:
	@echo "========== checking for availability of lss model =========="
	$(eval ARGS=$(shell cat $(CLUSTERS) | cut -d \  -f 3,4 )) # a list of p(z) and annuli
	bash helpers/model_lss.sh $(ARGS) $(CORES)
	@echo "========== finished with lss model =========="

### model software

src/resample_ell: src/resample_ell.cpp src/enfw/template_ell.h src/cosmology.h
	$(CPP)  src/resample_ell.cpp -o src/resample_ell $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)
src/resample_ell_g: src/resample_ell_g.cpp src/enfw/template_ell.h src/cosmology.h
	$(CPP)  src/resample_ell_g.cpp -o src/resample_ell_g $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)

src/resample_off: src/resample_off.cpp src/off/template_off.h src/cosmology.h
	$(CPP)  src/resample_off.cpp -o src/resample_off $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)
src/resample_off_g: src/resample_off_g.cpp src/off/template_off.h src/cosmology.h
	$(CPP)  src/resample_off_g.cpp -o src/resample_off_g $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)

src/resample_conc: src/resample_conc.cpp src/conc/template_conc.h src/cosmology.h
	$(CPP)  src/resample_conc.cpp -o src/resample_conc $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)
src/resample_conc_g: src/resample_conc_g.cpp src/conc/template_conc.h src/cosmology.h
	$(CPP)  src/resample_conc_g.cpp -o src/resample_conc_g $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)

src/resample_conc_g_test: src/resample_conc_g_test.cpp src/conc/template_conc.h src/cosmology.h
	$(CPP)  src/resample_conc_g_test.cpp -o src/resample_conc_g_test $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)

src/resample_corrh: src/resample_corrh.cpp src/corrh/template_corrh.h src/cosmology.h
	$(CPP)  src/resample_corrh.cpp -o src/resample_corrh $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)
src/resample_corrh_g: src/resample_corrh.cpp src/corrh/template_corrh.h src/cosmology.h
	$(CPP)  src/resample_corrh_g.cpp -o src/resample_corrh_g $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)

src/getmodel: src/getmodel.cpp src/nicaea_pz.h src/model/covariance.h src/cosmology.h src/enfw/enfw.h
	$(CPP) src/getmodel.cpp -o src/getmodel $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)
src/getmodel_g: src/getmodel_g.cpp src/nicaea_pz.h src/model/covariance.h src/cosmology.h src/enfw/enfw.h
	$(CPP)  src/getmodel_g.cpp -o src/getmodel_g $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)
src/getmodel_ds: src/getmodel_ds.cpp src/nicaea_pz.h src/model/covariance.h src/cosmology.h src/enfw/enfw.h
	$(CPP)  src/getmodel_ds.cpp -o src/getmodel_ds $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)

### filter library

src/filter/filter.o: src/filter/filter.cpp src/filter/filter.h
	$(CPP)  -o src/filter/filter.o $(INCLUDES) $(LIBFLAGS) -c src/filter/filter.cpp

#### clean: re-compile software afterwards

clean:
	rm -f src/filter/filter.o src/getmodel_g src/getmodel src/resample_corrh_g src/resample_corrh src/resample_conc_g src/resample_conc src/resample_ell_g src/resample_ell src/template_ell src/template_conc src/template_corrh src/template_corrh_combine
	$(MAKE) clean -C src/tinker
	$(MAKE) clean -C src/nicaea_2.5/Demo

#### forget about model (and re-do later, e.g. if you have changed your annuli definition)

clean_model:
	rm -rf model/*

#### pack templates and make available online (to be run by Daniel...)

pub: pub/templates_conc.tar.gz pub/templates_ell.tar.gz pub/ccv.tar.gz
	rsync -avv pub/templates_conc.tar.gz pub/templates_ell.tar.gz templates/corrh_*fits lut/W*tab lut/dndM*tab dgruen@moon.usm.uni-muenchen.de:/usr/web/users/dgruen/public_html/code/templates/
	rsync -avv pub/ccv-$(VERSION).tar.gz dgruen@moon.usm.uni-muenchen.de:/usr/web/users/dgruen/public_html/code/

pub/templates_conc.tar.gz: templates/conc/cov_100.fits
	tar czf pub/templates_conc.tar.gz templates/conc/*.fits

pub/templates_ell.tar.gz: templates/ell/cov_100.fits
	tar czf pub/templates_ell.tar.gz templates/ell/*.fits

pub/ccv.tar.gz:
	-git commit -a
	bash helpers/pack_software.sh $(VERSION)
