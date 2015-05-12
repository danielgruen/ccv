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

### edit these to have the right c++ compiler and include/library paths for tmv, blas, CCfits
CPP=g++
INCLUDES=-I ~/werc3/include 
LIBFLAGS=-L ~/werc3/lib -L ~/lib -lCCfits -lcfitsio
LIBFLAGS_TMV=-ltmv -ltmv_symband -lblas -lpthread

### this is how many processors you'd like to use; only worry about this for non-pre-computed redshifts
CORES=1

### this is the list of redshifts for which the model should be prepared
# PRE-COMPUTED REDSHIFTS: these are prepared already, templates will be downloaded so you can use them quickly
REDSHIFTS=0.24533 
#0.35 0.187 0.206 0.224 0.234 0.288 0.313 0.348 0.352 0.363 0.391 0.399 0.440 0.450 0.451 0.686 0.15 0.16 0.17 0.18 0.19 0.2 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.3 0.31 0.32 0.33 

# REDSHIFTS TO DO: these are computations in progress, will be put online as soon as they're finished
# 0.34 0.35 0.36 0.37 0.38 0.39 0.4 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.5 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.6

# you can always add your own redshifts to the list and the templates will be calculated (but that may take several processor-days)

### annuli definition file according to what you sent me
### simple format with N_annuli in the first line and then one line of theta_min theta_max each
#ANNULI=annuli_keiichi.tab  # Keiichi Umetsu's CLASH annuli
ANNULI=annuli_daniel.tab  # a set of annuli I've used in the paper

######## END INSTALLATION INSTRUCTIONS ########

VERSION=0.1

all: software lut templates model

##### all

software: lut_software template_software model_software tinker			# programs to calculate covariances

lut: lut_2pc lut_W lut_U lut_sigmam lut_dndM lut_rho0 lut_profiles		# look-up tables

templates: templates_corrh templates_conc templates_ell				# templates for intrinsic covariance components

model: model_corrh model_conc model_ell						# co-added and resampled temples according to some binning scheme


#### software

template_software: src/template_corrh src/template_corrh_combine src/template_conc src/template_ell

model_software: src/resample_ell src/resample_ell_g src/resample_conc src/resample_conc_g src/resample_conc_g src/resample_corrh src/resample_corrh_g src/getmodel src/getmodel_g

lut_software: src/calc_W 

tinker: src/mktinkerconf
	$(MAKE) -C src/tinker

src/mktinkerconf: src/mktinkerconf.cpp src/cosmology.h
	$(CPP) -o src/mktinkerconf src/mktinkerconf.cpp	

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


### lut_software

src/calc_W: src/calc_W.cpp src/corrh/corrh.h src/enfw/enfw.h src/profile/profile.h src/cosmology.h
	$(CPP) -fopenmp -o src/calc_W $(INCLUDES) $(LIBFLAGS) src/calc_W.cpp

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

templates_ell:
	@echo "========== checking for availability of ell templates =========="
	bash helpers/templates_ell.sh $(CORES)
	@echo "========== finished with ell templates =========="

### template_software

src/template_corrh: src/template_corrh.cpp src/corrh/template_corrh.h src/corrh/corrh.h src/enfw/enfw.h src/profile/profile.h src/cosmology.h
	$(CPP) -fopenmp src/template_corrh.cpp -o src/template_corrh $(INCLUDES) $(LIBFLAGS) 

src/template_corrh_combine: src/template_corrh_combine.cpp src/corrh/template_corrh.h src/corrh/corrh.h src/enfw/enfw.h src/profile/profile.h src/cosmology.h
	$(CPP) -fopenmp src/template_corrh_combine.cpp -o src/template_corrh_combine $(INCLUDES) $(LIBFLAGS)

src/template_conc: src/template_conc.cpp src/conc/template_conc.h src/enfw/enfw.h src/cosmology.h
	$(CPP) -fopenmp src/template_conc.cpp -o src/template_conc $(INCLUDES) $(LIBFLAGS) 

src/template_ell: src/template_ell.cpp src/enfw/template_ell.h src/enfw/enfw.h src/cosmology.h src/filter/filter.o
	$(CPP) -fopenmp src/template_ell.cpp src/filter/filter.o -o src/template_ell $(INCLUDES) $(LIBFLAGS) 


#### model

model_corrh: 
	@echo "========== checking for availability of corrh model =========="
	bash helpers/model_corrh.sh $(REDSHIFTS) $(ANNULI)
	@echo "========== finished with corrh model =========="

model_conc:
	@echo "========== checking for availability of conc model =========="
	bash helpers/model_conc.sh $(REDSHIFTS) $(ANNULI) $(CORES)
	@echo "========== finished with conc model =========="

model_ell:
	@echo "========== checking for availability of ell model =========="
	bash helpers/model_ell.sh $(REDSHIFTS) $(ANNULI) $(CORES)
	@echo "========== finished with ell model =========="


### model software

src/resample_ell: src/resample_ell.cpp src/enfw/template_ell.h src/cosmology.h
	$(CPP) -fopenmp src/resample_ell.cpp -o src/resample_ell $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)
src/resample_ell_g: src/resample_ell_g.cpp src/enfw/template_ell.h src/cosmology.h
	$(CPP) -fopenmp src/resample_ell_g.cpp -o src/resample_ell_g $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)

src/resample_conc: src/resample_conc.cpp src/conc/template_conc.h src/cosmology.h
	$(CPP) -fopenmp src/resample_conc.cpp -o src/resample_conc $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)
src/resample_conc_g: src/resample_conc_g.cpp src/conc/template_conc.h src/cosmology.h
	$(CPP) -fopenmp src/resample_conc_g.cpp -o src/resample_conc_g $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)

src/resample_corrh: src/resample_corrh.cpp src/corrh/template_corrh.h src/cosmology.h
	$(CPP) -fopenmp src/resample_corrh.cpp -o src/resample_corrh $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)
src/resample_corrh_g: src/resample_corrh.cpp src/corrh/template_corrh.h src/cosmology.h
	$(CPP) -fopenmp src/resample_corrh_g.cpp -o src/resample_corrh_g $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)

src/getmodel: src/getmodel.cpp src/model/covariance.h src/cosmology.h src/enfw/enfw.h
	$(CPP) -fopenmp src/getmodel.cpp -o src/getmodel $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)
src/getmodel_g: src/getmodel_g.cpp src/model/covariance.h src/cosmology.h src/enfw/enfw.h
	$(CPP) -fopenmp src/getmodel_g.cpp -o src/getmodel_g $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV)

### filter library

src/filter/filter.o: src/filter/filter.cpp src/filter/filter.h
	$(CPP) -fopenmp -o src/filter/filter.o $(INCLUDES) $(LIBFLAGS) -c src/filter/filter.cpp

#### forget about model (and re-do later, e.g. if you have changed your annuli definition)

forget_model:
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
	
