INCLUDES=-I ~/werc3/include 
LIBFLAGS=-L ~/werc3/lib -L ~/lib -lCCfits
LIBFLAGS_TMV=-ltmv -lblas -lpthread -ltmv_symband
CPP=g++
CORES=47
REDSHIFTS=0.24533 0.35 0.187 0.206 0.224 0.234 0.288 0.313 0.348 0.352 0.363 0.391 0.399 0.440 0.450 0.451 0.686
ANNULI=annuli_keiichi.tab

VERSION=0.1

all: software lut templates model

##### all

software: lut_software template_software model_software tinker			# programs to calculate covariances

lut: lut_2pc lut_W lut_U lut_sigmam lut_dndM lut_rho0 lut_profiles		# look-up tables

templates: templates_corrh templates_conc templates_ell				# templates for intrinsic covariance components

model: model_corrh model_conc model_ell						# co-added and resampled temples according to some binning scheme


#### software

template_software: src/template_corrh src/template_corrh_combine src/template_conc src/template_ell

model_software: src/resample_ell src/resample_conc src/resample_corrh src/getmodel

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
	$(CPP) -fopenmp -o src/template_corrh $(INCLUDES) $(LIBFLAGS) src/template_corrh.cpp

src/template_corrh_combine: src/template_corrh_combine.cpp src/corrh/template_corrh.h src/corrh/corrh.h src/enfw/enfw.h src/profile/profile.h src/cosmology.h
	$(CPP) -fopenmp -o src/template_corrh_combine $(INCLUDES) $(LIBFLAGS) src/template_corrh_combine.cpp

src/template_conc: src/template_conc.cpp src/conc/template_conc.h src/enfw/enfw.h src/cosmology.h
	$(CPP) -fopenmp -o src/template_conc $(INCLUDES) $(LIBFLAGS) src/template_conc.cpp

src/template_ell: src/template_ell.cpp src/enfw/template_ell.h src/enfw/enfw.h src/cosmology.h src/filter/filter.o
	$(CPP) -fopenmp -o src/template_ell $(INCLUDES) $(LIBFLAGS) src/template_ell.cpp src/filter/filter.o


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
	$(CPP) -o src/resample_ell $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV) src/resample_ell.cpp

src/resample_conc: src/resample_conc.cpp src/conc/template_conc.h src/cosmology.h
	$(CPP) -o src/resample_conc $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV) src/resample_conc.cpp

src/resample_corrh: src/resample_corrh.cpp src/corrh/template_corrh.h src/cosmology.h
	$(CPP) -o src/resample_corrh $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV) src/resample_corrh.cpp

src/getmodel: src/getmodel.cpp src/model/covariance.h src/cosmology.h src/enfw/enfw.h
	$(CPP) -o src/getmodel $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV) src/getmodel.cpp

### filter library

src/filter/filter.o: src/filter/filter.cpp src/filter/filter.h
	$(CPP) -fopenmp -o src/filter/filter.o $(INCLUDES) $(LIBFLAGS) -c src/filter/filter.cpp

### pack templates and make available online (to be run by Daniel...)

pub: pub/templates_conc.tar.gz pub/templates_ell.tar.gz pub/ccv.tar.gz
	rsync -avv pub/templates_conc.tar.gz pub/templates_ell.tar.gz templates/corrh_*fits lut/W*tab lut/dndM*tab dgruen@moon.usm.uni-muenchen.de:/usr/web/users/dgruen/public_html/code/templates/
	rsync -avv pub/ccv-$(VERSION).tar.gz dgruen@moon.usm.uni-muenchen.de:/usr/web/users/dgruen/public_html/code/

pub/templates_conc.tar.gz: templates/conc/cov_100.fits
	tar czf pub/templates_conc.tar.gz templates/conc/*.fits

pub/templates_ell.tar.gz: templates/ell/cov_100.fits
	tar czf pub/templates_ell.tar.gz templates/ell/*.fits

pub/ccv.tar.gz:
	git commit -a
	bash helpers/pack_software.sh $(VERSION)
	
