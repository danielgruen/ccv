INCLUDES=-I ~/werc3/include 
LIBFLAGS=-L ~/werc3/lib -L ~/lib -lCCfits
CPP=g++
CORES=15
REDSHIFT=0.24533
ZLONG=0.2453300 # must be 1.7 digits
ANNULI=annuli_keiichi.tab

LIBFLAGS_TMV=-ltmv -lblas -lpthread


all: software lut templates model

##### all

software: template_software model_software					# programs to calculate covariances

lut: lut_2pc lut_W lut_U lut_sigmam lut_bias lut_dndM lut_rho0 lut_profiles	# look-up tables

templates: templates_corrh templates_conc templates_ell				# templates for intrinsic covariance components

model: model_corrh model_conc model_ell						# co-added and resampled temples according to some binning scheme


#### software

template_software: src/template_corrh src/template_corrh_combine src/template_conc src/template_ell

model_software: src/resample_ell src/resample_conc src/resample_corrh


#### lut

lut_2pc:
	@echo "don't know how to make lut_2pc, but will ignore that"
lut_W:
	@echo "don't know how to make lut_W, but will ignore that"
lut_U:
	@echo "don't know how to make lut_U, but will ignore that"
lut_sigmam:
	@echo "don't know how to make lut_sigmam, but will ignore that"
lut_bias:
	@echo "don't know how to make lut_bias, but will ignore that"
lut_dndM:
	@echo "don't know how to make lut_dndM, but will ignore that"
lut_rho0:
	@echo "don't know how to make lut_rho0, but will ignore that"
lut_profiles:
	@echo "don't know how to make lut_profiles, but will ignore that"


#### templates

templates_corrh:
	@echo "========== checking for availability of corrh templates =========="
	bash helpers/templates_corrh.sh $(REDSHIFT) $(CORES)
	bash helpers/templates_corrh_combine.sh $(REDSHIFT)
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
	bash helpers/model_corrh.sh $(REDSHIFT) $(ANNULI)
	@echo "========== finished with corrh model =========="

model_conc:
	@echo "========== checking for availability of conc model =========="
	bash helpers/model_conc.sh $(REDSHIFT) $(ANNULI) $(CORES)
	@echo "========== finished with conc model =========="

model_ell:
	@echo "========== checking for availability of ell model =========="
	bash helpers/model_ell.sh $(REDSHIFT) $(ANNULI) $(CORES)
	@echo "========== finished with ell model =========="


### model software

src/resample_ell: src/resample_ell.cpp src/enfw/template_ell.h src/cosmology.h
	$(CPP) -o src/resample_ell $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV) src/resample_ell.cpp

src/resample_conc: src/resample_conc.cpp src/conc/template_conc.h src/cosmology.h
	$(CPP) -o src/resample_conc $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV) src/resample_conc.cpp

src/resample_corrh: src/resample_corrh.cpp src/corrh/template_corrh.h src/cosmology.h
	$(CPP) -o src/resample_corrh $(INCLUDES) $(LIBFLAGS) $(LIBFLAGS_TMV) src/resample_corrh.cpp

### filter library

src/filter/filter.o: src/filter/filter.cpp src/filter/filter.h
	$(CPP) -fopenmp -o src/filter/filter.o $(INCLUDES) $(LIBFLAGS) -c src/filter/filter.cpp