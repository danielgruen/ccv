
LIB = -lm 

CC = gcc
CFLAGS = -O2
EXEC = massfunction.x

OBJ_HOD = header.o main.o utility.o sigmac.o transfnc.o transfunc_file.o \
	nonlinear_power_spectrum.o \
	input_params.o \
	tf_eisenstein_hu.o \
	dark_matter_statistics.o halo_mass_function.o mstar.o growthfactor.o

OBJ_NR = nrutil.o qromo.o midpnt.o midinf.o polint.o splint.o spline.o \
	qtrap.o zbrent.o trapzd.o

OBJS01 = $(OBJ_HOD) $(OBJ_NR)

$(EXEC): $(OBJS01)
	$(CC) -o $@ $(OBJS01) $(LIB)
$(OBJS01):	header.h nrutil.h

print_sigmac: $(OBJS01) print_sigmac.c
	$(CC) -o print_sigmac $(OBJ_NR) header.o utility.o sigmac.o transfnc.o transfunc_file.o nonlinear_power_spectrum.o input_params.o tf_eisenstein_hu.o halo_mass_function.o mstar.o growthfactor.o print_sigmac.c $(LIB)

sigma_of_r: $(OBJS01) sigma_of_r.c
	$(CC) -o sigma_of_r $(OBJ_NR) header.o utility.o sigmac.o transfnc.o transfunc_file.o nonlinear_power_spectrum.o input_params.o tf_eisenstein_hu.o halo_mass_function.o mstar.o growthfactor.o sigma_of_r.c $(LIB)

clean:
	rm -f *.o
