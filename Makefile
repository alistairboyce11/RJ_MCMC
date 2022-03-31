#FFLAGS = -O2
#OBJECTS = minos_bran.o test_new_mineos.o
#
#myproject : $(OBJECTS)
#	gfortran -o myproject $(OBJECTS)
#
#test_new_mineos.o : src_mineos/test_new_mineos.f
#	gfortran -c src_mineos/test_new_mineos.f
#
#minos_bran.o : src_mineos/minos_bran.f
#	gfortran -c src_mineos/minos_bran.f
#
#
#
#F77c =  mpif90 -check all  # -check bounds
F77 =  mpif90 -O3 #-Wall -fbounds-check#-O3 #-check bounds
# F77c =  mpif90  -O3  -fbounds-check #-check bounds

# mpif90 -check all
#gfortran -fbounds-check 
#mpif90 -mcmodel=medium -shared-intel
CC = gcc # -O3

all  :  obj global 

	
# driver: params.h data_params.h
# 	$(F77) -o run RJ_MCMC_store.f90 \
# 	-L./ -lroutines -lm
#
# driver_2: params.h data_params.h
# 	$(F77) -o run RJ_MCMC_store_2.f90 \
# 	-L./ -lroutines -lm
#
test_combine: params.h data_params.h
	$(F77) -o run test_combine.f90 \
	-L./ -lroutines -lm
#
# test_new_mineos: params.h data_params.h
# 	$(F77) -o run test_new_mineos.f90 \
# 	-L./ -lroutines -lm
#
# test_dispersion_minos: params.h data_params.h
# 	$(F77) -o run test_dispersion_minos.f90 \
# 	-L./ -lroutines -lm
#
joint: params.h data_params.h
	$(F77) -o run RJ_MCMC_test_joint.f90 \
	-L./ -lroutines -lm

joint_prepare: params.h data_params.h
	$(F77) -o run RJ_MCMC_test_joint_prepare.f90 \
	-L./ -lroutines -lm

joint_invert: params.h data_params.h
	$(F77) -o run2 RJ_MCMC_test_joint_invert.f90 \
	-L./ -lroutines -lm

global: params.h data_params.h
	$(F77) -o run RJ_MCMC_test_nlay.f90 \
	-L./ -lroutines -lm

test_stuck: params.h data_params.h
	$(F77) -o run test_stuck_models.f90 \
	-L./ -lroutines -lm

test_timing: params.h data_params.h
	$(F77) -o run test_timing.f90 \
	-L./ -lroutines -lm

detail_stuck: params.h data_params.h
	$(F77) -o run detail_stuck_models.f90 \
	-L./ -lroutines -lm

synth: params.h data_params.h
	$(F77) -o run synth_profile.f90 \
	-L./ -lroutines -lm

obj: params.h data_params.h
	$(F77) -c src/whichcell_d.f90
	$(F77) -c src/combine.f90
	$(F77) -c src/minos_bran.f
	$(F77) -c src/dispersion_minos.f90
	$(F77) -c src/interp.f90
	$(F77) -c src/combine_linear.f90
#
#	# FOR ANIREC
#	$(F77) -c src/forward_anirec.f
#	$(CC) -c src/refft.c
#	$(F77) -c src/anirec.f
#	$(F77) -c src/eispack.f
#	$(F77) -c src/matrixops.f
	
#	ar -r libroutines.a  makevoro2c.o raymrx.o raydsp.o dispersion_2psi.o four1.o \
#	dis2qmodel.o whichcell_d.o refft.o anirec.o eispack.o forward_anirec.o misfitCC4.o minos_bran.o combine.o dispersion_minos.o \
#	matrixops.o readwrite.o spheror.o buildmod.o 
#	\rm ./*.o

	ar -r libroutines.a  whichcell_d.o minos_bran.o combine.o dispersion_minos.o interp.o combine_linear.o
	\rm ./*.o

#
# obj_combine: params.h data_params.h
# 	$(F77) -c src/combine.f90
# #
# #	# FOR ANIREC
# #	$(F77) -c src/forward_anirec.f
# #	$(CC) -c src/refft.c
# #	$(F77) -c src/anirec.f
# #	$(F77) -c src/eispack.f
# #	$(F77) -c src/matrixops.f
#
# #	ar -r libroutines.a  makevoro2c.o raymrx.o raydsp.o dispersion_2psi.o four1.o \
# #	dis2qmodel.o whichcell_d.o refft.o anirec.o eispack.o forward_anirec.o misfitCC4.o minos_bran.o combine.o dispersion_minos.o \
# #	matrixops.o readwrite.o spheror.o buildmod.o
# #	\rm ./*.o
#
# 	ar -r libroutines.a  combine.o
# 	\rm ./*.o
	
clean:
	/bin/rm *.o
