#-------------------------------------------------------------------------------
# AMD Makefile for compiling on Unix systems (for GNU make only)
#-------------------------------------------------------------------------------

all: ../Lib/libamd.a

include ../Make/Make.include

C = $(CC) $(CFLAGS) $(CONFIG) -I../Include

#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------

AMD = amd_aat amd_1 amd_2 amd_dump amd_postorder amd_post_tree amd_defaults \
	amd_order amd_control amd_info amd_valid amd_preprocess

INC = ../Include/amd.h amd_internal.h

#-------------------------------------------------------------------------------
# object files for each version
#-------------------------------------------------------------------------------

AMDI = $(addsuffix .o, $(subst amd_,amd_i_,$(AMD)))
AMDL = $(addsuffix .o, $(subst amd_,amd_l_,$(AMD)))

#-------------------------------------------------------------------------------
# compile each int and long routine (with no real/complex version)
#-------------------------------------------------------------------------------

amd_i_%.o: amd_%.c $(INC)
	$(C) -DDINT -c $< -o $@

amd_l_%.o: amd_%.c $(INC)
	$(C) -DDLONG -c $< -o $@

#-------------------------------------------------------------------------------
# Create the libamd.a library (C versions only)
#-------------------------------------------------------------------------------

../Lib/libamd.a: $(AMDI) $(AMDL)
	$(AR) cr ../Lib/libamd.a $^
	- $(RANLIB) ../Lib/libamd.a

#-------------------------------------------------------------------------------
# compile the Fortran versions and the libamdf77.a library
#-------------------------------------------------------------------------------

fortran: ../Lib/libamdf77.a

AMDF77 = amd.o amdbar.o

amd.o: amd.f
	$(F77) $(F77FLAGS) -c amd.f -o amd.o

amdbar.o: amdbar.f
	$(F77) $(F77FLAGS) -c amdbar.f -o amdbar.o

../Lib/libamdf77.a: $(AMDF77)
	$(AR) cr ../Lib/libamdf77.a $^
	- $(RANLIB) ../Lib/libamdf77.a

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

purge: clean
	- $(RM) ../Lib/libamd.a ../Lib/libamdf77.a
