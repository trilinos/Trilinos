#######################################################################
#
#  This makefile creates a library for distributed SuperLU.
#  The files are organized as follows:
#
#       ALLAUX  -- Auxiliary routines called from all precisions
#       DSLUSRC -- Double precision real serial SuperLU routines
#       DPLUSRC -- Double precision real parallel SuperLU routines
#       ZSLUSRC -- Double precision complex serial SuperLU routines
#       ZPLUSRC -- Double precision complex parallel SuperLU routines
#
#  The library can be set up to include routines for any combination
#  of the two precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make double
#       make double complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all two precisions.
#  The library is called
#       superlu.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#
#######################################################################
include ../make.inc
#
# Precision independent routines
#
ALLAUX 	= symbfact.o sp_ienv.o zfdperm.o \
	  etree.o sp_colorder.o get_perm_c.o \
	  mmd.o colamd.o comm.o memory.o util.o superlu_grid.o \
	  pxerbla.o superlu_timer.o

#
# Routines to permute large entries to diagonal
#
DLDPERM = dldperm.o mc64ad.o mc21.o #fd05.o

#
# Routines literally taken from SuperLU
#
DSLUSRC	= dlangs.o dgsequ.o dlaqgs.o dutil.o \
	  dmemory.o dsp_blas2.o dsp_blas3.o dreadhb.o
ZSLUSRC	= dcomplex.o dlamch.o slamch.o zlangs.o zgsequ.o zlaqgs.o \
	  zmemory.o zsp_blas2.o zsp_blas3.o zreadhb.o zutil.o

#
# Routines for double precision parallel SuperLU
DPLUSRC = ddistribute.o pdgstrf_irecv.o pdgstrs_Bglobal.o \
	  pdgstrs1.o pdgssvx_ABglobal.o pdgsrfs_ABXglobal.o pdgsmv.o
#
# Routines for double complex parallel SuperLU
ZPLUSRC = zldperm.o zdistribute.o pzgssvx_ABglobal.o pzgstrf_irecv.o \
	  pzgstrs1.o pzgstrs_Bglobal.o pzgsrfs_ABXglobal.o pzgsmv.o


all:  double complex16

double: $(DSLUSRC) $(DPLUSRC) $(DLDPERM) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) $(DSUPERLULIB) $(DSLUSRC) $(DPLUSRC) \
	$(DLDPERM) $(ALLAUX)
	$(RANLIB) $(DSUPERLULIB)

complex16: $(ZSLUSRC) $(ZPLUSRC) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) $(DSUPERLULIB) $(ZSLUSRC) $(ZPLUSRC) $(ALLAUX)
	$(RANLIB) $(DSUPERLULIB)


##################################
# Do not optimize these routines #
##################################
slamch.o: slamch.c ; $(CC) $(NOOPTS) $(CDEFS) -c $<
dlamch.o: dlamch.c ; $(CC) $(NOOPTS) $(CDEFS) -c $<
superlu_timer.o: superlu_timer.c ; $(CC) -c $(NOOPTS) $<


.c.o:
	$(CC) $(CFLAGS) $(CDEFS) $(BLASDEF) -c $< $(VERBOSE)

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

clean:	
	rm -f *.o
