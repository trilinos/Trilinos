#TOP LEVEL MAKEFILE 
#====================================================================
# ------------------------
# | CVS File Information |
# ------------------------
#
# $RCSfile$
#
# $Author$
#
# $Date$
#
# $Revision$
#
#====================================================================
#

##############################################################################
# Top-level Makefile for Zoltan library and Zoltan driver zdrive.
#
# To compile the library, type
#    gmake <optional flags> targetmachine
# where
#    targetmachine = sun, sunos, smos, or tflop.
#
# To compile the driver, type
#    gmake <optional flags> targetmachine-dr
#
# The optional flags are:
#    NO_PARMETIS=1   Use to compile and link a version of Zoltan that runs
#                    without ParMETIS.
#    NO_NEMESIS=1    Use to compile and link a version of zdrive that runs
#                    without Nemesis (i.e., reads only Chaco files).
#    YES_FORTRAN=1   Use to compile and link a version of Zoltan that can
#                    be called from Fortran.
#    YES_JOSTLE=1    Use to compile and link a version of Zoltan that runs
#                    with the JOSTLE library.
#
# Example:
# To compile a version of zdrive for SunOS 4 that doesn't use PARMETIS or
# NEMESIS, type
#    gmake NO_PARMETIS=1 NO_NEMESIS=1 sunos-dr
##############################################################################

GDEFINES	= -DDEBUG

# If NO_PARMETIS, compile the ParMETIS interface and link ParMETIS into zdrive.
ifdef NO_PARMETIS
GDEFINES       += -DLB_NO_PARMETIS
endif

# If NO_NEMESIS, don't compile or link Nemesis, Exodus and Netcdf into zdrive.
ifdef NO_NEMESIS
GDEFINES       += -DLB_NO_NEMESIS
endif

# If YES_FORTRAN, compile the Fortran interface into the library.
ifdef YES_FORTRAN
INC_FORTRAN	= inc_fortran
endif

# If YES_JOSTLE, compile the JOSTLE interface and link JOSTLE into zdrive.
ifdef YES_JOSTLE
GDEFINES       += -DLB_JOSTLE
endif

TARGET		= zdrive
LIB_TARGET	= zoltan
LIB_NAME	= libzoltan.a

##### Paths to third party libraries (not really part of Zoltan)  #####
  
# If you already have ParMetis on your machine, give the path here.
# Otherwise you first need to make ParMetis from the source in
# the directory parmetis/ParMetis (type "make" over there),
# and then move those libraries to the appropriate Obj_* directory. 
PARMETISLIB_PATH = .

# If you already have Jostle on your machine, give the path here.
# Otherwise, you first need to get Jostle from Chris Walshaw and 
# install it on your system.
JOSTLELIB_PATH = .

THIRD_PARTY_LIB_PATH = $(PARMETISLIB_PATH) $(JOSTLIB_PATH)

# Makefile for load-balancing library.
############################################################################
#                  SOURCE CODE
############################################################################
#
#   Paths to submodules
#
SUB_VPATH	= ../all ../lb ../par ../rcb ../comm ../oct ../driver ../ch \
		  ../params ../parmetis ../timer ../ha ../reftree
SUBINCPATH	= -I../all -I../lb -I../par -I../rcb -I../comm -I../oct \
		  -I../params -I../driver -I../ch -I../parmetis -I../timer \
                  -I../ha -I../reftree


ALL_CSRC	= all_allo.c free_structure.c
ALL_INC		= all_const.h all_allo_const.h

LB_CSRC		= lb.c lb_util.c 
LB_INC		= lb_const.h lb.h lbi_const.h lb_util_const.h

RCB_CSRC	= rcb.c rcb_util.c box_assign.c point_assign.c \
                  create_proc_list.c
RCB_INC		= rcb_const.h create_proc_list.h

OCT_CSRC	= msg.c octant.c oct_util.c octupdate.c dfs.c costs.c migoct.c \
		  output.c migreg.c hilbert.c
OCT_INC		= msg_const.h octant_const.h oct_util_const.h \
		  octupdate_const.h \
		  octupdate.h dfs_const.h dfs.h migoct_const.h \
		  migreg.h hilbert_const.h

PARAMS_CSRC	= set_param.c assign_param_vals.c check_param.c print_params.c \
		  key_params.c free_params.c bind_param.c 
PARAMS_INC	= params_const.h


COMM_CSRC	= comm_create.c comm_do.c comm_do_reverse.c comm_destroy.c \
		  invert_map.c
COMM_INC	= comm.h

PAR_CSRC	= par_util.c
PAR_INC		= par_const.h par.h

PARMETIS_CSRC	= parmetis_jostle.c 
PARMETIS_INC	= parmetis_jostle_const.h

TIMER_CSRC	= timer.c
TIMER_INC	= timer_const.h

HA_CSRC         = divide_machine.c get_processor_name.c \
                  build_machine_desc.c
HA_INC          = ha_const.h

REFTREE_CSRC    = reftree_build.c reftree_part.c reftree_hash.c
REFTREE_INC     = reftree_const.h

# object files for the Fortran interface

ifdef YES_FORTRAN
FOBJ_FILES	= lb_user_const.o ztypes.o fwrap.o cwrap.o
endif

# source code and defines for driver
DR_CSRC		= dr_main.c dr_err.c dr_input.c dr_par_util.c \
		  dr_util.c dr_exoII_io.c dr_elem.c dr_loadbal.c \
		  dr_migrate.c dr_output.c dr_chaco_io.c dr_elem_util.c \
		  dr_maps.c dr_eval.c
DR_INC		= dr_const.h dr_input_const.h dr_util_const.h \
		  dr_err_const.h dr_par_util_const.h dr_elem_const.h \
		  dr_elem_util_const.h dr_maps_const.h dr_eval_const.h

CH_CSRC		= ch_dist_graph.c ch_input_geom.c ch_input_graph.c \
                  ch_input_read.c ch_init_dist.c

CH_INC		= ch_input.h ch_input_const.h ch_init_dist_const.h

# libraries and paths
ROOT_ARCHD	= /Net/proj/zoltan/arch

ifndef NO_NEMESIS
NEM_LIBS	= -lnemIc -lexoIIc -lnetcdf
endif

ifndef NO_PARMETIS
PARMETIS_LIBS   = -lparmetis -lmetis
endif

ifdef YES_JOSTLE
JOSTLE_LIBS     = -ljostle
endif

THIRD_PARTY_LIBS = $(PARMETIS_LIBS) $(JOSTLE_LIBS)

# machine independent routines
INC		= $(GR_INC) $(ID_INC) $(VX_INC) $(ALL_INC) $(PAR_INC) \
		  $(RCB_INC) $(COMM_INC) $(OCT_INC) $(PARAMS_INC) \
		  $(PARMETIS_INC) $(TIMER_INC) $(HA_INC) $(REFTREE_INC)
#LIB_SRC_C	= $(GR_CSRC) $(ID_CSRC) $(VX_CSRC) $(ALL_CSRC) 
LIB_SRC_C	= $(LB_CSRC) $(ALL_CSRC) $(RCB_CSRC) $(PAR_CSRC) $(COMM_CSRC) \
		  $(OCT_CSRC) $(PARAMS_CSRC) $(PARMETIS_CSRC) $(TIMER_CSRC) \
                  $(HA_CSRC) $(REFTREE_CSRC)
LIB_SRC_OBJ	= $(LIB_SRC_C:.c=.o) 

#
# lint defines
#
#LINT		= alint # SunOs4
LINT		= lint
LINT_INC_PATH	= -Idriver -Ich -Iall -Ilb -Ipar -Ircb -Iparams\
		  -Icomm -Ioct -Iparmetis -Itimer -Iha -Ireftree\
		  -I/Net/local/lang/SC3.0.1/include \
		  -I/Net/local/lang/SC3.0.1/ansi_include \
		  -I/Net/local/mpi/include \
		  -I/Net/local/mpi/mpich/include \
                  -I/Net/local/mpi/build/solaris/ch_p4/include \
		  -I/Net/proj/zoltan/arch/sun/include
LINT_FLAGS	= $(LINT_INC_PATH) $(GDEFINES)
LINT_TARGET	= llib-l$(LIB_TARGET).ln
SUBDIRS		= driver ch all lb par rcb comm oct parmetis params timer ha \
		  reftree
VPATH		= $(SUBDIRS)

##############################################################################
# sunmos defines
##############################################################################
SMOS_CC		= sicc
SMOS_DBG_FLAGS	= -g
SMOS_OPT_FLAGS	= -O
SMOS_RANLIB	= touch
SMOS_AR		= ar860 r

SMOS_CSRC	= $(LIB_SRC_C)
SMOS_INC	= $(INC)
SMOS_OBJ	= $(SMOS_CSRC:.c=.o)

SMOS_INC_PATH	= $(SUBINCPATH)

SMOS_CCOMP	= $(SMOS_CC) $(SMOS_DBG_FLAGS) $(SMOS_INC_PATH)

SMOS_DEFS	= $(GDEFINES) -DSMOS

# driver stuff
SMOS_ARCHD	= $(ROOT_ARCHD)/smos
SMOS_DR_CSRC	= $(DR_CSRC) $(CH_CSRC)
SMOS_DR_OBJ	= $(SMOS_DR_CSRC:.c=.o)
SMOS_DR_IPATH	= -I$(SMOS_ARCHD)/include
#SMOS_MPI_LIB	= -L/Net/local/mpi/lib/paragon/sunmos -lmpi
SMOS_MPI_LIB	= -L/Net/mp/mpi/lib/paragon/ch_nx -lmpi
SMOS_TP_LPATH  = -L$(THIRD_PARTY_LIB_PATH) 
SMOS_DR_LIBS	= -L. -l$(LIB_TARGET) -L$(SMOS_ARCHD)/lib $(NEM_LIBS) \
		   $(SMOS_TP_LPATH) $(THIRD_PARTY_LIBS) $(SMOS_MPI_LIB)

##############################################################################
# Intel P6 defines
##############################################################################
#
TFLOP_CC	= pgcc -cougar
TFLOP_DBG_FLAGS	= -g
TFLOP_OPT_FLAGS	= -O3 -Minfo=all
TFLOP_RANLIB	= xranlib
TFLOP_AR	= xar r

TFLOP_CSRC	= $(LIB_SRC_C)
TFLOP_INC	= $(INC)
TFLOP_OBJ	= $(TFLOP_CSRC:.c=.o)

TFLOP_INC_PATH	= $(SUBINCPATH)

TFLOP_CCOMP	= $(TFLOP_CC) $(TFLOP_INC_PATH) $(TFLOP_DBG_FLAGS)

TFLOP_DEFS	= $(GDEFINES) -DTFLOP

# driver stuff
TFLOP_ARCHD	= $(ROOT_ARCHD)/tflop
TFLOP_DR_CSRC	= $(DR_CSRC) $(CH_CSRC)
TFLOP_DR_OBJ	= $(TFLOP_DR_CSRC:.c=.o)
TFLOP_DR_IPATH	= -I$(TFLOP_ARCHD)/include
TFLOP_MPI_LIB	= -lmpi
TFLOP_TP_LPATH = -L$(THIRD_PARTY_LIB_PATH)
TFLOP_DR_LIBS	= -L. -l$(LIB_TARGET) -L$(TFLOP_ARCHD)/lib $(NEM_LIBS) \
	          $(TFLOP_TP_LPATH) $(THIRD_PARTY_LIBS) $(TFLOP_MPI_LIB)

#
##############################################################################
# SUN defines (Solaris)
##############################################################################
#SUN_CC          = acc -vc -Xa -DUSE_STDARG -DHAVE_STDLIB_H=1 \
#                -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1 \
#                -DUSE_STDARG=1 -DMALLOC_RET_VOID=1
#SUN_CC		= purify -follow-child-processes=yes -cache-dir=/tmp/purify gcc
SUN_CC		= gcc
SUN_DBG_FLAGS	= -g 
SUN_OPT_FLAGS	= -O
SUN_RANLIB	= ranlib
SUN_AR		= ar r

SUN_CSRC	= $(LIB_SRC_C)
SUN_INC		= $(INC)
SUN_OBJ		= $(SUN_CSRC:.c=.o)

SUN_INC_PATH	= -I/Net/local/lang/SC3.0.1/include \
		  -I/Net/local/lang/SC3.0.1/ansi_include \
		  -I/Net/local/mpi/include \
                  -I/Net/local/mpi/build/solaris/ch_p4/include \
		  $(SUBINCPATH)

SUN_CCOMP	= $(SUN_CC) $(SUN_INC_PATH) $(SUN_DBG_FLAGS)

SUN_TP_LPATH   = -L$(THIRD_PARTY_LIB_PATH)

SUN_DEFS	= $(GDEFINES)
# driver stuff
SUN_ARCHD	= $(ROOT_ARCHD)/solaris
SUN_MPI_LIBS	= -L/Net/local/mpi/build/solaris/ch_p4/lib \
		  -lpmpich -lmpich -lsocket -lnsl -laio \
		  -L/Net/local/gnu/lib/gcc-lib/sparc-sun-solaris2.6/2.8.1 -lgcc

SUN_DR_CSRC	= $(DR_CSRC) $(CH_CSRC)
SUN_DR_OBJ	= $(SUN_DR_CSRC:.c=.o)
SUN_DR_IPATH	= -I$(SUN_ARCHD)/include

SUN_DR_LIBS	= -L. -l$(LIB_TARGET) -L$(SUN_ARCHD)/lib $(NEM_LIBS) \
	          $(SUN_TP_LPATH) $(THIRD_PARTY_LIBS) $(SUN_MPI_LIBS) -lm 
		  
##############################################################################
# SUNOS Definitions
##############################################################################

#SUNOS_CC          = purify -follow-child-processes=yes acc-2.0.1 -vc -Xa -DUSE_STDARG -DHAVE_STDLIB_H=1 \
#                -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1 \
#                -DUSE_STDARG=1 -DMALLOC_RET_VOID=1

SUNOS_CC          = acc-2.0.1 -vc -Xa -DUSE_STDARG -DHAVE_STDLIB_H=1 \
                -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1 \
                -DUSE_STDARG=1 -DMALLOC_RET_VOID=1
SUNOS_DBG_FLAGS   = -g
SUNOS_OPT_FLAGS   = -O4
SUNOS_RANLIB      = ranlib
SUNOS_AR          = ar r

SUNOS_CSRC        = $(LIB_SRC_C)
SUNOS_INC         = $(INC)
SUNOS_OBJ         = $(SUNOS_CSRC:.c=.o)
SUNOS_INC_PATH	= -I/Net/local/lang/SC2.0.1/include \
		  -I/Net/local/lang/SC2.0.1/ansi_include \
		  -I/Net/local/mpi/mpich/include \
                  -I/Net/local/mpi/mpich/lib/sun4/ch_p4 \
		  $(SUBINCPATH)
SUNOS_TP_LPATH   = -L$(THIRD_PARTY_LIB_PATH)
SUNOS_DEFS	= $(GDEFINES)

SUNOS_CCOMP	= $(SUNOS_CC) $(SUNOS_INC_PATH) $(SUNOS_DBG_FLAGS)

# driver stuff
SUNOS_ARCHD	= $(ROOT_ARCHD)/sun
SUNOS_MPI_LIBS	= -L/Net/local/mpi/mpich/lib/sun4/ch_p4 -lpmpich -lmpich \
		  -L/Net/local/gnu/lib/gcc-lib/sparc-sun-sunos4.1.3_U1/2.6.3 \
		  -lgcc

SUNOS_DR_CSRC	= $(DR_CSRC) $(CH_CSRC)
SUNOS_DR_OBJ	= $(SUNOS_DR_CSRC:.c=.o)
SUNOS_DR_IPATH	= -I$(SUNOS_ARCHD)/include

SUNOS_DR_LIBS	= -L. -l$(LIB_TARGET) -L$(SUNOS_ARCHD)/lib $(NEM_LIBS) \
	          $(SUNOS_TP_LPATH) $(THIRD_PARTY_LIBS) $(SUNOS_MPI_LIBS) -lm 

##############################################################################
# CPlant Definitions
##############################################################################
CPLANT_CC          = /Net/mp/cplant/current/bin/cc
CPLANT_DBG_FLAGS   = -g
CPLANT_OPT_FLAGS   = -O
CPLANT_RANLIB      = ranlib
CPLANT_AR          = ar r

CPLANT_CSRC        = $(LIB_SRC_C)
CPLANT_INC         = $(INC)
CPLANT_OBJ         = $(CPLANT_CSRC:.c=.o)
CPLANT_INC_PATH	= -I/Net/mp/cplant/current/include \
		  $(SUBINCPATH)
CPLANT_TP_LPATH   = -L$(THIRD_PARTY_LIB_PATH)
CPLANT_DEFS	= $(GDEFINES)

CPLANT_CCOMP	= $(CPLANT_CC) $(CPLANT_INC_PATH) $(CPLANT_DBG_FLAGS)

# driver stuff
CPLANT_ARCHD	= $(ROOT_ARCHD)/alpha
CPLANT_MPI_LIBS	= -L/Net/mp/cplant/current/lib -lmpi

CPLANT_DR_CSRC	= $(DR_CSRC) $(CH_CSRC)
CPLANT_DR_OBJ	= $(CPLANT_DR_CSRC:.c=.o)
CPLANT_DR_IPATH	= -I$(CPLANT_ARCHD)/include

CPLANT_DR_LIBS	= -L. -l$(LIB_TARGET) -L$(CPLANT_ARCHD)/lib $(NEM_LIBS) \
	          $(CPLANT_TP_LPATH) $(THIRD_PARTY_LIBS) $(CPLANT_MPI_LIBS) -lm

##############################################################################
#
# Optional parameters for the NIST targets.  These will not be documented
# anywhere but here, unless they find more widespread use.  Valid values
# are the character strings listed after the description of the parameter.
#
# MPI_LIB       Determines which MPI library to use if more than one is avail
#   lam
#   mpich  (default)
#
# F90_COMPILER  Determines which f90 compiler to use on a multicompiler machine
#   absoft
#   fujitsu
#   lahey
#   nag    (default)
#   nas
#   pgi
#   psr
#
##############################################################################

##############################################################################
# NIST SUN defines (Solaris)
##############################################################################
NISTSUN_F90 = f90
NISTSUN_FFLAGS = -g
NISTSUN_F90CFLAGS = -DFMANGLE=UNDERSCORE
NISTSUN_SPPR_HEAD = spprinc.sun
NISTSUN_FARG = farg_typical
NISTSUN_ZOLTAN_MOD_DIR = -M../fort

NISTSUN_CC		= gcc
NISTSUN_DBG_FLAGS	= -g 
NISTSUN_OPT_FLAGS	= -O
NISTSUN_RANLIB		= ranlib
NISTSUN_AR		= ar r

NISTSUN_CSRC		= $(LIB_SRC_C)
NISTSUN_INC		= $(INC)
NISTSUN_OBJ		= $(NISTSUN_CSRC:.c=.o)

ifeq "$(MPI_LIB)" "lam"
# For LAM
NISTSUN_INC_PATH	= -I$(LAMHOME)/include \
		  	$(SUBINCPATH)

NISTSUN_MPI_LIBS	= -L$(LAMHOME)/lib \
		  	-lmpi -largs -ltrillium -ltstdio -lt -lsocket \
		  	-L/itl/links/gnu/lib/gcc-lib/sparc-sun-solaris2.7/2.95.1 -lgcc
# end For LAM
else

# For MPICH (default if MPI_LIB is not defined)
NISTSUN_INC_PATH	= -I$(MPIR_HOME)/include \
		  	$(SUBINCPATH)

NISTSUN_MPI_LIBS	= -L$(MPIR_HOME)/build/solaris/ch_p4/lib \
		  	-lpmpich -lmpich -lsocket -lnsl -laio \
		  	-L/itl/links/gnu/lib/gcc-lib/sparc-sun-solaris2.7/2.95.1 -lgcc
# end For MPICH, my compilation
endif

NISTSUN_CCOMP		= $(NISTSUN_CC) $(NISTSUN_INC_PATH) \
			$(NISTSUN_DBG_FLAGS) $(NISTSUN_F90CFLAGS)

NISTSUN_TP_LPATH   	= -L$(THIRD_PARTY_LIB_PATH)

NISTSUN_DEFS		= $(GDEFINES) 

# driver stuff
NISTSUN_ARCHD		=

NISTSUN_DR_CSRC		= $(DR_CSRC) $(CH_CSRC)
NISTSUN_DR_OBJ		= $(NISTSUN_DR_CSRC:.c=.o)
NISTSUN_DR_IPATH	= -I$(NISTSUN_ARCHD)/include

NISTSUN_DR_LIBS		= -L. -l$(LIB_TARGET) -L$(NISTSUN_ARCHD)/lib \
			 $(NEM_LIBS) $(NISTSUN_TP_LPATH) $(THIRD_PARTY_LIBS) \
			 $(NISTSUN_MPI_LIBS) -lm 

##############################################################################
# NIST Linux defines
##############################################################################

# Standard definitions for this section

# Compiler dependent macros

ifeq "$(F90_COMPILER)" "absoft"
NISTLINUX_F90 = abf90
NISTLINUX_FFLAGS = -g -YEXT_NAMES=LCS -YCFRL=1 -YEXT_SFX=_ -lU77
NISTLINUX_F90CFLAGS = -DFMANGLE=UNDERSCORE
NISTLINUX_SPPR_HEAD = spprinc.absoft
NISTLINUX_FARG = farg_typical
NISTLINUX_ZOLTAN_MOD_DIR = -p ../fort
else
ifeq "$(F90_COMPILER)" "fujitsu"
NISTLINUX_F90 = f95
NISTLINUX_FFLAGS = -g -Am -X9 -static-flib
NISTLINUX_F90CFLAGS = -DFMANGLE=UNDERSCORE -DFUJITSU
NISTLINUX_SPPR_HEAD = spprinc.most
NISTLINUX_FARG = farg_typical
NISTLINUX_ZOLTAN_MOD_DIR = -I../fort
else
ifeq "$(F90_COMPILER)" "lahey"
NISTLINUX_F90 = lf95
NISTLINUX_FFLAGS = -g --staticlink
NISTLINUX_F90CFLAGS = -DFMANGLE=UNDERSCORE -DFUJITSU
NISTLINUX_SPPR_HEAD = spprinc.most
NISTLINUX_FARG = farg_typical
NISTLINUX_ZOLTAN_MOD_DIR = -I../fort
else
ifeq "$(F90_COMPILER)" "nas"
NISTLINUX_F90 = nasf95
NISTLINUX_FFLAGS = -g -nocheckint
NISTLINUX_F90CFLAGS = -DFMANGLE=UNDERSCORE
NISTLINUX_SPPR_HEAD = spprinc.nas
NISTLINUX_FARG = farg_nas
NISTLINUX_ZOLTAN_MOD_DIR =
else
ifeq "$(F90_COMPILER)" "pgi"
NISTLINUX_F90 = pgf90
NISTLINUX_FFLAGS =
NISTLINUX_F90CFLAGS = -DFMANGLE=UNDERSCORE -DPGI
NISTLINUX_SPPR_HEAD = spprinc.most
NISTLINUX_FARG = farg_typical
NISTLINUX_ZOLTAN_MOD_DIR = -I../fort
else
ifeq "$(F90_COMPILER)" "psr"
NISTLINUX_F90 = psrf90
NISTLINUX_FFLAGS = -g
NISTLINUX_F90CFLAGS = -DFMANGLE=UNDERSCORE
NISTLINUX_SPPR_HEAD = spprinc.most
NISTLINUX_FARG = farg_typical
# not sure, -H may only be for include files, not module information files
NISTLINUX_ZOLTAN_MOD_DIR = -H../fort
else
# default is NAG
NISTLINUX_F90 = nagf95
NISTLINUX_FFLAGS = -g -mismatch -gline -C
NISTLINUX_F90CFLAGS = -DFMANGLE=UNDERSCORE
NISTLINUX_SPPR_HEAD = spprinc.most
NISTLINUX_FARG = farg_nagf95
NISTLINUX_ZOLTAN_MOD_DIR = -I../fort
endif
endif
endif
endif
endif
endif
NISTLINUX_CC		= gcc
NISTLINUX_DBG_FLAGS	= -g 
NISTLINUX_OPT_FLAGS	= -O
NISTLINUX_RANLIB	= ranlib
NISTLINUX_AR		= ar r

NISTLINUX_CSRC		= $(LIB_SRC_C)
NISTLINUX_INC		= $(INC)
NISTLINUX_OBJ		= $(NISTLINUX_CSRC:.c=.o)

ifeq "$(MPI_LIB)" "lam"
# For LAM
NISTLINUX_INC_PATH	= -I$(LAMHOME)/include \
		  	$(SUBINCPATH)

NISTLINUX_MPI_LIBS	= -L$(LAMHOME)/lib \
		  	-lmpi -largs -ltrillium -ltstdio -lt \
		  	-L/usr/lib/gcc-lib/i386-redhat-linux/egcs-2.91.66 -lgcc
# end For LAM

else

# For MPICH (default if MPI_LIB is not defined)
NISTLINUX_INC_PATH	= -I$(MPIR_HOME)/include \
		  	$(SUBINCPATH)

NISTLINUX_MPI_LIBS	= -L$(MPIR_HOME)/build/LINUX/ch_p4/lib \
		  	-lpmpich -lmpich \
		  	-L/usr/lib/gcc-lib/i386-redhat-linux/egcs-2.91.66 -lg2c
# end For MPICH
endif

NISTLINUX_CCOMP		= $(NISTLINUX_CC) $(NISTLINUX_INC_PATH) \
			$(NISTLINUX_DBG_FLAGS) $(NISTLINUX_F90CFLAGS)

NISTLINUX_TP_LPATH   	= -L$(THIRD_PARTY_LIB_PATH)

NISTLINUX_DEFS		= $(GDEFINES)

# driver stuff
NISTLINUX_ARCHD		=

NISTLINUX_DR_CSRC	= $(DR_CSRC) $(CH_CSRC)
NISTLINUX_DR_OBJ	= $(NISTLINUX_DR_CSRC:.c=.o)
NISTLINUX_DR_IPATH	= -I$(NISTLINUX_ARCHD)/include

NISTLINUX_DR_LIBS	= -L. -l$(LIB_TARGET) -L$(NISTLINUX_ARCHD)/lib \
			 $(NEM_LIBS) $(NISTLINUX_TP_LPATH) $(THIRD_PARTY_LIBS) \
			 $(NISTLINUX_MPI_LIBS) -lm 

###########################################################################
#      HELP HELP HELP HELP HELP
###########################################################################
help:
	@echo ' '
	@echo 'To build the Zoltan library only, type "gmake target";'
	@echo 'to build the Zoltan library and the driver zdrive,'
	@echo 'type "gmake target-dr", where target is one of the following:'
	@echo ' '
	@echo '      smos       (for Intel Paragon)'
	@echo '      sun        (for Solaris with MPI; same as solaris)'
	@echo '      solaris    (for Solaris with MPI; same as sun)'
	@echo '      sunos      (for SunOS4 with MPI)'
	@echo '      tflop      (for Intel/Sandia ASCI Red)'
	@echo '      cplant     (for Computational Plant)'
	@echo ' '
	@echo 'Accepted compilation flags are'
	@echo '      NO_PARMETIS=1  (to exclude ParMetis from library)'
	@echo '      NO_NEMESIS=1   (to exclude Nemesis from driver)'
	@echo '      YES_FORTRAN=1  (to compile fortran interface)'
	@echo '      YES_JOSTLE=1   (to include Jostle in the library)'
	@echo ' '

###########################################################################
#      TARGETS FOR DIFFERENT MACHINES
###########################################################################

smos:
	@if [ ! -d Obj_smos ]; \
	then mkdir Obj_smos; \
	fi
	@rm -f Obj_smos/Makefile
	@cp Makefile_sub Obj_smos/Makefile
	@cd Obj_smos; \
	$(MAKE) $(LIB_NAME) "CC=$(SMOS_CCOMP)"\
	"LIB_NAME=$(LIB_NAME)" "INC_PATH=$(SMOS_INC_PATH)"\
	"AR=$(SMOS_AR)" "OBJ_FILES=$(SMOS_OBJ)" \
	"DEFS=$(SMOS_DEFS)" "VPATH=$(SUB_VPATH)" \
	"RANLIB=$(SMOS_RANLIB)"

smos-dr:
	@if [ ! -d Obj_smos ]; \
	then mkdir Obj_smos; \
	fi
	@rm -f Obj_smos/Makefile
	@cp Makefile_sub Obj_smos/Makefile
	@cd Obj_smos; \
	$(MAKE) $(TARGET) "TARGET=$(TARGET)" \
	"CC=$(SMOS_CCOMP) $(SMOS_DR_IPATH)" "LIB_NAME=$(LIB_NAME)" \
	"AR=$(SMOS_AR)" "OBJ_FILES=$(SMOS_OBJ)" \
	"DEFS=$(SMOS_DEFS)" "VPATH=$(SUB_VPATH)" \
	"RANLIB=$(SMOS_RANLIB)" "MAIN_OBJ=$(SMOS_DR_OBJ)" \
	"LNK_LIBS=$(SMOS_DR_LIBS)"
	@if [ ! -f test/$(TARGET).smos ]; \
	then ln -s ../Obj_smos/$(TARGET) test/$(TARGET).smos; \
	fi

tflop:
	@if [ ! -d Obj_tflop ]; \
	then mkdir Obj_tflop; \
	fi
	@rm -f Obj_tflop/Makefile
	@cp Makefile_sub Obj_tflop/Makefile
	@cd Obj_tflop; \
	$(MAKE) $(LIB_NAME) "CC=$(TFLOP_CCOMP)"\
	"LIB_NAME=$(LIB_NAME)" "INC=$(TFLOP_INC)"\
	"AR=$(TFLOP_AR)" "OBJ_FILES=$(TFLOP_OBJ)" \
	"DEFS=$(TFLOP_DEFS)" "VPATH=$(SUB_VPATH)" \
	"RANLIB=$(TFLOP_RANLIB)"

tflop-dr:
	@if [ ! -d Obj_tflop ]; \
	then mkdir Obj_tflop; \
	fi
	@rm -f Obj_tflop/Makefile
	@cp Makefile_sub Obj_tflop/Makefile
	@cd Obj_tflop; \
	$(MAKE) $(TARGET) "TARGET=$(TARGET)" \
	"CC=$(TFLOP_CCOMP) $(TFLOP_DR_IPATH)" "LIB_NAME=$(LIB_NAME)" \
	"AR=$(TFLOP_AR)" "OBJ_FILES=$(TFLOP_OBJ)" \
	"DEFS=$(TFLOP_DEFS)" "VPATH=$(SUB_VPATH)" \
	"RANLIB=$(TFLOP_RANLIB)" "MAIN_OBJ=$(TFLOP_DR_OBJ)" \
	"LNK_LIBS=$(TFLOP_DR_LIBS)"
	@if [ ! -f test/$(TARGET).tflop ]; \
	then ln -s ../Obj_tflop/$(TARGET) test/$(TARGET).tflop; \
	fi

solaris: sun

sun:
	@if [ ! -d Obj_sun ]; \
	then mkdir Obj_sun; \
	fi
	@rm -f Obj_sun/Makefile
	@cp Makefile_sub Obj_sun/Makefile
	@cd Obj_sun; \
	$(MAKE) $(LIB_NAME) "CC=$(SUN_CCOMP)"\
	"LIB_NAME=$(LIB_NAME)" "INC=$(SUN_INC)"\
	"AR=$(SUN_AR)" "OBJ_FILES=$(SUN_OBJ)" \
	"DEFS=$(SUN_DEFS)" "VPATH=$(SUB_VPATH)" \
	"RANLIB=$(SUN_RANLIB)"
#	@cscope -b $(MAIN) $(SUN_CSRC) $(SUN_INC)

solaris-dr: sun-dr

sun-dr:
	@if [ ! -d Obj_sun ]; \
	then mkdir Obj_sun; \
	fi
	@rm -f Obj_sun/Makefile
	@cp Makefile_sub Obj_sun/Makefile
	@cd Obj_sun; \
	$(MAKE) $(TARGET) "TARGET=$(TARGET)" \
	"CC=$(SUN_CCOMP) $(SUN_DR_IPATH)" "LIB_NAME=$(LIB_NAME)" \
	"AR=$(SUN_AR)" "OBJ_FILES=$(SUN_OBJ)" \
	"DEFS=$(SUN_DEFS)" "VPATH=$(SUB_VPATH)" \
	"RANLIB=$(SUN_RANLIB)" "MAIN_OBJ=$(SUN_DR_OBJ)" \
	"LNK_LIBS=$(SUN_DR_LIBS)"
	@if [ ! -f test/$(TARGET).sun ]; \
	then ln -s ../Obj_sun/$(TARGET) test/$(TARGET).sun; \
	fi

sunos:
	@if [ ! -d Obj_sunos ]; \
	then mkdir Obj_sunos; \
	fi
	@rm -f Obj_sunos/Makefile
	@cp Makefile_sub Obj_sunos/Makefile
	@cd Obj_sunos; \
	$(MAKE) $(LIB_NAME) "CC=$(SUNOS_CCOMP)"\
	"LIB_NAME=$(LIB_NAME)" "INC=$(SUNOS_INC)"\
	"AR=$(SUNOS_AR)" "OBJ_FILES=$(SUNOS_OBJ)" \
	"DEFS=$(SUNOS_DEFS)" "VPATH=$(SUB_VPATH)" \
	"RANLIB=$(SUNOS_RANLIB)"
#	@cscope -b $(MAIN) $(SUNOS_CSRC) $(SUNOS_INC)

sunos-dr:
	@if [ ! -d Obj_sunos ]; \
	then mkdir Obj_sunos; \
	fi
	@rm -f Obj_sunos/Makefile
	@cp Makefile_sub Obj_sunos/Makefile
	@cd Obj_sunos; \
	$(MAKE) $(TARGET) "TARGET=$(TARGET)" \
	"CC=$(SUNOS_CCOMP) $(SUNOS_DR_IPATH)" "LIB_NAME=$(LIB_NAME)" \
	"AR=$(SUNOS_AR)" "OBJ_FILES=$(SUNOS_OBJ)" \
	"DEFS=$(SUNOS_DEFS)" "VPATH=$(SUB_VPATH)" \
	"RANLIB=$(SUNOS_RANLIB)" "MAIN_OBJ=$(SUNOS_DR_OBJ)" \
	"LNK_LIBS=$(SUNOS_DR_LIBS)"
	@if [ ! -f test/$(TARGET).sunos ]; \
	then ln -s ../Obj_sunos/$(TARGET) test/$(TARGET).sunos; \
	fi

cplant:
	@if [ ! -d Obj_cplant ]; \
	then mkdir Obj_cplant; \
	fi
	@rm -f Obj_cplant/Makefile
	@cp Makefile_sub Obj_cplant/Makefile
	@cd Obj_cplant; \
	$(MAKE) $(LIB_NAME) "CC=$(CPLANT_CCOMP)"\
	"LIB_NAME=$(LIB_NAME)" "INC=$(CPLANT_INC)"\
	"AR=$(CPLANT_AR)" "OBJ_FILES=$(CPLANT_OBJ)" \
	"DEFS=$(CPLANT_DEFS)" "VPATH=$(SUB_VPATH)" \
	"RANLIB=$(CPLANT_RANLIB)"

cplant-dr:
	@if [ ! -d Obj_cplant ]; \
	then mkdir Obj_cplant; \
	fi
	@rm -f Obj_cplant/Makefile
	@cp Makefile_sub Obj_cplant/Makefile
	@cd Obj_cplant; \
	$(MAKE) $(TARGET) "TARGET=$(TARGET)" \
	"CC=$(CPLANT_CCOMP) $(CPLANT_DR_IPATH)" "LIB_NAME=$(LIB_NAME)" \
	"AR=$(CPLANT_AR)" "OBJ_FILES=$(CPLANT_OBJ)" \
	"DEFS=$(CPLANT_DEFS)" "VPATH=$(SUB_VPATH)" \
	"RANLIB=$(CPLANT_RANLIB)" "MAIN_OBJ=$(CPLANT_DR_OBJ)" \
	"LNK_LIBS=$(CPLANT_DR_LIBS)"
	@if [ ! -f test/$(TARGET).cplant ]; \
	then ln -s ../Obj_cplant/$(TARGET) test/$(TARGET).cplant; \
	fi

nistsun:
	@if [ ! -d Obj_nistsun ]; \
	then mkdir Obj_nistsun; \
	fi
	@rm -f Obj_nistsun/Makefile
	@cp Makefile_sub Obj_nistsun/Makefile
	@cd Obj_nistsun; \
	$(MAKE) $(LIB_NAME) "CC=$(NISTSUN_CCOMP)"\
	"LIB_NAME=$(LIB_NAME)" "INC=$(NISTSUN_INC)"\
	"AR=$(NISTSUN_AR)" "OBJ_FILES=$(NISTSUN_OBJ)" \
	"DEFS=$(NISTSUN_DEFS)" "VPATH=$(SUB_VPATH)" \
	"RANLIB=$(NISTSUN_RANLIB)" "INC_FORTRAN=$(INC_FORTRAN)" \
	"FOBJ_FILES=$(FOBJ_FILES)" "F90=$(NISTSUN_F90) $(NISTSUN_FFLAGS)" \
	"SPPR_HEAD=$(NISTSUN_SPPR_HEAD)"

nistsun-dr:
	@if [ ! -d Obj_nistsun ]; \
	then mkdir Obj_nistsun; \
	fi
	@rm -f Obj_nistsun/Makefile
	@cp Makefile_sub Obj_nistsun/Makefile
	@cd Obj_nistsun; \
	$(MAKE) $(TARGET) "TARGET=$(TARGET)" \
	"CC=$(NISTSUN_CCOMP) $(NISTSUN_DR_IPATH)" "LIB_NAME=$(LIB_NAME)" \
	"AR=$(NISTSUN_AR)" "OBJ_FILES=$(NISTSUN_OBJ)" \
	"DEFS=$(NISTSUN_DEFS)" "VPATH=$(SUB_VPATH)" \
	"RANLIB=$(NISTSUN_RANLIB)" "INC_FORTRAN=$(INC_FORTRAN)" \
	"FOBJ_FILES=$(FOBJ_FILES)" "MAIN_OBJ=$(NISTSUN_DR_OBJ)" \
	"LNK_LIBS=$(NISTSUN_DR_LIBS)"
	@if [ ! -f test/$(TARGET).nistsun ]; \
	then ln -s ../Obj_nistsun/$(TARGET) test/$(TARGET).nistsun; \
	fi

nistsun-fdr:
	@cd fdriver; \
	$(MAKE) zfdrive "F90=$(NISTSUN_F90)" "FFLAGS=$(NISTSUN_FFLAGS)" \
	"CC=$(NISTSUN_CCOMP)" "FARG=$(NISTSUN_FARG)" \
        "ZOLTAN_MOD_DIR=$(NISTSUN_ZOLTAN_MOD_DIR)" \
	"ZOLTAN_OBJ_DIR=../Obj_nistsun" "MPIINC=$(NISTSUN_INC_PATH)" \
	"MPILIB=$(NISTSUN_MPI_LIBS)" "TARGET_ARCH=nistsun"

nistlinux:
	@if [ ! -d Obj_nistlinux ]; \
	then mkdir Obj_nistlinux; \
	fi
	@rm -f Obj_nistlinux/Makefile
	@cp Makefile_sub Obj_nistlinux/Makefile
	@cd Obj_nistlinux; \
	$(MAKE) $(LIB_NAME) "CC=$(NISTLINUX_CCOMP)"\
	"LIB_NAME=$(LIB_NAME)" "INC=$(NISTLINUX_INC)"\
	"AR=$(NISTLINUX_AR)" "OBJ_FILES=$(NISTLINUX_OBJ)" \
	"DEFS=$(NISTLINUX_DEFS)" "VPATH=$(SUB_VPATH)" \
	"RANLIB=$(NISTLINUX_RANLIB)" "INC_FORTRAN=$(INC_FORTRAN)" \
	"FOBJ_FILES=$(FOBJ_FILES)" "F90=$(NISTLINUX_F90) $(NISTLINUX_FFLAGS)" \
	"SPPR_HEAD=$(NISTLINUX_SPPR_HEAD)"

nistlinux-dr:
	@if [ ! -d Obj_nistlinux ]; \
	then mkdir Obj_nistlinux; \
	fi
	@rm -f Obj_nistlinux/Makefile
	@cp Makefile_sub Obj_nistlinux/Makefile
	@cd Obj_nistlinux; \
	$(MAKE) $(TARGET) "TARGET=$(TARGET)" \
	"CC=$(NISTLINUX_CCOMP) $(NISTLINUX_DR_IPATH)" "LIB_NAME=$(LIB_NAME)" \
	"AR=$(NISTLINUX_AR)" "OBJ_FILES=$(NISTLINUX_OBJ)" \
	"DEFS=$(NISTLINUX_DEFS)" "VPATH=$(SUB_VPATH)" \
	"RANLIB=$(NISTLINUX_RANLIB)" "INC_FORTRAN=$(INC_FORTRAN)" \
	"FOBJ_FILES=$(FOBJ_FILES)" "MAIN_OBJ=$(NISTLINUX_DR_OBJ)" \
	"LNK_LIBS=$(NISTLINUX_DR_LIBS)"
	if [ ! -f test/$(TARGET).nistlinux ]; \
	then ln -s ../Obj_nistlinux/$(TARGET) test/$(TARGET).nistlinux; \
	fi

nistlinux-fdr:
	@cd fdriver; \
	$(MAKE) zfdrive "F90=$(NISTLINUX_F90)" "FFLAGS=$(NISTLINUX_FFLAGS)" \
	"CC=$(NISTLINUX_CCOMP)" "FARG=$(NISTLINUX_FARG)" \
	"ZOLTAN_MOD_DIR=$(NISTLINUX_ZOLTAN_MOD_DIR)" \
	"ZOLTAN_OBJ_DIR=../Obj_nistlinux" "MPIINC=$(NISTLINUX_INC_PATH)" \
        "MPILIB=$(NISTLINUX_MPI_LIBS)" "TARGET_ARCH=nistlinux"

#
############################################################################
#          UTILITY TARGETS
############################################################################
#

$(LINT_TARGET): $(LIB_SRC_C)
	@echo "Generating lint file for library..."
	$(LINT) $(LINT_FLAGS) $? -o $(LIB_TARGET)

lint: $(DR_CSRC) $(CH_CSRC) $(LINT_TARGET)
	@echo "running lint..."
	$(LINT) $(LINT_FLAGS) $?


printout:
	cat $(MAIN) $(LIB_SRC_C) > source.c
	tps -c source.c | lpr
	-rm source.c
	tps  $(INC) | lpr

clean:
	-for a in $(SUBDIRS); do cd $$a; \
	 $(RM) -f *.o *.d *.d@ zdrive zdrive.pure; \
	 cd ..; done;
	 cd fort; $(MAKE) clean;
	 cd fdriver; $(MAKE) clean; exit 0
