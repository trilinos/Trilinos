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
#
# Example:
# To compile a version of zdrive for SunOS 4 that doesn't use PARMETIS or
# NEMESIS, type
#    gmake NO_PARMETIS=1 NO_NEMESIS=1 sunos-dr
##############################################################################

GDEFINES	= -DDEBUG

# If NO_PARMETIS, don't compile ParMetis into the library.
ifdef NO_PARMETIS
GDEFINES       += -DLB_NO_PARMETIS
endif

# If NO_NEMESIS, don't compile Nemesis, Exodus and Netcdf into the library.
ifdef NO_NEMESIS
GDEFINES       += -DLB_NO_NEMESIS
endif

TARGET		= zdrive
LIB_TARGET	= zoltan
LIB_NAME	= libzoltan.a

SUPP_TARGS	= "smos sun solaris sunos tflop"

# Paths to third party libraries (not really part of Zoltan)
# 
# If you already have ParMetis on your machine, give the path here.
# Otherwise you need to type "make parmetis" to compile 
# the source that was distributed with Zoltan.
PARMETISLIB_PATH = parmetis/ParMetis

# Makefile for load-balancing library.
############################################################################
#                  SOURCE CODE
############################################################################
#
#   Paths to submodules
#
SUB_VPATH	= ../all ../lb ../par ../rcb ../comm ../oct ../driver ../ch \
		  ../parmetis
SUBINCPATH	= -I../all -I../lb -I../par -I../rcb -I../comm -I../oct \
		  -I../driver -I../ch -I../parmetis

#source code from SCCS projects

GR_CSRC		= gr.c gr_hash.c gr_list.c gr_tree.c gr_type.c gr_util.c
GR_INC		= gr_type.h gr_const.h gr_hash_const.h gr_list_const.h \
		  gr_tree_const.h gr_type_const.h \
		  gr_util_const.h

ID_CSRC		= id.c id_util.c
ID_INC		= id_const.h id_util.h id_util_const.h

VX_CSRC		= vx.c
VX_INC		= vx_const.h

ALL_CSRC	= all_allo.c
ALL_INC		= all_const.h all_allo_const.h

LB_CSRC		= lb.c lb_util.c 
LB_INC		= lb_const.h lb.h lbi_const.h lb_util_const.h

RCB_CSRC	= rcb.c rcb_util.c
RCB_INC		= rcb_const.h

OCT_CSRC	= msg.c octant.c util.c octupdate.c dfs.c costs.c migoct.c \
		  output.c migreg.c hilbert.c
OCT_INC		= msg_const.h octant_const.h util_const.h octupdate_const.h \
		  octupdate.h dfs_const.h dfs.h migoct_const.h \
		  migreg.h hilbert_const.h

COMM_CSRC	= comm.c
COMM_INC	= comm.h

PAR_CSRC	= par_util.c
PAR_INC		= par_const.h par.h

PARMETIS_CSRC	= parmetis_part.c parmetis_util.c 
PARMETIS_INC	= parmetis_const.h

# source code and defines for driver
DR_CSRC		= dr_main.c dr_err.c dr_input.c dr_par_util.c \
		  dr_util.c dr_exoII_io.c dr_elem.c dr_loadbal.c \
		  dr_migrate.c dr_output.c dr_chaco_io.c
DR_INC		= dr_const.h dr_input_const.h dr_util_const.h \
		  dr_err_const.h dr_par_util_const.h dr_elem_const.h

CH_CSRC		= ch_dist_graph.c ch_input_geom.c ch_input_graph.c \
                  ch_input_read.c

CH_INC		= ch_input.h ch_input_const.h

# libraries and paths
ROOT_ARCHD	= /Net/pecos/rf_crew/arch

ifndef NO_NEMESIS
NEM_LIBS	= -lnemIc -lexoIIc -lnetcdf
endif

ifndef NO_PARMETIS
PARMETIS_LIBS   = -lparmetis -lmetis
endif

# machine independent routines
INC		= $(GR_INC) $(ID_INC) $(VX_INC) $(ALL_INC) $(PAR_INC) \
		  $(RCB_INC) $(COMM_INC) $(OCT_INC) $(PARMETIS_INC)
#LIB_SRC_C	= $(GR_CSRC) $(ID_CSRC) $(VX_CSRC) $(ALL_CSRC) 
LIB_SRC_C	= $(LB_CSRC) $(ALL_CSRC) $(RCB_CSRC) $(PAR_CSRC) $(COMM_CSRC) \
		  $(OCT_CSRC) $(PARMETIS_CSRC)
LIB_SRC_OBJ	= $(LIB_SRC_C:.c=.o) 

#
# lint defines
#
#LINT		= alint # SunOs4
LINT		= lint
LINT_INC_PATH	= -Idriver -Ich -Iall -Ilb -Ipar -Ircb \
		  -Icomm -Ioct -Iparmetis \
		  -I/Net/local/lang/SC3.0.1/include \
		  -I/Net/local/lang/SC3.0.1/ansi_include \
		  -I/Net/local/mpi/include \
		  -I/Net/local/mpi/mpich/include \
                  -I/Net/local/mpi/build/solaris/ch_p4/include \
		  -I/Net/pecos/rf_crew/arch/sun/include
LINT_FLAGS	= $(LINT_INC_PATH) $(GDEFINES)
LINT_TARGET	= llib-l$(LIB_TARGET).ln
SUBDIRS		= driver ch all lb par rcb comm oct parmetis
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

SMOS_DEFS	= $(GDEFINES)

# driver stuff
SMOS_ARCHD	= $(ROOT_ARCHD)/smos
SMOS_DR_CSRC	= $(DR_CSRC) $(CH_CSRC)
SMOS_DR_OBJ	= $(SMOS_DR_CSRC:.c=.o)
SMOS_DR_IPATH	= -I$(SMOS_ARCHD)/include
SMOS_MPI_LIB	= -L/Net/local/mpi/lib/paragon/sunmos -lmpi
SMOS_PAR_LPATH  = -L$(PARMETISLIB_PATH)/Obj_smos \
                  -L../$(PARMETISLIB_PATH)/Obj_smos 
SMOS_DR_LIBS	= -L. -l$(LIB_TARGET) -L$(SMOS_ARCHD)/lib $(NEM_LIBS) \
		   $(SMOS_MPI_LIB) $(SMOS_PAR_LPATH) $(PARMETIS_LIBS)

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

TFLOP_DEFS	= $(GDEFINES)

# driver stuff
TFLOP_ARCHD	= $(ROOT_ARCHD)/tflop
TFLOP_DR_CSRC	= $(DR_CSRC) $(CH_CSRC)
TFLOP_DR_OBJ	= $(TFLOP_DR_CSRC:.c=.o)
TFLOP_DR_IPATH	= -I$(TFLOP_ARCHD)/include
TFLOP_MPI_LIB	= -lmpi
TFLOP_PAR_LPATH = -L$(PARMETISLIB_PATH)/Obj_tflop \
                  -L../$(PARMETISLIB_PATH)/Obj_tflop 
TFLOP_DR_LIBS	= -L. -l$(LIB_TARGET) -L$(TFLOP_ARCHD)/lib $(NEM_LIBS) \
	          $(TFLOP_MPI_LIB) $(TFLOP_PAR_LPATH) $(PARMETIS_LIBS)

#
##############################################################################
# SUN defines (Solaris)
##############################################################################
#SUN_CC          = acc -vc -Xa -DUSE_STDARG -DHAVE_STDLIB_H=1 \
#                -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1 \
#                -DUSE_STDARG=1 -DMALLOC_RET_VOID=1
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

SUN_CCOMP	= $(SUN_CC) $(SUN_INC_PATH) $(SUN_OPT_FLAGS)

SUN_PAR_LPATH   = -L$(PARMETISLIB_PATH)/Obj_sun \
                  -L../$(PARMETISLIB_PATH)/Obj_sun 

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
	          $(SUN_PAR_LPATH) $(PARMETIS_LIBS) $(SUN_MPI_LIBS) -lm 
		  
##############################################################################
# SUNOS Definitions
##############################################################################

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
SUNOS_PAR_LPATH   = -L$(PARMETISLIB_PATH)/Obj_sunos \
                  -L../$(PARMETISLIB_PATH)/Obj_sunos
SUNOS_DEFS	= $(GDEFINES)

SUNOS_CCOMP	= $(SUNOS_CC) $(SUNOS_INC_PATH) $(SUNOS_OPT_FLAGS)

# driver stuff
SUNOS_ARCHD	= $(ROOT_ARCHD)/sun
SUNOS_MPI_LIBS	= -L/Net/local/mpi/mpich/lib/sun4/ch_p4 -lpmpich -lmpich \
		  -L/Net/local/gnu/lib/gcc-lib/sparc-sun-sunos4.1.3_U1/2.6.3 \
		  -lgcc

SUNOS_DR_CSRC	= $(DR_CSRC) $(CH_CSRC)
SUNOS_DR_OBJ	= $(SUNOS_DR_CSRC:.c=.o)
SUNOS_DR_IPATH	= -I$(SUNOS_ARCHD)/include

SUNOS_DR_LIBS	= -L. -l$(LIB_TARGET) -L$(SUNOS_ARCHD)/lib $(NEM_LIBS) \
	          $(SUNOS_PAR_LPATH) $(PARMETIS_LIBS) $(SUNOS_MPI_LIBS) -lm 

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
CPLANT_PAR_LPATH   = -L$(PARMETISLIB_PATH)/Obj_cplant \
                  -L../$(PARMETISLIB_PATH)/Obj_cplant
CPLANT_DEFS	= $(GDEFINES)

CPLANT_CCOMP	= $(CPLANT_CC) $(CPLANT_INC_PATH) $(CPLANT_OPT_FLAGS)

# driver stuff
CPLANT_ARCHD	= $(ROOT_ARCHD)/alpha
CPLANT_MPI_LIBS	= -L/Net/mp/cplant/current/lib -lmpi

CPLANT_DR_CSRC	= $(DR_CSRC) $(CH_CSRC)
CPLANT_DR_OBJ	= $(CPLANT_DR_CSRC:.c=.o)
CPLANT_DR_IPATH	= -I$(CPLANT_ARCHD)/include

CPLANT_DR_LIBS	= -L. -l$(LIB_TARGET) -L$(CPLANT_ARCHD)/lib $(NEM_LIBS) \
	          $(CPLANT_PAR_LPATH) $(PARMETIS_LIBS) $(CPLANT_MPI_LIBS) -lm

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
#
############################################################################
#          UTILITY TARGETS
############################################################################
#

default:
	@echo " "
	@echo "usage: make <target>"
	@echo "where target is one of:"
	@for targ in ${SUPP_TARGS}; \
	do echo "  $$targ"; \
	done
	@echo " "

$(LINT_TARGET): $(LIB_SRC_C)
	@echo "Generating lint file for library..."
	$(LINT) $(LINT_FLAGS) $? -y -o $(LIB_TARGET)

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
	 cd ..; done; exit 0
