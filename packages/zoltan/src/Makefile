##############################################################################
# Zoltan Dynamic Load-Balancing Library for Parallel Applications            #
# Copyright (c) 2000, Sandia National Laboratories.                          #
# For more info, see the README file in the top-level Zoltan directory.      # 
##############################################################################
##############################################################################
# CVS File Information
#    $RCSfile$
#    $Author$
#    $Date$
#    $Revision$
##############################################################################

##############################################################################
# Top-level Makefile for Zoltan library and Zoltan driver zdrive.
#
# To compile the library, type
#    gmake ZOLTAN_ARCH=<targetmachine> zoltan
# where
#    targetmachine = solaris, cplant, smos, or tflop.
#
# To compile zdrive, type
#    gmake <optional flags> zdrive
#
# The optional flags are:
#    YES_FORTRAN=1   Use to compile and link a version of Zoltan that can
#                    be called from Fortran.
# Example:
# To compile a version of Zoltan for tflop including the F90 interface, type
#    gmake YES_FORTRAN=1 ZOLTAN_ARCH=tflop zoltan
##############################################################################

TARGET=zdrive
LIB_NAME=libzoltan.a
LIB_TARGET=zoltan

##############################################################################
# Process architecture flag
##############################################################################

ifndef ZOLTAN_ARCH
ZOLTAN_ARCH	= generic
endif

ifndef ZOLTAN_CONFIG
ZOLTAN_CONFIG = Utilities/Config/Config.$(ZOLTAN_ARCH)
endif

include $(ZOLTAN_CONFIG)

ZOLTAN_OBJ_DIR = Obj_$(ZOLTAN_ARCH)

ifndef ZOLTAN_LIB_DIR
ZOLTAN_LIB_DIR = $(ZOLTAN_OBJ_DIR)
endif

##############################################################################
# Process input flags
##############################################################################
GDEFINES	= -DDEBUG

# If YES_FORTRAN, compile the Fortran interface into the library.
ifdef YES_FORTRAN
  INC_FORTRAN	= inc_fortran
endif

##############################################################################
# Process third-party library information.
##############################################################################

THIRD_PARTY_CSRC =
THIRD_PARTY_VPATH =
THIRD_PARTY_LIBS =
THIRD_PARTY_LIBPATH = 
THIRD_PARTY_INCPATH =

# If PARMETIS_LIBPATH is defined, compile the ParMETIS interface and 
# link ParMETIS into zdrive.
ifdef PARMETIS_LIBPATH
  GDEFINES     += -DLB_PARMETIS
  THIRD_PARTY_LIBS += -lparmetis -lmetis
  THIRD_PARTY_LIBPATH += $(PARMETIS_LIBPATH)
endif

# If JOSTLE_LIBPATH is defined, compile the JOSTLE interface and 
# link JOSTLE into zdrive.
ifdef JOSTLE_LIBPATH
  GDEFINES     += -DLB_JOSTLE
  THIRD_PARTY_LIBS += -ljostle
  THIRD_PARTY_INCPATH += $(JOSTLE_INCPATH)
  THIRD_PARTY_LIBPATH += $(JOSTLE_LIBPATH)
  THIRD_PARTY_CSRC += jmpilib.c
  THIRD_PARTY_VPATH += $(JOSTLE_SRCPATH)
endif

# If NEMESIS_LIBPATH is defined, link Nemesis, Exodus and Netcdf into zdrive.
ifdef NEMESIS_LIBPATH
  GDEFINES += -DLB_NEMESIS
  THIRD_PARTY_LIBS += -lnemIc -lexoIIc -lnetcdf
  THIRD_PARTY_LIBPATH += $(NEMESIS_LIBPATH)
  THIRD_PARTY_INCPATH += $(NEMESIS_INCPATH)
endif


############################################################################
#                  OBJECT DIRECTORIES
############################################################################

ALL_OBJ_DIRS := $(wildcard Obj_*)

############################################################################
#                  SOURCE CODE
############################################################################
#
#   Paths to submodules
#
SUB_VPATH	= ../all ../lb ../par ../rcb ../oct \
		  ../params ../driver ../ch ../parmetis ../timer \
                  ../ha ../reftree \
		  ../Utilities/Communication ../Utilities/Memory
SUBINCPATH	= -I../all -I../lb -I../par -I../rcb -I../comm -I../oct \
		  -I../params -I../driver -I../ch -I../parmetis -I../timer \
                  -I../ha -I../reftree \
		  -I../Utilities/Communication -I../Utilities/Memory
SUBDIRS		= driver ch all lb par rcb comm oct parmetis params timer ha \
		  reftree Utilities/Communication Utilities/Memory


ALL_CSRC	= all_allo.c free_structure.c
ALL_INC		= all_const.h all_allo_const.h

LB_CSRC		= lb_balance.c lb_eval.c lb_init.c lb_migrate.c lb_set_fn.c \
		  lb_set_method.c lb_struct.c lb_util.c 
LB_INC		= lb_const.h lb.h lbi_const.h lb_util_const.h

RCB_CSRC	= rcb.c rcb_util.c box_assign.c point_assign.c \
                  create_proc_list.c irb.c irb_util.c \
                  inertial1d.c inertial2d.c inertial3d.c
RCB_INC		= rcb_const.h create_proc_list.h irb_const.h

OCT_CSRC	= msg.c octant.c oct_util.c octupdate.c dfs.c costs.c migoct.c \
		  output.c migreg.c hilbert.c
OCT_INC		= msg_const.h octant_const.h oct_util_const.h \
		  octupdate_const.h \
		  octupdate.h dfs_const.h dfs.h migoct_const.h \
		  migreg.h hilbert_const.h

PARAMS_CSRC	= set_param.c assign_param_vals.c check_param.c print_params.c \
		  key_params.c free_params.c bind_param.c 
PARAMS_INC	= params_const.h

PAR_CSRC	= par_util.c
PAR_INC		= par_const.h par.h

PARMETIS_CSRC	= parmetis_jostle.c verify_graph.c scatter_graph.c
PARMETIS_INC	= parmetis_jostle_const.h

TIMER_CSRC	= timer.c
TIMER_INC	= timer_const.h

HA_CSRC         = divide_machine.c get_processor_name.c \
                  build_machine_desc.c
HA_INC          = ha_const.h

REFTREE_CSRC    = reftree_build.c reftree_part.c reftree_hash.c
REFTREE_INC     = reftree_const.h

# source code and defines for zdrive
DR_CSRC		= dr_main.c dr_err.c dr_input.c dr_par_util.c \
		  dr_util.c dr_exoII_io.c dr_elem.c dr_loadbal.c \
		  dr_migrate.c dr_output.c dr_chaco_io.c dr_elem_util.c \
		  dr_maps.c dr_eval.c
DR_INC		= dr_const.h dr_input_const.h dr_util_const.h \
		  dr_err_const.h dr_par_util_const.h dr_elem_const.h \
		  dr_elem_util_const.h dr_maps_const.h dr_eval_const.h

CH_CSRC		= ch_dist_graph.c ch_input_geom.c ch_input_graph.c \
                  ch_input_read.c ch_init_dist.c

CH_INC		= ch_input_const.h ch_init_dist_const.h


# object files for the Fortran interface

ifdef YES_FORTRAN
FOBJ_FILES	= lb_user_const.o ztypes.o fwrap.o cwrap.o
endif



##############################################################################
#  Architecture-independent definitions
##############################################################################

ZOLTAN_INC_PATH	= $(SUBINCPATH) $(INCLUDE_PATH) $(THIRD_PARTY_INCPATH)
ZOLTAN_CSRC    += $(LB_CSRC) $(ALL_CSRC) $(RCB_CSRC) $(PAR_CSRC) \
		  $(OCT_CSRC) $(PARAMS_CSRC) $(PARMETIS_CSRC) \
		  $(TIMER_CSRC) $(HA_CSRC) $(REFTREE_CSRC) $(THIRD_PARTY_CSRC)
ZOLTAN_VPATH	= $(SUB_VPATH) $(THIRD_PARTY_VPATH)
ZOLTAN_OBJ	= $(ZOLTAN_CSRC:.c=.o) 

ZDRIVE_CSRC 	= $(DR_CSRC) $(CH_CSRC)
ZDRIVE_OBJ	= $(ZDRIVE_CSRC:.c=.o)
ZDRIVE_LIBS	= -L. -l$(LIB_TARGET) \
		  -L../Utilities/$(ZOLTAN_OBJ_DIR) -lzoltan_comm -lzoltan_mem \
	          $(THIRD_PARTY_LIBPATH) $(THIRD_PARTY_LIBS) \
		  $(MPI_LIBPATH) $(MPI_LIB) -lm

CCOMP		= $(CC) $(ZOLTAN_INC_PATH) $(CFLAGS) $(F90CFLAGS) 
F90COMP		= $(F90) $(CFLAGS) $(FFLAGS)

############################################################################
#                  LINT DEFINITIONS KDDKDD
############################################################################
#
# lint defines
#
LINT		= lint
LINT_INC_PATH	= -Idriver -Ich -Iall -Ilb -Ipar -Ircb -Iparams\
		  -Icomm -Ioct -Iparmetis -Itimer -Iha -Ireftree\
		  -I$(ZOLTAN_INC_PATH)
#		  -I/Net/local/lang/SC3.0.1/include \
#		  -I/Net/local/lang/SC3.0.1/ansi_include \
#		  -I/Net/local/mpi/include \
#		  -I/Net/local/mpi/mpich/include \
#                  -I/Net/local/mpi/build/solaris/ch_p4/include \
#		  -I/Net/proj/zoltan/arch/sun/include
LINT_FLAGS	= $(LINT_INC_PATH) $(GDEFINES)
LINT_TARGET	= llib-l$(LIB_TARGET).ln
VPATH		= $(SUBDIRS) $(THIRD_PARTY_VPATH)

###########################################################################
#      HELP HELP HELP HELP HELP  
###########################################################################
help:
	@echo ' '
	@echo 'To build the Zoltan library only, type '
	@echo '   gmake ZOLTAN_ARCH=desired_arch zoltan'
	@echo 'where desired_arch is the architecture on which Zoltan'
	@echo 'will run.  The file Makefiles/Makefile.desired_arch '
	@echo 'should exist.  See Makefiles/Makefile.generic for '
	@echo 'instructions on creating Makefiles/Makefile.desired_arch'
	@echo 'and Makefiles/Makefile.* for examples.'
	@echo ' '
	@echo 'To build the the test driver zdrive, type'
	@echo '   gmake ZOLTAN_ARCH=desired_arch zdrive'
	@echo 'where desired_arch is defined as above.'
	@echo ' '
	@echo 'To include the F90 interface in Zoltan, use the gmake
	@echo 'command line argument'
	@echo '      YES_FORTRAN=1  (to compile fortran interface)'
	@echo ' '
	@echo 'Example:'
	@echo '   gmake ZOLTAN_ARCH=tflop YES_FORTRAN=1 zoltan'
	@echo 'builds the tflop version of the Zoltan library'
	@echo 'including the F90 interface.'


###########################################################################
#      TARGETS FOR DIFFERENT MACHINES
###########################################################################


zoltan:
	@cd Utilities; \
	$(MAKE) libzoltan_mem.a; \
	$(MAKE) libzoltan_comm.a  
	@if [ ! -d $(ZOLTAN_OBJ_DIR) ]; 	\
	then mkdir $(ZOLTAN_OBJ_DIR); 		\
	fi
	@rm -f $(ZOLTAN_OBJ_DIR)/Makefile
	@cp Makefile_sub $(ZOLTAN_OBJ_DIR)/Makefile
	@cd $(ZOLTAN_OBJ_DIR); 			\
	$(MAKE) $(LIB_NAME) 			\
	"LIB_NAME=$(LIB_NAME)"			\
	"CC=$(CCOMP)"				\
	"VPATH=$(ZOLTAN_VPATH)"			\
	"OBJ_FILES=$(ZOLTAN_OBJ)" 		\
	"DEFS=$(GDEFINES) $(DEFS)" 		\
	"INC_FORTRAN=$(INC_FORTRAN)" 		\
	"F90=$(F90COMP)" 			\
	"LOCAL_F90=$(LOCAL_F90)" 		\
	"FOBJ_FILES=$(FOBJ_FILES)" 		\
	"F90_MODULE_PREFIX=$(F90_MODULE_PREFIX)"\
	"SPPR_HEAD=$(SPPR_HEAD)" 		\
	"AR=$(AR)" 				\
	"RANLIB=$(RANLIB)" 			\
	"ZOLTAN_OBJ_DIR=$(ZOLTAN_OBJ_DIR)" 
	@if [ $(ZOLTAN_OBJ_DIR) != $(ZOLTAN_LIB_DIR) ]; \
	then  						\
          if [ ! -d $(ZOLTAN_LIB_DIR) ];		\
	  then mkdir $(ZOLTAN_LIB_DIR);			\
	  fi;						\
	  mv -f $(ZOLTAN_OBJ_DIR)/$(LIB_NAME) $(ZOLTAN_LIB_DIR); \
	fi

zdrive:
	@cd Utilities; \
	$(MAKE) libzoltan_mem.a; \
	$(MAKE) libzoltan_comm.a  
	@if [ ! -d $(ZOLTAN_OBJ_DIR) ]; 	\
	then mkdir $(ZOLTAN_OBJ_DIR); 		\
	fi
	@rm -f $(ZOLTAN_OBJ_DIR)/Makefile
	@cp Makefile_sub $(ZOLTAN_OBJ_DIR)/Makefile
	@cd $(ZOLTAN_OBJ_DIR); 			\
	$(MAKE) $(TARGET)			\
	"TARGET=$(TARGET)"			\
	"LIB_NAME=$(LIB_NAME)"			\
	"CC=$(CCOMP)" 				\
	"VPATH=$(ZOLTAN_VPATH)"			\
	"OBJ_FILES=$(ZOLTAN_OBJ)"		\
	"DEFS=$(GDEFINES) $(DEFS)" 		\
	"INC_FORTRAN=$(INC_FORTRAN)" 		\
	"F90=$(F90COMP)"			\
	"LOCAL_F90=$(LOCAL_F90)" 		\
	"FOBJ_FILES=$(FOBJ_FILES)" 		\
	"F90_MODULE_PREFIX=$(F90_MODULE_PREFIX)"\
	"SPPR_HEAD=$(SPPR_HEAD)" 		\
	"AR=$(AR)" 				\
	"RANLIB=$(RANLIB)" 			\
	"MAIN_OBJ=$(ZDRIVE_OBJ)" 		\
	"LNK_LIBS=$(ZDRIVE_LIBS)" 		\
	"ZOLTAN_OBJ_DIR=$(ZOLTAN_OBJ_DIR)" 

zfdrive: 
	$(MAKE) YES_FORTRAN=1 ZOLTAN_ARCH=$(ZOLTAN_ARCH) zoltan
	@cd fdriver; 				\
	$(MAKE) zfdrive 			\
	"CC=$(CCOMP)"				\
	"ZOLTAN_INC_PATH=$(ZOLTAN_INC_PATH)"	\
	"F90=$(F90COMP)" 			\
	"F90_MODULE_PREFIX=$(F90_MODULE_PREFIX)"\
	"FARG=$(FARG)" 				\
        "LNK_LIBS=$(ZDRIVE_LIBS)"		\
	"ZOLTAN_ARCH=$(ZOLTAN_ARCH)"		\
	"ZOLTAN_OBJ_DIR=$(ZOLTAN_OBJ_DIR)" 

zscript:
	@echo "Building zoltan_script..."
	@if [ -f zoltan_script ]; \
	then /bin/rm zoltan_script; \
	fi
	@echo "#!/bin/sh" > zoltan_script
	@echo "unalias rm" >> zoltan_script
	@echo DEFS='"$$1"' >> zoltan_script
	@echo CC='"$$2"' >> zoltan_script
	@echo CFLAGS='"$$3"' >> zoltan_script
	@echo F90CFLAGS='"$$4"' >> zoltan_script
	@echo INCLUDE_PATH='"$$5"' >> zoltan_script
	@echo PARMETIS_FLAG='"$${6:+-DLB_PARMETIS}"' >> zoltan_script
	@echo JOSTLE_FLAG='"$${7:+-DLB_JOSTLE}"' >> zoltan_script
	@echo JOSTLE_INCPATH='"$$8"' >> zoltan_script
	@echo AR='"$$9"' >> zoltan_script
	@echo shift 9 >> zoltan_script
	@echo RANLIB='"$$1"' >> zoltan_script
	@echo ZOLTAN_ARCH='"$${2:-generic}"' >> zoltan_script
	@echo ZOLTAN_LIB_DIR='"$$3"' >> zoltan_script
	@echo CCOMP='"$$CC $$CFLAGS $$F90CFLAGS $$INCLUDE_PATH $$JOSTLE_INCPATH"' >> zoltan_script
	@echo GDEFINES='"-DDEBUG $$PARMETIS_FLAG $$JOSTLE_FLAG"' >> zoltan_script
	@echo ZOLTAN_OBJ_DIR='"Obj_$$ZOLTAN_ARCH"' >> zoltan_script
	@echo 'if [ ! -d $$ZOLTAN_OBJ_DIR ]; ' >> zoltan_script
	@echo 'then mkdir $$ZOLTAN_OBJ_DIR ;' >> zoltan_script
	@echo 'fi' >> zoltan_script
	@echo 'rm -f $$ZOLTAN_OBJ_DIR/*' >> zoltan_script
	@echo 'cd $$ZOLTAN_OBJ_DIR' >> zoltan_script
	@echo "echo 'Zoltan:  Compiling Zoltan/all'" >> zoltan_script
	@for iii in $(ALL_CSRC)  ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../all/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/lb'" >> zoltan_script
	@for iii in $(LB_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../lb/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/rcb'" >> zoltan_script
	@for iii in $(RCB_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../rcb/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/oct'" >> zoltan_script
	@for iii in $(OCT_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../oct/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/params'" >> zoltan_script
	@for iii in $(PARAMS_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../params/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/comm'" >> zoltan_script
	@for iii in $(COMM_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../comm/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/par'" >> zoltan_script
	@for iii in $(PAR_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../par/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/parmetis'" >> zoltan_script
	@for iii in $(PARMETIS_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../parmetis/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/timer'" >> zoltan_script
	@for iii in $(TIMER_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../timer/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/ha'" >> zoltan_script
	@for iii in $(HA_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../ha/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/reftree'" >> zoltan_script
	@for iii in $(REFTREE_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../reftree/$$iii >> zoltan_script; \
	 done;
	@echo "echo Zoltan:  Linking $(LIB_NAME)" >> zoltan_script
	@echo '$$AR' $(LIB_NAME) $(OBJ_FILES) $(FOBJ_FILES) >> zoltan_script
	@echo '$$RANLIB' $(LIB_NAME) >> zoltan_script
	@echo 'if [ $$ZOLTAN_LIB_DIR ] ; ' >> zoltan_script
	@echo 'then mv -f $(LIB_NAME) $$ZOLTAN_LIB_DIR ; ' >> zoltan_script
	@echo 'fi' >> zoltan_script
	@echo "cd .." >> zoltan_script
	@chmod +x zoltan_script

#
############################################################################
#          UTILITY TARGETS
############################################################################
#

$(LINT_TARGET): $(ZOLTAN_CSRC)
	@echo "Generating lint file for library..."
	$(LINT) $(LINT_FLAGS)  $? -o $(LIB_TARGET)

lint: $(ZDRIVE_CSRC) $(LINT_TARGET)
	@echo "running lint..."
	$(LINT) $(LINT_FLAGS) "VPATH=$(LINT_VPATH) $?

printout:
	cat $(MAIN) $(ZOLTAN_CSRC) > source.c
	tps -c source.c | lpr
	-rm source.c
	tps  $(INC) | lpr

clean:
	-for a in $(SUBDIRS); do \
	 cd $$a; \
         echo $$a; \
	 $(RM) -f *.o *.d *.d@ $(TARGET) $(TARGET).pure $(LIB_NAME); \
	 cd ..; \
	 done;
	-for a in $(ALL_OBJ_DIRS); do \
         echo $$a; \
	 cd fort; $(MAKE) "ZOLTAN_OBJ_DIR=$$a" clean; cd ..; \
	 cd fdriver; $(MAKE) "ZOLTAN_OBJ_DIR=$$a" clean; cd ..; \
	 cd $$a; $(RM) -f *.o *.d *.d@ $(TARGET) $(TARGET).pure $(LIB_NAME); cd ..; \
	 done;
	 $(RM) $(LINT_TARGET)
	 exit 0
