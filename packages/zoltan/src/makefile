##############################################################################
# Zoltan Library for Parallel Applications                                   #
# Copyright (c) 2000,2001,2002, Sandia National Laboratories.                #
# For more info, see the README file in the top-level Zoltan directory.      # 
##############################################################################
##############################################################################
# CVS File Information
#    $RCSfile$
#    $Author$
#    $Date$
#    Revision: 1.107 $
##############################################################################

##############################################################################
# Top-level makefile for Zoltan library, test drivers, and graphical tools
#
# To compile the library, type
#    gmake <optional flags> ZOLTAN_ARCH=<targetmachine> zoltan
# where
#    targetmachine = the suffix of one of the configuration files Config.*
#                    in Zoltan/Utilities/Config.
#
# The optional flags are:
#    YES_FORTRAN=1   Use to compile and link a version of Zoltan that can
#                    be called from Fortran.
#
# To compile the C, Fortran, or C++ test driver:
#    gmake ZOLTAN_ARCH=<> zdrive
#    gmake ZOLTAN_ARCH=<> zfdrive
#    gmake ZOLTAN_ARCH=<> zCPPdrive
#
# Two graphical tools for visualizing geometric partitions: 
#
# gmake ZOLTAN_ARCH=<> vtk_view (needs VTK from www.vtk.org and OpenGL)
# gmake ZOLTAN_ARCH=<> vtk_write (needs VTK and Mesa GL from www.mesa.org)
#
# Example:
# To compile a version of Zoltan for tflop including the F90 interface, type
#    gmake YES_FORTRAN=1 ZOLTAN_ARCH=tflop zoltan
##############################################################################

ZTARGET=zdrive 
LIB_NAME=libzoltan.a     # if change LIB_NAME, must change ZOLTAN_LIB 
                         # in Utilities/makefile
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

ifndef ZOLTAN_OBJ_DIR
ZOLTAN_OBJ_DIR = Obj_$(ZOLTAN_ARCH)
endif

ifndef ZOLTAN_LIB_DIR
ZOLTAN_LIB_DIR = $(ZOLTAN_OBJ_DIR)
endif

ifeq ($(ZOLTAN_LIB_DIR), $(ZOLTAN_OBJ_DIR))
ZOLTAN_SUBLIB_DIR = ../$(ZOLTAN_OBJ_DIR)
else
ZOLTAN_SUBLIB_DIR = $(ZOLTAN_LIB_DIR)
endif

# If ZOLTAN_OLD_INTERFACE is defined, builds zdrive executable using the
# old LB_* interface.  Source code from driver_old and fdriver_old is used.
ifdef ZOLTAN_OLD_INTERFACE
ZOLTAN_OLD = _old
else
ZOLTAN_OLD = 
endif

##############################################################################
# Process input flags
##############################################################################
# Default values for inclusion of methods.
ZOLTAN_OCT=1
ZOLTAN_HIER=1

GDEFINES	= -DTRILINOS_NO_CONFIG_H

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
  GDEFINES     += -DZOLTAN_PARMETIS
  THIRD_PARTY_LIBS += -lparmetis -lmetis
  THIRD_PARTY_LIBPATH += $(PARMETIS_LIBPATH)
  THIRD_PARTY_INCPATH += $(PARMETIS_INCPATH)
endif

# If SCOTCH_LIBPATH is defined, compile the SCOTCH interface and 
# link SCOTCH into zdrive.
ifdef SCOTCH_LIBPATH
  GDEFINES     += -DZOLTAN_SCOTCH
  THIRD_PARTY_LIBS +=  -lptscotch -lscotch -lscotcherr
  THIRD_PARTY_LIBPATH += $(SCOTCH_LIBPATH)
  THIRD_PARTY_INCPATH += $(SCOTCH_INCPATH)
endif	


# If PATOH_LIBPATH is defined, compile the PaToH interface and 
# link PaToH into zdrive.
# KDD -- does not work with zscript yet.
ifdef PATOH_LIBPATH
  GDEFINES     += -DZOLTAN_PATOH
  THIRD_PARTY_LIBS += -lpatoh 
  THIRD_PARTY_LIBPATH += $(PATOH_LIBPATH)
  THIRD_PARTY_INCPATH += $(PATOH_INCPATH)
endif

# If PARKWAY_LIBPATH is defined, compile the ParKway interface and 
# link ParKway into zdrive.
# KDD -- does not work with zscript yet.
ifdef PARKWAY_LIBPATH
  GDEFINES     += -DZOLTAN_PARKWAY
  THIRD_PARTY_LIBS += -lparkway -lpmpich++ -lstdc++
  THIRD_PARTY_LIBPATH += $(PARKWAY_LIBPATH)
  THIRD_PARTY_INCPATH += $(PARKWAY_INCPATH)
endif


# If NEMESIS_LIBPATH is defined, link Nemesis, Exodus and Netcdf into zdrive.
ifdef NEMESIS_LIBPATH
  GDEFINES += -DZOLTAN_NEMESIS
  THIRD_PARTY_LIBS += -lnemIc -lexoIIv2c -lnetcdf
  THIRD_PARTY_LIBPATH += $(NEMESIS_LIBPATH)
  THIRD_PARTY_INCPATH += $(NEMESIS_INCPATH)
endif

# If DRUM_LIBPATH is defined, set include path and link DRUM into zdrive
ifdef DRUM_LIBPATH
  GDEFINES += -DZOLTAN_DRUM
  THIRD_PARTY_LIBS += $(DRUM_LIBNAME)
  THIRD_PARTY_LIBPATH += $(DRUM_LIBPATH)
  THIRD_PARTY_INCPATH += $(DRUM_INCPATH)
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
SUB_VPATH	= ../include ../Utilities/shared \
		  ../zz ../all ../lb ../par ../rcb ../coloring ../oct ../hsfc \
		  ../phg  ../simple ../util \
		  ../params ../driver$(ZOLTAN_OLD) ../ch ../parmetis ../timer \
                  ../order ../ha ../reftree ../hier  \
		  ../Utilities/Communication ../Utilities/Memory \
                  ../Utilities/DDirectory ../Utilities/Timer \
		  ../postprocessing/ordering
SUBINCPATH	= -I../include -I../Utilities/shared \
		  -I../zz -I../all -I../lb -I../par -I../rcb \
		  -I../coloring -I../simple -I../oct  \
		  -I../order -I../hsfc -I../phg  -I../util \
		  -I../params -I../driver$(ZOLTAN_OLD) -I../ch \
		  -I../parmetis -I../timer  \
                  -I../ha -I../hier -I../reftree \
		  -I../Utilities/Communication -I../Utilities/Memory \
                  -I../Utilities/DDirectory -I../Utilities/Timer
SUBDIRS		= driver$(ZOLTAN_OLD) ch zz all lb order par rcb \
		  coloring oct phg simple util \
                  hsfc parmetis params hier  \
		  timer ha  reftree include postprocessing/ordering
UTIL_SUBDIRS	= Utilities/shared Utilities/Communication \
		  Utilities/Memory Utilities/DDirectory Utilities/Timer

ZZ_CSRC		= zz_coord.c zz_obj_list.c zz_struct.c \
		  zz_init.c zz_set_fn.c zz_util.c zz_gen_files.c \
                  zz_hash.c zz_heap.c zz_sort.c zz_rand.c
ZZ_INC		= zz_const.h zz_id_const.h zz_util_const.h zz_rand.h zz_sort.h

ALL_CSRC	= all_allo.c 
ALL_INC		= all_const.h all_allo_const.h

LB_CSRC		= lb_balance.c lb_eval.c lb_free.c lb_copy.c lb_init.c \
		  lb_invert.c lb_migrate.c \
		  lb_set_fn.c lb_set_method.c  \
		  lb_point_assign.c lb_box_assign.c lb_part2proc.c \
                  lb_set_part_sizes.c lb_remap.c
LB_INC		= lb_const.h lb_init_const.h

ORDER_CSRC	= order.c order_struct.c order_tools.c perm.c
ORDER_INC	= order_const.h 

SIMPLE_CSRC	= block.c random.c
SIMPLE_INC	= 

RCB_CSRC	= rcb.c rcb_util.c rcb_box.c box_assign.c point_assign.c \
                  create_proc_list.c rib.c rib_util.c \
                  inertial1d.c inertial2d.c inertial3d.c shared.c
RCB_INC		= rcb_const.h create_proc_list.h rib_const.h shared_const.h

COLOR_CSRC      = coloring.c g2l_hash.c color_test.c
COLOR_INC       = coloring.h g2l_hash.h color_test.c

ifdef ZOLTAN_OCT
OCT_CSRC	= msg.c octant.c oct_util.c octupdate.c dfs.c costs.c migoct.c \
		  output.c migreg.c migtags.c octree.c rootlist.c oct_plot.c
OCT_INC		= msg_const.h octant_const.h oct_util_const.h \
		  octupdate_const.h \
		  octupdate.h dfs_const.h dfs.h migoct_const.h \
		  migreg.h hsfc_hilbert_const.h migtags_const.h octree_const.h \
		  rootlist_const.h SFC.h
GDEFINES       += -DZOLTAN_OCT
else
OCT_CSRC	=
OCT_INC		=
endif

HSFC_CSRC 	= hsfc_hilbert.c hsfc.c hsfc_box_assign.c hsfc_point_assign.c
HSFC_INC 	= hsfc_hilbert_const.h hsfc_const.h hsfc.h

PHG_CSRC 	= phg.c phg_hypergraph.c phg_build.c phg_build_calls.c \
		  phg_coarse.c phg_distrib.c phg_gather.c \
		  phg_match.c phg_plot.c \
		  phg_refinement.c phg_serialpartition.c \
		  phg_Vcycle.c phg_rdivide.c phg_util.c \
		  phg_scale.c phg_comm.c phg_order.c \
                  phg_parkway.c phg_patoh.c
PHG_INC 	= phg.h phg_util.h phg_hypergraph.h 

PARAMS_CSRC	= set_param.c assign_param_vals.c check_param.c print_params.c \
		  key_params.c free_params.c bind_param.c
PARAMS_INC	= params_const.h

PAR_CSRC	= par_average.c par_bisect.c par_median.c par_median_randomized.c \
		  par_sync.c par_stats.c \
                  par_tflops_special.c
PAR_INC		= par_const.h par_average_const.h par_median_const.h \
		  par_bisect_const.h \
                  par_tflops_special_const.c

PARMETIS_CSRC	= third_library.c preprocessing.c postprocessing.c build_graph.c verify_graph.c \
                  scatter_graph.c
PARMETIS_INC	= third_library_const.h
ifdef PARMETIS_LIBPATH
PARMETIS_CSRC  += parmetis_interface.c
endif
ifdef SCOTCH_LIBPATH
PARMETIS_CSRC  += scotch_interface.c
endif


TIMER_CSRC	= timer_params.c
TIMER_INC	= timer_const.h

HA_CSRC         = divide_machine.c get_processor_name.c \
                  build_machine_desc.c ha_drum.c
HA_INC          = ha_const.h ha_drum.h

REFTREE_CSRC    = reftree_build.c reftree_part.c reftree_hash.c \
                  reftree_coarse_path.c
REFTREE_INC     = reftree_const.h

ifdef ZOLTAN_HIER
HIER_CSRC       = hier.c hier_free_struct.c
HIER_INC        = hier.h
GDEFINES       += -DZOLTAN_HIER
else
HIER_CSRC       = 
HIER_INC        = 
endif


# source code and defines for zdrive and zCPPdrive
DR_BOTH_SRC	= dr_err.c dr_input.c dr_par_util.c \
		  dr_util.c dr_elem.c dr_setfixed.c \
		  dr_output.c dr_chaco_io.c dr_random_io.c \
                  dr_elem_util.c dr_eval.c dr_gnuplot.c \
		  dr_hg_io.c dr_hg_readfile.c dr_mm_readfile.c dr_mmio.c

DR_C_SRC	= dr_main.c dr_exoII_io.c dr_loadbal.c \
		  dr_migrate.c dr_maps.c dr_dd.c dr_param_file.c order_test.c

DR_CPP_SRC	= dr_mainCPP.cpp dr_exoII_ioCPP.cpp dr_loadbalCPP.cpp \
		  dr_migrateCPP.cpp dr_mapsCPP.cpp dr_ddCPP.cpp \
		  dr_param_fileCPP.cpp

DR_BOTH_INC	= dr_const.h dr_input_const.h dr_util_const.h \
		  dr_err_const.h dr_par_util_const.h dr_elem_const.h \
		  dr_elem_util_const.h dr_maps_const.h dr_eval_const.h \
                  zoltan_hg_readfile_const.h dr_mmio.h dr_dd.h

ifdef GZIP_LIBPATH
GDEFINES       += -DZOLTAN_GZIP
DR_BOTH_SRC    += dr_compress.c
DR_BOTH_INC    += dr_compress_const.h
THIRD_PARTY_LIBS += -lz 
THIRD_PARTY_LIBPATH += $(GZIP_LIBPATH)
THIRD_PARTY_INCPATH += $(GZIP_INCPATH)
endif

DR_CPP_INC	= zoltan_cpp.h zoltan_comm_cpp.h zoltan_dd_cpp.h \
		  dr_param_fileCPP.h

CH_CSRC		= ch_dist_graph.c ch_input_assign.c ch_input_geom.c \
		  ch_input_graph.c ch_input_read.c ch_init_dist.c

CH_INC		= ch_input_const.h ch_init_dist_const.h

# source code and defines for vtk_view

VIEW_SRC        = vtk_view.cpp dr_input.c dr_err.c

VIEW_INC        = dr_input_const.h


# object files for the Fortran interface

ifdef YES_FORTRAN
FOBJ_FILES	= zoltan_user_data.o lb_user_const.o ztypes.o fwrap.o cwrap.o
endif



##############################################################################
#  Architecture-independent definitions
##############################################################################

ZOLTAN_INC_PATH	= $(SUBINCPATH) $(INCLUDE_PATH) $(MPI_INCPATH) \
		  $(THIRD_PARTY_INCPATH)
ZOLTAN_CSRC    += $(ZZ_CSRC) $(LB_CSRC) $(ALL_CSRC) $(RCB_CSRC) $(PAR_CSRC) \
		  $(COLOR_CSRC) $(SIMPLE_CSRC) \
		  $(ORDER_CSRC) $(OCT_CSRC) $(HSFC_CSRC) \
		  $(PARAMS_CSRC) $(PARMETIS_CSRC) $(PHG_CSRC) $(HG_CSRC) \
		  $(TIMER_CSRC) $(HA_CSRC) $(REFTREE_CSRC) \
		  $(HIER_CSRC) $(THIRD_PARTY_CSRC) $(SIMPLE_CSRC) 
ZOLTAN_VPATH	= $(SUB_VPATH) $(THIRD_PARTY_VPATH)
ZOLTAN_OBJ	= $(ZOLTAN_CSRC:.c=.o)

ZDRIVE_CSRC 	= $(DR_C_SRC) $(DR_BOTH_SRC) $(CH_CSRC)
ZDRIVE_OBJ	= $(ZDRIVE_CSRC:.c=.o)
ZDRIVE_LIBS	= -L. -l$(LIB_TARGET) \
	          $(THIRD_PARTY_LIBPATH) $(THIRD_PARTY_LIBS) \
		  $(MPI_LIBPATH) $(MPI_LIB) -lm 
#		  -lzoltan_comm -lzoltan_mem -lzoltan_dd \

ZDRIVE_DEPS	= libzoltan_comm.a libzoltan_mem.a libzoltan_dd.a

ZCPPDRIVE_CSRC 	= $(DR_CPP_SRC) $(DR_BOTH_SRC) $(CH_CSRC)
ZCPPDRIVE_OBJ	= $(DR_BOTH_SRC:.c=.O) $(CH_CSRC:.c=.O) $(DR_CPP_SRC:.cpp=.o)
ZCPPDRIVE_LIBS	= -L. -l$(LIB_TARGET) \
	          $(THIRD_PARTY_LIBPATH) $(THIRD_PARTY_LIBS) \
		  $(MPI_LIBPATH) $(MPI_LIB) -lm 
ZCPPDRIVE_DEPS	= $(ZDRIVE_DEPS)


#Tool for computing ordering quality
ifdef SCOTCH_LIBPATH
OPP_CSRC	= mm_input.c order_eval.c
OPP_INC		= mm_input.h order_eval.h
OPP_OBJ		= $(OPP_CSRC:.c=.o)
OPP_LIBS	= $(SCOTCH_LIBPATH) -lscotch -lscotcherrcom -lcommon -lrt
OPP_INCLUDEPATH = $(SCOTCH_INCPATH)
endif


VIEW_OBJ        = vtk_view.o dr_input.o  dr_err.o
WRITE_OBJ        = vtk_write.o dr_input.o  dr_err.o

VIEW_LIBS       = -L. $(VTK_LIBPATH) $(MPI_LIBPATH) $(X_LIBPATH) $(GL_LIBPATH) \
   -lvtkParallel -lvtkHybrid -lvtkIO -lvtkRendering -lvtkGraphics -lvtkImaging \
   -lvtkFiltering -lvtkCommon \
   $(GL_LIBS) $(X_LIBS)  \
   -lvtkftgl -lvtkfreetype \
   -lvtkpng -lvtktiff -lvtkzlib -lvtkjpeg -lvtkexpat \
   -lpthread -lvtkexoIIc -lvtkNetCDF $(MPI_LIB) -ldl

OFFSCREEN_VIEW_LIBS = \
   -L. $(VTK_OFFSCREEN_LIBPATH) $(X_LIBPATH) $(MPI_LIBPATH) $(MESA_LIBPATH) \
   -lvtkHybrid -lvtkParallel -lvtkRendering -lvtkIO -lvtkGraphics -lvtkImaging \
   $(MESA_LIBS) $(X_LIBS) -lvtkftgl -lvtkfreetype \
   -lvtkFiltering -lvtkCommon \
   -lvtkpng -lvtktiff -lvtkzlib -lvtkjpeg -lvtkexpat \
   -lpthread -lvtkexoIIc -lvtkNetCDF $(MPI_LIB) -ldl


CCOMP		= $(CC) $(ZOLTAN_INC_PATH) $(CFLAGS) $(F90CFLAGS)
F90COMP		= $(F90) $(CFLAGS) $(FFLAGS)

ifdef CPPC
CPPCOMP		= $(CPPC) $(CFLAGS) $(CPPFLAGS)
else
# See if C compiler will work with C++ files
CPPCOMP		= $(CC) $(CFLAGS) 
endif

############################################################################
#                  LINT DEFINITIONS KDDKDD
############################################################################
#
# lint defines
# MPI Paths must be explicitly stated; see Utilities/Config/Config.gcc
# for an example.
#
LINT		= lint
LINT_INC_PATH	= -Iinclude -Idriver$(ZOLTAN_OLD) -Ich -Izz -Iall -Ilb \
		  -Iorder -Ipar -Ircb -Icoloring -Iparams \
		  -Ioct -Ihsfc -Iparmetis -Itimer -Iha -Ireftree \
		  -Iphg -Ihier \
		  -IUtilities/shared -IUtilities/Memory \
		  -IUtilities/Communication -IUtilities/DDirectory \
		  -IUtilities/Timer \
		  -I$(ZOLTAN_INC_PATH)
LINT_FLAGS	= $(LINT_INC_PATH) $(GDEFINES)
VPATH		= $(SUBDIRS) $(UTIL_SUBDIRS) $(THIRD_PARTY_VPATH)

###########################################################################
#      HELP HELP HELP HELP HELP
###########################################################################
help:
	@echo ' '
	@echo 'To build the Zoltan library only, type '
	@echo '   gmake ZOLTAN_ARCH=desired_arch zoltan'
	@echo ' '
	@echo 'where desired_arch is the architecture on which Zoltan'
	@echo 'will run.  The file Utilities/Config/Config.desired_arch'
	@echo 'should exist and may need to be adapted to your system.'
	@echo 'See  Utilities/Config/Config.generic for '
	@echo 'instructions on creating Config.desired_arch.'
	@echo ' '
	@echo 'To build the C, Fortran or C++ test driver, type'
	@echo '   gmake ZOLTAN_ARCH=desired_arch zdrive'
	@echo '   gmake ZOLTAN_ARCH=desired_arch zfdrive'
	@echo '   gmake ZOLTAN_ARCH=desired_arch zCPPdrive'
	@echo 'where desired_arch is defined as above.'
	@echo ' '
	@echo 'To include the F90 interface in Zoltan, use the gmake'
	@echo 'command line argument'
	@echo '      YES_FORTRAN=1  (to compile fortran interface)'
	@echo ' '
	@echo 'Example:'
	@echo '   gmake ZOLTAN_ARCH=tflop YES_FORTRAN=1 zoltan'
	@echo 'builds the tflop version of the Zoltan library'
	@echo 'including the F90 interface.'
	@echo ' '
	@echo 'The graphical tools require VTK and OpenGL or Mesa.'
	@echo 'To compile the display tool:'
	@echo '   gmake ZOLTAN_ARCH=desired_arch vtk_view '
	@echo 'To compile the off-screen image writing tool:'
	@echo '   gmake ZOLTAN_ARCH=desired_arch vtk_write '
	@echo ' '
	@echo 'To build examples, go to the examples directory and'
	@echo '   gmake ZOLTAN_ARCH=desired_arch all'



###########################################################################
#      TARGETS FOR DIFFERENT MACHINES
###########################################################################


zoltan:
	@if [ ! -d $(ZOLTAN_LIB_DIR) ];		\
	then mkdir $(ZOLTAN_LIB_DIR);		\
	fi;
	@if [ ! -d $(ZOLTAN_OBJ_DIR) ]; 	\
	then mkdir $(ZOLTAN_OBJ_DIR); 		\
	fi
	@cd Utilities; \
	$(MAKE) \
	"ZOLTAN=$(LIB_TARGET)" \
	"ZOLTAN_LIB_DIR=$(ZOLTAN_SUBLIB_DIR)" \
	"ZOLTAN_LIB=$(LIB_NAME)"  $(LIB_TARGET)
	@rm -f $(ZOLTAN_OBJ_DIR)/makefile
	@cp makefile_sub $(ZOLTAN_OBJ_DIR)/makefile
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
	"F90_MODULE_PREFIX=$(F90_MODULE_PREFIX)" \
	"SPPR_HEAD=$(SPPR_HEAD)" 		\
	"AR=$(AR)" 				\
	"RANLIB=$(RANLIB)" 			\
	"ZOLTAN_CONFIG=../$(ZOLTAN_CONFIG)"     \
	"ZOLTAN_OBJ_DIR=$(ZOLTAN_OBJ_DIR)"
	@if [ $(ZOLTAN_OBJ_DIR) != $(ZOLTAN_LIB_DIR) ]; 		\
	then cp -f $(ZOLTAN_OBJ_DIR)/$(LIB_NAME) $(ZOLTAN_LIB_DIR); 	\
	fi


ifdef SCOTCH_LIBPATH
order_eval:
	@rm -f $(ZOLTAN_OBJ_DIR)/makefile
	@cp makefile_sub $(ZOLTAN_OBJ_DIR)/makefile
	@cd $(ZOLTAN_OBJ_DIR); 			\
	$(MAKE) order_eval                      \
	"ZTARGET=order_eval"			\
	"CC=$(CCOMP)" 				\
	"VPATH=$(ZOLTAN_VPATH)"			\
	"DEFS=$(OPP_INCLUDEPATH)"		\
	"OBJ_FILES=$(OPP_OBJ)"			\
	"MAIN_OBJ=$(OPP_OBJ)" 			\
	"LNK_LIBS=$(OPP_LIBS)" 			\
	"ZOLTAN_CONFIG=../$(ZOLTAN_CONFIG)"     \
	"ZOLTAN_OBJ_DIR=$(ZOLTAN_OBJ_DIR)"
endif

vtk_view:
	@rm -f $(ZOLTAN_OBJ_DIR)/makefile
	@cp makefile_sub $(ZOLTAN_OBJ_DIR)/makefile
	@cd $(ZOLTAN_OBJ_DIR); 			\
	$(MAKE) vtk_view                        \
	"ZTARGET=vtk_view"			\
	"CC=$(CPPCOMP) "  	                \
	"MAIN_OBJ=$(VIEW_OBJ)"		        \
	"OBJ_FILES=$(VIEW_OBJ)"		        \
	"VPATH=$(ZOLTAN_VPATH)"                 \
	"DEFS=$(DEFS) $(SUBINCPATH) $(INCLUDE_PATH) $(VTK_INCPATH) $(GL_INCPATH)" \
	"LNK_LIBS=$(VIEW_LIBS)" 		\
	"ZOLTAN_CONFIG=../$(ZOLTAN_CONFIG)"     \
	"ZOLTAN_OBJ_DIR=$(ZOLTAN_OBJ_DIR)"

vtk_write:
	@rm -f $(ZOLTAN_OBJ_DIR)/makefile
	@cp makefile_sub $(ZOLTAN_OBJ_DIR)/makefile
	@cd util ; cp -f vtk_view.cpp vtk_write.cpp ; 
	@cd $(ZOLTAN_OBJ_DIR); 			\
	$(MAKE) vtk_write                       \
	"ZTARGET=vtk_write"			\
	"CC=$(CPPCOMP) "  	                \
	"OBJ_FILES=$(WRITE_OBJ)"	        \
	"VPATH=$(ZOLTAN_VPATH)"                 \
	"DEFS=$(DEFS) $(SUBINCPATH) $(INCLUDE_PATH) $(VTK_OFFSCREEN_INCPATH) $(MESA_INCPATH) -DOUTPUT_TO_FILE" \
	"MAIN_OBJ=$(WRITE_OBJ)"  		\
	"LNK_LIBS=$(OFFSCREEN_VIEW_LIBS)" 	\
	"ZOLTAN_CONFIG=../$(ZOLTAN_CONFIG)"     \
	"ZOLTAN_OBJ_DIR=$(ZOLTAN_OBJ_DIR)"

zdrive:
	$(MAKE) zoltan
	@cd $(ZOLTAN_OBJ_DIR); 			\
	$(MAKE) $(ZTARGET)			\
	"ZTARGET=$(ZTARGET)"			\
	"LIB_NAME=$(LIB_NAME)"			\
	"CC=$(CCOMP)" 				\
	"VPATH=$(ZOLTAN_VPATH)"			\
	"OBJ_FILES=$(ZOLTAN_OBJ)"		\
	"DEFS=$(GDEFINES) $(DEFS)" 		\
	"MAIN_OBJ=$(ZDRIVE_OBJ)" 		\
	"LNK_LIBS=$(ZDRIVE_LIBS)" 		\
	"DEP_LIBS=$(ZDRIVE_DEPS)"		\
	"ZOLTAN_CONFIG=../$(ZOLTAN_CONFIG)"     \
	"ZOLTAN_OBJ_DIR=$(ZOLTAN_OBJ_DIR)"

zCPPdrive:
	$(MAKE) zoltan
	@cd $(ZOLTAN_OBJ_DIR); 			\
	$(MAKE) zCPPdrive			\
	"ZTARGET=zCPPdrive"			\
	"LIB_NAME=$(LIB_NAME)"			\
	"VPATH=$(ZOLTAN_VPATH)"			\
	"CC=$(CPPCOMP)"				\
	"OBJ_FILES=$(ZOLTAN_OBJ)"		\
	"DEFS=$(GDEFINES) $(DEFS) $(ZOLTAN_INC_PATH) " 		\
	"MAIN_OBJ=$(ZCPPDRIVE_OBJ)" 		\
	"LNK_LIBS=$(ZCPPDRIVE_LIBS)" 		\
	"DEP_LIBS=$(ZCPPDRIVE_DEPS)"		\
	"ZOLTAN_CONFIG=../$(ZOLTAN_CONFIG)"     \
	"ZOLTAN_OBJ_DIR=$(ZOLTAN_OBJ_DIR)"

zfdrive:
	$(MAKE) YES_FORTRAN=1 ZOLTAN_ARCH=$(ZOLTAN_ARCH) zoltan
	@cd fdriver$(ZOLTAN_OLD); 				\
	$(MAKE) zfdrive 			\
	"CC=$(CCOMP)"				\
	"ZOLTAN_INC_PATH=$(ZOLTAN_INC_PATH)"	\
	"F90=$(F90COMP)" 			\
	"F90_MODULE_PREFIX=$(F90_MODULE_PREFIX)" \
	"FARG=$(FARG)" 				\
        "LNK_LIBS=$(ZDRIVE_LIBS)"		\
	"DEP_LIBS=$(ZDRIVE_DEPS)"		\
	"ZOLTAN_ARCH=$(ZOLTAN_ARCH)"		\
	"ZOLTAN_CONFIG=../$(ZOLTAN_CONFIG)"     \
	"ZOLTAN_OBJ_DIR=$(ZOLTAN_OBJ_DIR)"

zscript:
	cd Utilities; \
	$(MAKE) zscript
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
	@echo PARMETIS_FLAG='"$${6:+-DZOLTAN_PARMETIS}"' >> zoltan_script
	@echo PARMETIS_INCPATH='"$$7"' >> zoltan_script
	@echo JOSTLE_FLAG='"$${8:+-DZOLTAN_JOSTLE}"' >> zoltan_script
	@echo JOSTLE_INCPATH='"$$9"' >> zoltan_script
	@echo shift 9 >> zoltan_script
	@echo AR='"$$1"' >> zoltan_script
	@echo RANLIB='"$$2"' >> zoltan_script
	@echo ZOLTAN_ARCH='"$${3:-generic}"' >> zoltan_script
	@echo ZOLTAN_LIB_DIR='"$$4"' >> zoltan_script
	@echo MPI_INCPATH='"$$5"' >> zoltan_script
	@echo CCOMP='"$$CC $$CFLAGS $$F90CFLAGS $$INCLUDE_PATH $$MPI_INCPATH $$PARMETIS_INCPATH $$JOSTLE_INCPATH"' >> zoltan_script
	@echo GDEFINES='"$$PARMETIS_FLAG $$JOSTLE_FLAG"' >> zoltan_script
	@echo ZOLTAN_OBJ_DIR='"Obj_$$ZOLTAN_ARCH"' >> zoltan_script
	@echo 'if [ ! -d $$ZOLTAN_OBJ_DIR ]; ' >> zoltan_script
	@echo 'then mkdir $$ZOLTAN_OBJ_DIR ;' >> zoltan_script
	@echo 'fi' >> zoltan_script
	@echo 'rm -f $$ZOLTAN_OBJ_DIR/*' >> zoltan_script
	@echo 'cp -f Utilities/$$ZOLTAN_OBJ_DIR/$(LIB_NAME) $$ZOLTAN_OBJ_DIR' >> zoltan_script
	@echo 'cd $$ZOLTAN_OBJ_DIR' >> zoltan_script
	@echo "echo 'Zoltan:  Compiling Zoltan/zz'" >> zoltan_script
	@for iii in $(ZZ_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../zz/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/all'" >> zoltan_script
	@for iii in $(ALL_CSRC)  ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../all/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/lb'" >> zoltan_script
	@for iii in $(LB_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../lb/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/order'" >> zoltan_script
	@for iii in $(ORDER_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../order/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/rcb'" >> zoltan_script
	@for iii in $(RCB_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../rcb/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/coloring'" >> zoltan_script
	@for iii in $(COLOR_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../coloring/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/oct'" >> zoltan_script
	@for iii in $(OCT_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../oct/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/hsfc'" >> zoltan_script
	@for iii in $(HSFC_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../hsfc/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/phg'" >> zoltan_script
	@for iii in $(PHG_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../phg/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/params'" >> zoltan_script
	@for iii in $(PARAMS_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../params/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/par'" >> zoltan_script
	@for iii in $(PAR_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../par/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/parmetis'" >> zoltan_script
	@for iii in $(PARMETIS_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) '$$PARMETIS_INCPATH' -c '$$GDEFINES $$DEFS' ../parmetis/$$iii >> zoltan_script; \
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
	@echo "echo 'Zoltan:  Compiling Zoltan/hier'" >> zoltan_script
	@for iii in $(HIER_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../hier/$$iii >> zoltan_script; \
	 done;
	@echo "echo 'Zoltan:  Compiling Zoltan/simple'" >> zoltan_script
	@for iii in $(SIMPLE_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../simple/$$iii >> zoltan_script; \
	 done;
	@echo "echo Zoltan:  Linking $(LIB_NAME)" >> zoltan_script
	@echo '$$AR' $(LIB_NAME) $(ZOLTAN_OBJ) $(FOBJ_FILES) >> zoltan_script
	@echo '$$RANLIB' $(LIB_NAME) >> zoltan_script
	@echo 'if [ $$ZOLTAN_LIB_DIR ] ; ' >> zoltan_script
	@echo 'then cp -f $(LIB_NAME) $$ZOLTAN_LIB_DIR ; ' >> zoltan_script
	@echo 'fi' >> zoltan_script
	@echo "cd .." >> zoltan_script
	@chmod +x zoltan_script

# makefile target to compile ZoltanComponent. This assumesthat the CCA_SPEC_DIR and
# DCCAFE_DIR definitions in Utilities/Config/Config.<arch> exist.

zoltan_component : zoltan zdrive
	cd ZoltanComponent ; ./configure --with-ZoltanArch=$(ZOLTAN_ARCH) ; make ; cd ..

#
############################################################################
#          UTILITY TARGETS
############################################################################
#

# For lint, MPI Paths must be explicitly stated; see 
# Utilities/Config/Config.gcc for an example.

lint: $(ZDRIVE_CSRC) $(ZOLTAN_CSRC)
	@echo "running lint..."
	$(LINT) $(LINT_FLAGS) $?

printout:
	cat $(MAIN) $(ZOLTAN_CSRC) > source.c
	tps -c source.c | lpr
	-rm source.c
	tps  $(INC) | lpr

clean:
	@-for a in $(SUBDIRS) ; do \
	 cd $$a; \
         echo $$a; \
	 $(RM) -f *.o *.d *.d@ $(ZTARGET) $(ZTARGET).pure $(LIB_NAME); \
	 cd ..; \
	 done;
	@cd Utilities; $(MAKE) clean; cd .. ; 
	@echo $(ZOLTAN_OBJ_DIR);
	@$(RM) -rf $(ZOLTAN_OBJ_DIR)/*; 
	@echo "fort"
	@cd fort; $(MAKE) "ZOLTAN_OBJ_DIR=$(ZOLTAN_OBJ_DIR)" clean; cd ..;
	@echo "fdriver"
	@cd fdriver; $(MAKE) "ZOLTAN_OBJ_DIR=$(ZOLTAN_OBJ_DIR)" clean; cd ..;
	@echo "fdriver_old"
	@cd fdriver_old; $(MAKE) "ZOLTAN_OBJ_DIR=$(ZOLTAN_OBJ_DIR)" clean; cd ..;
	 exit 0

allclean:
	@-for a in $(SUBDIRS); do \
	 cd $$a; \
         echo $$a; \
	 $(RM) -f *.o *.d *.d@ $(ZTARGET) $(ZTARGET).pure $(LIB_NAME); \
	 cd ..; \
	 done;
	@-for a in $(ALL_OBJ_DIRS); do \
         echo $$a; \
	 cd fort; $(MAKE) "ZOLTAN_OBJ_DIR=$$a" clean; cd ..; \
	 cd fdriver; $(MAKE) "ZOLTAN_OBJ_DIR=$$a" clean; cd ..; \
	 cd fdriver_old; $(MAKE) "ZOLTAN_OBJ_DIR=$$a" clean; cd ..; \
	 $(RM) -rf $$a/*; \
	 done;
	@cd Utilities; $(MAKE) allclean; cd .. ; 
#	 if [ -d ZoltanComponent ]; \
#	 then cd ZoltanComponent ; $(MAKE) clean ; cd .. ; \
#	 fi;
	 exit 0

cleantest:
	 @if [ -d ../test ];         \
	 then cd ../test; \
	 touch ch_simple/zdrive.inp; \
	 touch ch_simple/output/junk; \
	 touch nem_ti_4k/pio_01/.blot.; \
	 /bin/rm -f */zdrive.inp; \
	 /bin/rm -fr */output; \
         /bin/rm -f */pio*/*.blot.*; \
	 fi;
	 exit 0
