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
GDEFINES = -DDEBUG 

#
# Makefile for load-balancing library.
############################################################################
#                  SOURCE CODE
############################################################################
#
#   Paths to submodules
#
SUBVPATH=
SUBINCPATH=
SUBDIRS= Obj_sgi10k Obj_smos Obj_sun Obj_tflop 

#source code from SCCS projects

GR_CSRC=gr.c gr_hash.c gr_list.c gr_tree.c gr_type.c gr_util.c
GR_INC= gr_type.h gr_const.h gr_hash_const.h gr_list_const.h \
	gr_tree_const.h gr_type_const.h \
        gr_util_const.h

ID_CSRC=id.c id_util.c
ID_INC= id_const.h id_util.h id_util_const.h

VX_CSRC=vx.c
VX_INC= vx_const.h

ALL_CSRC=all_allo.c
ALL_INC=all_const.h all_allo_const.h

LB_CSRC=lb.c
LB_INC=lb_const.h lb.h lbi_const.h

# machine independent routines
INC          = $(GR_INC) $(ID_INC) $(VX_INC) $(ALL_INC) 
#LIB_SRC_C    = $(GR_CSRC) $(ID_CSRC) $(VX_CSRC) $(ALL_CSRC) 
LIB_SRC_C    = $(LB_CSRC) $(ALL_CSRC) 
LIB_SRC_OBJ  = $(LIB_SRC_C:.c=.o) 

#library paths for different machines

LIB_PATH    =   /Net/hondo/raid1/rf_crew/arch
TFLOP_PATH  =   $(LIB_PATH)/tflop
SMOS_PATH   =   $(LIB_PATH)/smos
SUN_PATH    =   $(LIB_PATH)/sun
SGI10K_PATH =   $(LIB_PATH)/sgi10k



SUN_INC_PATH     =   -I$(SUN_PATH)/include  -I/Net/local/lang/SC2.0.1/include \
			-I/Net/local/lang/SC2.0.1/ansi_include $(SUBINCPATH)
SGI10K_INC_PATH  =   -I$(SGI10K_PATH)/include $(SUBINCPATH)
SMOS_INC_PATH    =   -I$(SMOS_PATH)/include $(SUBINCPATH)
TFLOP_INC_PATH   =   -I$(TFLOP_PATH)/include $(SUBINCPATH)


###########################################################################
#      TARGETS FOR DIFFERENT MACHINES
###########################################################################
#
smos:
	@cd Obj_smos; \
	$(MAKE) dlb "SMOS_OBJ=$(LIB_SRC_OBJ)" "INC=$(SMOS_INC)" "TARGET=dlb" \
	"LIB=$(LIB)" "SMOS_LIB_PATH=$(SMOS_LIB_PATH)" \
	"SMOS_INC_PATH=$(SMOS_INC_PATH)" "VPATH=$(SUBVPATH)" \
        "GDEFINES=$(GDEFINES)"

tflop:
	@cd Obj_tflop; \
	$(MAKE) dlb "TFLOP_OBJ=$(LIB_SRC_OBJ)" "INC=$(TFLOP_INC)" \
	"TARGET=dlb" "LIB=$(LIB)" "TFLOP_LIB_PATH=$(TFLOP_LIB_PATH)" \
	"TFLOP_INC_PATH=$(TFLOP_INC_PATH)" "VPATH=$(SUBVPATH)" \
        "GDEFINES=$(GDEFINES)"

sun:
	@cd Obj_sun; \
	$(MAKE) dlb "SUN_OBJ=$(LIB_SRC_OBJ)" "INC=$(INC)" \
	"TARGET=dlb" "LIB=$(LIB) $(SUN_BLASLIB)" \
	"SUN_LIB_PATH=$(SUN_LIB_PATH)" \
	"SUN_INC_PATH=$(SUN_INC_PATH)" "VPATH=$(SUBVPATH)" \
        "GDEFINES=$(GDEFINES)"
#	@cscope -b $(MAIN) $(SUN_CSRC) $(INC)

sgi10k:
	@cd Obj_sgi10k; \
	$(MAKE) dlb "SGI_OBJ=$(LIB_SRC_OBJ)" "INC=$(INC)" "TARGET=dlb" \
	"LIB=$(LIB)" "SGI_LIB_PATH=$(SGI10K_LIB_PATH)" \
	"SGI_INC_PATH=$(SGI10K_INC_PATH)" "VPATH=$(SUBVPATH)" \
        "GDEFINES=$(GDEFINES)"

#
############################################################################
#          UTILITY TARGETS
############################################################################
#

printout:
	cat $(MAIN) $(LIB_SRC_C) > source.c
	tps -c source.c | lpr
	-rm source.c
	tps  $(INC) | lpr

clean:
	-for a in $(SUBDIRS); do cd $$a; \
	 $(RM) *.o *.d *.d@ dlb dlb.pure; \
	 cd ..; done; exit 0
