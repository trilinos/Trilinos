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
# Top-level Makefile for Zoltan utilities.
#
# Utilities currently supported:
#    Memory management utility (libzoltan_mem.a)
#    Unstructured communication utility (libzoltan_comm.a)
#
# To compile a utility, type
#    gmake ZOLTAN_ARCH=<targetmachine> <Utility name>
# where
#    <targetmachine> = solaris, gcc, sgi64, cplant, or tflop,
# and 
#    <Utility name> can be a library name 
#        libzoltan_comm.a  
#        libzoltan_mem.a
#    or a driver name
#        zdrive_comm (or comm, for short)
#        zdrive_mem  (or mem, for short)
#
##############################################################################

##############################################################################
# Process architecture flag
##############################################################################

ifndef ZOLTAN_ARCH
ZOLTAN_ARCH	= generic
endif

ifndef CONFIG
CONFIG = Config/Config.$(ZOLTAN_ARCH)
endif

include $(CONFIG)

OBJ_DIR = Obj_$(ZOLTAN_ARCH)
ALL_OBJ_DIRS := $(wildcard Obj_*)

ifndef LIB_DIR
LIB_DIR = $(OBJ_DIR)
endif

############################################################################
#                  SOURCE CODE
############################################################################
#
#   Paths to submodules
#
SUB_VPATH	= ../Memory ../Communication
SUBINCPATH	= -I../Memory -I../Communication
SUBDIRS		= Memory Communication

############################################################################
# MEMORY MANAGEMENT PACKAGE

MEM_CSRC	= mem.c
MEM_INC		= mem_const.h
MEM_LIB		= libzoltan_mem.a

MEM_DRIVER	= zdrive_mem
MEM_DRIVER_CSRC	= mem_main.c
MEM_DRIVER_LIBS = $(MEM_LIB) $(MPI_LIBPATH) $(MPI_LIB) 

############################################################################
# COMMUNICATION PACKAGE

COMM_CSRC	= comm_create.c comm_do.c comm_do_reverse.c comm_destroy.c \
		  comm_invert_map.c
COMM_INC	= comm_const.h
COMM_LIB	= libzoltan_comm.a

COMM_DRIVER	= zdrive_comm
COMM_DRIVER_CSRC= comm_main.c
COMM_DRIVER_LIBS= $(COMM_LIB) $(MEM_LIB) $(MPI_LIBPATH) $(MPI_LIB) 


##############################################################################
#  Architecture-independent definitions
##############################################################################

INC_PATH	= $(SUBINCPATH) $(INCLUDE_PATH) 

CCOMP		= $(CC) $(INC_PATH) $(CFLAGS) 

############################################################################
#                  LINT DEFINITIONS KDDKDD
############################################################################
#
# lint defines
#
LINT		= lint
LINT_INC_PATH	= -IMemory -ICommunication \
		  -I$(INC_PATH)
#		  -I/Net/local/lang/SC3.0.1/include \
#		  -I/Net/local/lang/SC3.0.1/ansi_include \
#		  -I/Net/local/mpi/include \
#		  -I/Net/local/mpi/mpich/include \
#                  -I/Net/local/mpi/build/solaris/ch_p4/include \
#		  -I/Net/proj/zoltan/arch/sun/include
LINT_FLAGS	= $(LINT_INC_PATH) $(GDEFINES)
LINT_TARGET	= llib-l$(LIB_TARGET).ln
VPATH		= $(SUBDIRS) 

###########################################################################
#      HELP HELP HELP HELP HELP  
###########################################################################
help:
	@echo  ' '
	@echo  'Top-level Makefile for Zoltan utilities.'
	@echo  ' '
	@echo  'Utilities currently supported:'
	@echo  '   Memory management utility ('$(MEM_LIB)')'
	@echo  '   Unstructured communication utility ('$(COMM_LIB)')'
	@echo  ' '
	@echo  'To compile a utility, type'
	@echo  '   gmake ZOLTAN_ARCH=<targetmachine> <Utility name>'
	@echo  'where'
	@echo  '   <targetmachine> = solaris, gcc, sgi64, cplant, or tflop,'
	@echo  'and '
	@echo  '   <Utility name> can be a library name '
	@echo  '        '$(COMM_LIB)
	@echo  '        '$(MEM_LIB)
	@echo  'or a driver name'
	@echo  '        '$(COMM_DRIVER)' (or comm, for short)'
	@echo  '        '$(MEM_DRIVER)' (or mem, for short)'
	@echo  ' '
	@echo  'The file Config/Config.<targetmachine> '
	@echo  'should exist.  See Config/Config.generic for '
	@echo  'instructions on creating Config/Config.<targetmachine>.'
	@echo  ' '
	@echo  'Example: '
	@echo  '   gmake ZOLTAN_ARCH=tflop comm '
	@echo  'builds '$(COMM_DRIVER)' for the tflop (ASCI Red) machine.'
	@echo  ' '


###########################################################################
#      TARGETS FOR DIFFERENT MACHINES
###########################################################################

$(MEM_LIB):
	@if [ ! -d $(OBJ_DIR) ]; 		\
	then mkdir $(OBJ_DIR); 			\
	fi
	@rm -f $(OBJ_DIR)/Makefile
	@cp Makefile_sub $(OBJ_DIR)/Makefile
	@cd $(OBJ_DIR); 			\
	$(MAKE) $(MEM_LIB) 			\
	"LIB_NAME=$(MEM_LIB)"			\
	"CC=$(CCOMP)"				\
	"VPATH=$(SUB_VPATH)"			\
	"OBJ_FILES=$(MEM_CSRC:.c=.o)"		\
	"DEFS=$(GDEFINES) $(DEFS)" 		\
	"AR=$(AR)" 				\
	"RANLIB=$(RANLIB)" 			

zdrive_mem mem:
	@if [ ! -d $(OBJ_DIR) ]; 		\
	then mkdir $(OBJ_DIR); 			\
	fi
	$(MAKE) $(MEM_LIB) 
	@rm -f $(OBJ_DIR)/Makefile
	@cp Makefile_sub $(OBJ_DIR)/Makefile
	@cd $(OBJ_DIR); 			\
	$(MAKE) zdrive_mem 			\
	"DRIVER_NAME=zdrive_mem"		\
	"CC=$(CCOMP)"				\
	"VPATH=$(SUB_VPATH)"			\
	"OBJ_FILES=$(MEM_DRIVER_CSRC:.c=.o)" 	\
	"DEFS=$(GDEFINES) $(DEFS)" 		\
	"LNK_LIBS=$(MEM_DRIVER_LIBS)"

$(COMM_LIB):
	@if [ ! -d $(OBJ_DIR) ]; 		\
	then mkdir $(OBJ_DIR); 			\
	fi
	@rm -f $(OBJ_DIR)/Makefile
	@cp Makefile_sub $(OBJ_DIR)/Makefile
	@cd $(OBJ_DIR); 			\
	$(MAKE) $(COMM_LIB)			\
	"LIB_NAME=$(COMM_LIB)"			\
	"CC=$(CCOMP)"				\
	"VPATH=$(SUB_VPATH)"			\
	"OBJ_FILES=$(COMM_CSRC:.c=.o)" 		\
	"DEFS=$(GDEFINES) $(DEFS)" 		\
	"AR=$(AR)" 				\
	"RANLIB=$(RANLIB)" 			

zdrive_comm comm:
	@if [ ! -d $(OBJ_DIR) ]; 		\
	then mkdir $(OBJ_DIR); 			\
	fi
	$(MAKE) $(MEM_LIB)
	$(MAKE) $(COMM_LIB)
	@rm -f $(OBJ_DIR)/Makefile
	@cp Makefile_sub $(OBJ_DIR)/Makefile
	@cd $(OBJ_DIR); 			\
	$(MAKE) zdrive_comm 			\
	"DRIVER_NAME=zdrive_comm"		\
	"CC=$(CCOMP)"				\
	"VPATH=$(SUB_VPATH)"			\
	"OBJ_FILES=$(COMM_DRIVER_CSRC:.c=.o)" 	\
	"DEFS=$(GDEFINES) $(DEFS)" 		\
	"LNK_LIBS=$(COMM_DRIVER_LIBS)"
	
#
############################################################################
#         OTHER TARGETS
############################################################################
#

$(LINT_TARGET): $(ZOLTAN_CSRC)
	@echo "Generating lint file for library..."
	$(LINT) $(LINT_FLAGS)  $? -o $(LIB_TARGET)

lint: $(ZDRIVE_CSRC) $(LINT_TARGET)
	@echo "running lint..."
	$(LINT) $(LINT_FLAGS) "VPATH=$(LINT_VPATH) $?

clean:
	-for a in $(SUBDIRS); do \
	 cd $$a; \
         echo $$a; \
	 $(RM) -f *.o *.d *.d@ lib*; \
	 cd ..; \
	 done;
	-for a in $(ALL_OBJ_DIRS); do \
         echo $$a; \
	 $(RM) -rf $$a; \
	 done;
	 $(RM) $(LINT_TARGET)
	 exit 0
