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
#    Distributed directory utility (libzoltan_dd.a)
#
# To compile a utility, type
#    gmake ZOLTAN_ARCH=<targetmachine> <Utility name>
# where
#    <targetmachine> = suffix of one of the configuration files
# Config.* in Zoltan/Utilities/Config
# and 
#    <Utility name> can be a library name 
#        libzoltan_comm.a  
#        libzoltan_mem.a
#        libzoltan_dd.a
#    or a driver name
#        zdrive_comm (or comm, for short)
#        zdrive_mem  (or mem, for short)
#        zdrive_dd   (or dd, for short)
#
##############################################################################

##############################################################################
# Process architecture flag
##############################################################################

ifndef ZOLTAN_ARCH
ZOLTAN_ARCH	= generic
endif

ifndef ZOLTAN_CONFIG
ZOLTAN_CONFIG = Config/Config.$(ZOLTAN_ARCH)
endif

include $(ZOLTAN_CONFIG)

OBJ_DIR = Obj_$(ZOLTAN_ARCH)
ALL_OBJ_DIRS := $(wildcard Obj_*)

ifdef ZOLTAN_LIB_DIR
LIB_DIR = $(ZOLTAN_LIB_DIR)
else
LIB_DIR = $(OBJ_DIR)
endif

ifndef ZOLTAN_LIB
ZOLTAN_LIB = libzoltan.a
endif

############################################################################
#                  SOURCE CODE
############################################################################
#
#   Paths to submodules
#
SUB_VPATH	= ../../include ../shared \
		  ../Memory ../Communication ../DDirectory 
SUBINCPATH	= -I../../include -I../shared \
		  -I../Memory -I../Communication -I../DDirectory 
SUBDIRS		= Memory Communication DDirectory shared

############################################################################
# MEMORY MANAGEMENT PACKAGE

MEM_CSRC	= mem.c
MEM_INC		= zoltan_mem.h
MEM_LIB		= libzoltan_mem.a
MEM_TARGET	= zoltan_mem

MEM_DRIVER	= zdrive_mem
MEM_DRIVER_CSRC	= mem_main.c
MEM_DRIVER_LIBS = -L. -l$(MEM_TARGET) $(MPI_LIBPATH) $(MPI_LIB) 
MEM_DRIVER_DEP  = $(MEM_LIB)

############################################################################
# COMMUNICATION PACKAGE

COMM_CSRC	= comm_create.c comm_do.c comm_do_reverse.c comm_destroy.c \
		  comm_invert_map.c comm_exchange_sizes.c comm_resize.c \
                  comm_sort_ints.c
COMM_INC	= comm.h zoltan_comm.h
COMM_LIB	= libzoltan_comm.a
COMM_TARGET     = zoltan_comm

COMM_DRIVER	= zdrive_comm
COMM_DRIVER_CSRC= comm_main.c
COMM_DRIVER_LIBS= -L. -l$(COMM_TARGET) -l$(MEM_TARGET) $(MPI_LIBPATH) $(MPI_LIB) 
COMM_DRIVER_DEP = $(COMM_LIB) $(MEM_LIB)

############################################################################
# UTILITIES
UTIL_CSRC		= zoltan_id.c zoltan_align.c
UTIL_INC		= zoltan_id.h

  
############################################################################
# DISTRIBUTED DIRECTORY PACKAGE

DD_CSRC		= DD_Create.c DD_Destroy.c DD_Find.c DD_Remove.c DD_Update.c \
                  DD_Set_Hash_Fn.c DD_Hash2.c DD_Stats.c DD_Print.c\
                  DD_Set_Neighbor_Hash_Fn1.c DD_Set_Neighbor_Hash_Fn2.c \
		  $(UTIL_CSRC)
DD_INC		= zoltan_dd.h
DD_LIB		= libzoltan_dd.a
DD_TARGET	= zoltan_dd

DD_DRIVER	= zdrive_dd
DD_DRIVER_CSRC	= DD_Main.c
DD_DRIVER_LIBS 	= -L. -l$(DD_TARGET) -l$(COMM_TARGET) -l$(MEM_TARGET) \
                  $(MPI_LIBPATH) $(MPI_LIB)
DD_DRIVER_DEP  	= $(DD_LIB) $(COMM_LIB) $(MEM_LIB) $(ID_OBJ)


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
LINT_INC_PATH	= -Ishared -IMemory -ICommunication -IDDirectory $(INC_PATH)
LINT_FLAGS	= $(LINT_INC_PATH) $(GDEFINES)
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
	@echo  '   Distributed directory utility ('$(DD_LIB)')'
	@echo  ' '
	@echo  'To compile a utility, type'
	@echo  '   gmake ZOLTAN_ARCH=<targetmachine> <Utility name>'
	@echo  'where'
	@echo  '   <targetmachine> = the suffix of one of the configuration'
	@echo  '                     files Config.* in Zoltan/Utilities/Config.'
	@echo  'and '
	@echo  '   <Utility name> can be a library name '
	@echo  '        '$(COMM_LIB)
	@echo  '        '$(MEM_LIB)
	@echo  '        '$(DD_LIB)
	@echo  'or a driver name'
	@echo  '        '$(COMM_DRIVER)' (or comm, for short)'
	@echo  '        '$(MEM_DRIVER)' (or mem, for short)'
	@echo  '        '$(DD_DRIVER)' (or dd, for short)'
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
	@if [ ! -d $(LIB_DIR) ]; 		\
	then mkdir $(LIB_DIR); 			\
	fi
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
	"RANLIB=$(RANLIB)" 			\
	"ZOLTAN_CONFIG=../$(ZOLTAN_CONFIG)"     \
        "BUILD_DIR=./"
	@if [ $(OBJ_DIR) != $(LIB_DIR) ];            \
	then cp -f $(OBJ_DIR)/$(MEM_LIB) $(LIB_DIR); \
	fi


zdrive_mem mem:
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
	"ZOLTAN_CONFIG=../$(ZOLTAN_CONFIG)"     \
	"LNK_LIBS=$(MEM_DRIVER_LIBS)"		\
	"DEP_LIBS=$(MEM_DRIVER_DEP)"

$(COMM_LIB):
	@if [ ! -d $(LIB_DIR) ]; 		\
	then mkdir $(LIB_DIR); 			\
	fi
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
	"RANLIB=$(RANLIB)" 			\
	"ZOLTAN_CONFIG=../$(ZOLTAN_CONFIG)"     \
        "BUILD_DIR=./"
	@if [ $(OBJ_DIR) != $(LIB_DIR) ];                \
	then cp -f $(OBJ_DIR)/$(COMM_LIB) $(LIB_DIR);    \
	fi

zdrive_comm comm: 
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
	"ZOLTAN_CONFIG=../$(ZOLTAN_CONFIG)"     \
	"LNK_LIBS=$(COMM_DRIVER_LIBS)"		\
	"DEP_LIBS=$(COMM_DRIVER_DEP)"
	
$(DD_LIB):
	@if [ ! -d $(LIB_DIR) ]; 		\
	then mkdir $(LIB_DIR); 			\
	fi
	@if [ ! -d $(OBJ_DIR) ]; 		\
	then mkdir $(OBJ_DIR); 			\
	fi
	@rm -f $(OBJ_DIR)/Makefile
	@cp Makefile_sub $(OBJ_DIR)/Makefile
	@cd $(OBJ_DIR); 			\
	$(MAKE) $(DD_LIB) 			\
	"LIB_NAME=$(DD_LIB)"			\
	"CC=$(CCOMP)"				\
	"VPATH=$(SUB_VPATH)"			\
	"OBJ_FILES=$(DD_CSRC:.c=.o)"		\
	"DEFS=$(GDEFINES) $(DEFS)" 		\
	"AR=$(AR)" 				\
	"RANLIB=$(RANLIB)" 			\
	"ZOLTAN_CONFIG=../$(ZOLTAN_CONFIG)"     \
        "BUILD_DIR=./"
	@if [ $(OBJ_DIR) != $(LIB_DIR) ];            \
	then cp -f $(OBJ_DIR)/$(DD_LIB) $(LIB_DIR); \
	fi


zdrive_dd dd:
	$(MAKE) $(MEM_LIB)
	$(MAKE) $(COMM_LIB)
	$(MAKE) $(DD_LIB) 
	@rm -f $(OBJ_DIR)/Makefile
	@cp Makefile_sub $(OBJ_DIR)/Makefile
	@cd $(OBJ_DIR); 			\
	$(MAKE) zdrive_dd 			\
	"DRIVER_NAME=zdrive_dd"		\
	"CC=$(CCOMP)"				\
	"VPATH=$(SUB_VPATH)"			\
	"OBJ_FILES=$(DD_DRIVER_CSRC:.c=.o)" 	\
	"DEFS=$(GDEFINES) $(DEFS)" 		\
	"ZOLTAN_CONFIG=../$(ZOLTAN_CONFIG)"     \
	"LNK_LIBS=$(DD_DRIVER_LIBS)"		\
	"DEP_LIBS=$(DD_DRIVER_DEP)"

$(ZOLTAN):
	$(MAKE) $(MEM_LIB)
	$(MAKE) $(COMM_LIB)
	$(MAKE) $(DD_LIB) 
	@rm -f $(OBJ_DIR)/Makefile
	@cp Makefile_sub $(OBJ_DIR)/Makefile
	@cd $(OBJ_DIR); 			\
	$(MAKE) $(ZOLTAN_LIB) 			\
	"LIB_NAME=$(ZOLTAN_LIB)"		\
	"OBJ_FILES=$(DD_CSRC:.c=.o) $(MEM_CSRC:.c=.o) $(COMM_CSRC:.c=.o)" \
	"AR=$(AR)" 				\
	"RANLIB=$(RANLIB)" 			\
	"ZOLTAN_CONFIG=../$(ZOLTAN_CONFIG)"     \
        "BUILD_DIR=../../$(OBJ_DIR)/"

zscript:
	@echo "Building Utilities zoltan_script..."
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
	@echo GDEFINES='"$$PARMETIS_FLAG $$JOSTLE_FLAG"' >> zoltan_script
	@echo ZOLTAN_OBJ_DIR='"Obj_$$ZOLTAN_ARCH"' >> zoltan_script
	@echo 'if [ ! -d $$ZOLTAN_OBJ_DIR ]; ' >> zoltan_script
	@echo 'then mkdir $$ZOLTAN_OBJ_DIR ;' >> zoltan_script
	@echo 'fi' >> zoltan_script
	@echo 'rm -f $$ZOLTAN_OBJ_DIR/*' >> zoltan_script
	@echo 'cd $$ZOLTAN_OBJ_DIR' >> zoltan_script
	@echo "echo 'Zoltan:  Compiling Utilities/Memory'" >> zoltan_script
	@for iii in $(MEM_CSRC)  ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../Memory/$$iii >> zoltan_script; \
	 done;
	@echo "echo Zoltan:  Linking $(MEM_LIB)" >> zoltan_script
	@echo '$$AR' $(MEM_LIB) $(MEM_CSRC:.c=.o) >> zoltan_script
	@echo '$$RANLIB' $(MEM_LIB) >> zoltan_script
	@echo 'if [ $$ZOLTAN_LIB_DIR ] ; ' >> zoltan_script
	@echo 'then cp -f $(MEM_LIB) $$ZOLTAN_LIB_DIR ; ' >> zoltan_script
	@echo 'else cp -f $(MEM_LIB) ../../$$ZOLTAN_OBJ_DIR ; ' >> zoltan_script
	@echo 'fi' >> zoltan_script
	@echo "echo 'Zoltan:  Compiling Utilities/Communication'" >> zoltan_script
	@for iii in $(COMM_CSRC) ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../Communication/$$iii >> zoltan_script; \
	 done;
	@echo "echo Zoltan:  Linking $(COMM_LIB)" >> zoltan_script
	@echo '$$AR' $(COMM_LIB) $(COMM_CSRC:.c=.o) >> zoltan_script
	@echo '$$RANLIB' $(COMM_LIB) >> zoltan_script
	@echo 'if [ $$ZOLTAN_LIB_DIR ] ; ' >> zoltan_script
	@echo 'then cp -f $(COMM_LIB) $$ZOLTAN_LIB_DIR ; ' >> zoltan_script
	@echo 'else cp -f $(COMM_LIB) ../../$$ZOLTAN_OBJ_DIR ; ' >> zoltan_script
	@echo 'fi' >> zoltan_script
	@echo "echo 'Zoltan:  Compiling Utilities/DDirectory'" >> zoltan_script
	@for iii in $(DD_CSRC)  ; do \
	 echo '$$CCOMP' $(SUBINCPATH) -c '$$GDEFINES $$DEFS' ../DDirectory/$$iii >> zoltan_script; \
	 done;
	@echo "echo Zoltan:  Linking $(DD_LIB)" >> zoltan_script
	@echo '$$AR' $(DD_LIB) $(DD_CSRC:.c=.o) >> zoltan_script
	@echo '$$RANLIB' $(DD_LIB) >> zoltan_script
	@echo 'if [ $$ZOLTAN_LIB_DIR ] ; ' >> zoltan_script
	@echo 'then cp -f $(DD_LIB) $$ZOLTAN_LIB_DIR ; ' >> zoltan_script
	@echo 'else cp -f $(DD_LIB) ../../$$ZOLTAN_OBJ_DIR ; ' >> zoltan_script
	@echo 'fi' >> zoltan_script
	@echo "cd .." >> zoltan_script
	@chmod +x zoltan_script
#
############################################################################
#         LINT TARGETS
############################################################################
#
# Note:  The path for MPI include files must be explicitly specified for
#        lint targets; using ZOLTAN_ARCH=gcc with the gmake command 
#        automatically defines these MPI paths.  For example,
#           gmake ZOLTAN_ARCH=gcc lint_comm

lint_mem: $(MEM_CSRC) $(MEM_DRIVER_CSRC) 
	@echo "running lint..."
	$(LINT) $(LINT_FLAGS) $?

lint_comm: $(MEM_CSRC) $(COMM_CSRC) $(COMM_DRIVER_CSRC) 
	@echo "running lint..."
	$(LINT) $(LINT_FLAGS) $?

lint_dd: $(MEM_CSRC) $(COMM_CSRC) $(DD_CSRC) $(DD_DRIVER_CSRC) 
	@echo "running lint..."
	$(LINT) -I../include $(LINT_FLAGS)  $?

#
############################################################################
#         OTHER TARGETS
############################################################################
#

clean:
	-for a in $(SUBDIRS); do \
	 cd $$a; \
         echo $$a; \
	 $(RM) -f *.o *.d *.d@ lib*; \
	 cd ..; \
	 done;
	-for a in $(ALL_OBJ_DIRS); do \
         echo $$a; \
	 $(RM) -rf $$a/*; \
	 done;
	 $(RM) $(LINT_TARGET)
	 exit 0
