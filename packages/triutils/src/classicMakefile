include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

DEFINES= -D$(TRILINOS_ARCH) -D$(EPETRA_COMM) $(EPETRA_DEFAULT_DEFINES)

INCLUDES = $(ARCH_INCLUDES) $(TRILUTILS_INCLUDES) $(EPETRA_INCLUDES)

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) 
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS=$(CXXFLAGS)
#=======================================================================
# TRIUTILS source files
#=======================================================================


TRIUTILS_CC = \
           Trilinos_Util_duscr_vbr.cpp	\
           Trilinos_Util_dusds_vbr.cpp	\
           Trilinos_Util_distrib_vbr_matrix.cpp	\
           Trilinos_Util_read_hb.cpp \
           Trilinos_Util_svbrres.cpp \
           Trilinos_Util_dusmm.cpp	\
           iohb.cpp	\
           Trilinos_Util_scscmv.cpp   \
           Trilinos_Util_write_vec.cpp \
           Trilinos_Util_create_vbr.cpp \
           Trilinos_Util_msr2vbr.cpp	\
           Trilinos_Util_scscres.cpp \
           Trilinos_Util_distrib_msr_matrix.cpp	\
           Trilinos_Util_read_coo.cpp	\
           Trilinos_Util_smsrres.cpp \
           Trilinos_Util_csrcsc.cpp \
           Trilinos_Util_csrmsr.cpp \
           Trilinos_Util_ssrcsr.cpp \
           Trilinos_Util_coocsr.cpp  \
           Trilinos_Util_ReadHb2Epetra.cpp \
           Trilinos_Util_ReadHpc2Epetra.cpp	\
           Trilinos_Util_ReadHb2EpetraVbr.cpp

TRIUTILS_C = 
TRIUTILS_F =

#=======================================================================
# TRIUTILS include files
#=======================================================================

TRIUTILS_INC =

TRIUTILS_OBJ = $(TRIUTILS_CC:.cpp=.$(OBJ_EXT)) \
                    $(TRIUTILS_C:.c=.$(OBJ_EXT)) \
                    $(TRIUTILS_F:.f=.$(OBJ_EXT))

#$(LIBTRIUTILS): $(TRIUTILS_OBJ)
LIBTRIUTILS: $(TRIUTILS_OBJ)
	mkdir -p $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)
	$(AR) $(ARFLAGS) $(LIB_WIN_FLAGS)$(LIBTRIUTILS) $(TRIUTILS_OBJ)
# following ranlib commented out by rbl
#	$(RANLIB) $(LIBTRIUTILS)

#
# dependencies for 'f' files (none at this time)
#
#include ../../build/depends.epetra

clean:
	echo "junk" > dummy.$(OBJ_EXT)
	rm -f *.$(OBJ_EXT)  $(LIBTRIUTILS) *~
