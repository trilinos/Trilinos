include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

DEFINES= -D$(TRILINOS_ARCH) $(EPETRA_ARCH_DEFINES) -D$(EPETRA_COMM) \
         -DBTF_VERBOSE -DZOLTAN_ORDER

EPETRA_TRANSFORM_INCLUDES = -I$(TRILINOS_HOME)/contrib/epetra_transform/src
ZOLTAN_CPP_INCLUDES = -I$(TRILINOS_HOME)/external/Zoltan_C++
INCLUDES = $(ARCH_INCLUDES) $(EPETRA_INCLUDES) $(ZOLTAN_CPP_INCLUDES) $(EPETRA_TRANSFORM_INCLUDES)

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
LDFLAGS=$(ARCH_LDFLAGS)

#=======================================================================
# EPETRA_TRANSFORM source files
#=======================================================================

LIBEPETRA_TRANSFORM = $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libepetra_trans.a

EPETRA_TRANSFORM_CC = EDT_MultiVector_View.cpp \
                      ET_CrsGraph_BTF.cpp \
                      EDT_CrsGraph_Overlap.cpp \
                      EDT_CrsGraph_Transpose.cpp \
                      EDT_CrsGraph_MapColoring.cpp \
                      EDT_CrsGraph_View.cpp \
                      EDT_CrsMatrix_View.cpp 
#                      EDT_CrsGraph_Zoltan.cpp \
#                      EDT_CrsGraph_ZoltanOrder.cpp \
#                      Epetra_ZoltanQuery.cpp

EPETRA_TRANSFORM_F =  btf/concmp.f \
                      btf/genbtf.f \
                      btf/corsum.f \
                      btf/finsum.f \
                      btf/fnrsum.f \
                      btf/izero.f \
                      btf/maxmatch.f \
                      btf/rectblk.f \
                      btf/mmc13e.f \
                      btf/prtivs.f \
                      btf/strchk.f \
                      btf/mattrans.f


EPETRA_TRANSFORM_OBJ = $(EPETRA_TRANSFORM_CC:.cpp=.o) \
                       $(EPETRA_TRANSFORM_F:.f=.o)

$(LIBEPETRA_TRANSFORM): $(EPETRA_TRANSFORM_OBJ)
	mkdir -p $(TRILINOS_LIB_PATH)
	$(AR) $(ARFLAGS) $(LIBEPETRA_TRANSFORM) $(EPETRA_TRANSFORM_OBJ)
	$(RANLIB) $(LIBEPETRA_TRANSFORM)

clean:
	@rm -rf $(TRILINOS_TARGET)
	@rm -f *.o btf/*.o $(LIBEPETRA_TRANSFORM) *~ 
