# To run makefile:
#    1) set environment variable TRILINOS_ARCH to sgi, sun, tflop, or pclinux.
#       Other machines require an appropriate makefile.$(TRILINOS_ARCH) file.
#    2) Set TRILINOS_COMM to SERIAL or MPI
#    3) (Optional) Set TRILINOS_ID to make unique version for same 
#       architecture and communication mode.
#
#    4) Make and run the epetra esi unit tests by typing 'make'.
#
include $(TRILINOS_HOME)/build/TARGET

TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)

LIBEPETRA_ESI= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libepetra_esi.a
LIBEPETRA= $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libepetra.a
LIBAZTEC=$(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libaztecoo.a
LIBY12M=$(TRILINOS_HOME)/external/y12m/lib/$(TRILINOS_TARGET)/liby12m.a

include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

ESI_PATH = $(TRILINOS_HOME)/external/esi/cxx/include
EPETRA_ESI_PATH = $(TRILINOS_HOME)/adaptors/epetra_esi/src
AZTEC_PATH = $(TRILINOS_HOME)/packages/aztecoo/src

# Petra communication defines
EPETRA_COMM_SERIAL          = SERIAL
EPETRA_COMM_MPI             = EPETRA_MPI
EPETRA_COMM                 = $(EPETRA_COMM_$(TRILINOS_COMM))

DEFINES= -D$(TRILINOS_ARCH) $(EPETRA_ARCH_DEFINES) -D$(EPETRA_COMM)
#  -DEPETRA_ESI_INCLUDE_IMPLEMENTATION

INCLUDES = $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/packages/epetra/src $(BLAS_INCLUDES) -I$(AZTEC_PATH) -I$(ESI_PATH) -I$(EPETRA_ESI_PATH)

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS=$(CXXFLAGS)

purify=purify -best-effort -cache-dir=$(HOME)/purify-cache -chain-length=15
purify=


LIB_PATHS= $(LIBEPETRA_ESI) $(LIBAZTEC) $(LIBEPETRA) $(LIBY12M) $(LIBLAPACK) $(LIBBLAS) 

#=======================================================================
# Petra test source files
#=======================================================================

TEST_CXX = \
	epetra_esi_object.cpp \
	epetra_esi_indexspace.cpp \
	epetra_esi_vector.cpp \
	epetra_esi_crsmatrix_1.cpp \
	epetra_esi_crsmatrix_2.cpp \
	aztec_esi_solver.cpp

#=======================================================================
# TEST include files
#=======================================================================

TEST_INC = 

TEST_OBJ = $(TEST_CXX:.cpp=.o)

TARGETS = $(TEST_OBJ:.o=.$(TRILINOS_COMM).x)

default: tests

tests:
	@for i in $(TARGETS) ; do \
	  $(MAKE) LINKER="$(LINKER)" CXX="$(CXX)" \
		 ARCH_CXXFLAGS="$(ARCH_CXXFLAGS)" $$i ; \
	done
	(run_epetra_esi_tests)

.SUFFIXES:
.SUFFIXES: .cpp .o .MPI.x .SERIAL.x

.cpp.o:
	@echo " "
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.cpp.$(TRILINOS_COMM).x:
	@echo " "
	$(purify) $(LINKER) $(CXXFLAGS) \
	   -o $@ $< \
	   $(LIB_PATHS) $(ARCH_LIBS) $(LIBMPI)
	@rm -f *.o */*.o */*/*.o

.o.$(TRILINOS_COMM).x:
	@echo " "
	@echo "Making " $@
	$(purify) $(LINKER) $(ARCH_FLAGS) \
	   -o $@ $< \
	   $(LIB_PATHS) $(ARCH_LIBS) $(LIBMPI)
	@rm -f *.o */*.o */*/*.o

clean:
	@rm -f PI* *.o  *~ *.MPI.x *.SERIAL.x core *.out
	@rm -rf SunWS_cache *.ti

