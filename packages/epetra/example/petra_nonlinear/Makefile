include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

DEFINES= -D$(TRILINOS_ARCH) $(PETRA_ARCH_DEFINES) -D$(PETRA_COMM) \
         -DIFPACK -DWITHSPBLAS

INCLUDES = $(ARCH_INCLUDES) -I$(TRILINOS_HOME)/packages/epetra/src \
           -I$(TRILINOS_HOME)/packages/aztecoo/src

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS=$(CXXFLAGS)


LIB_PATHS= $(LIBPETRA) $(LIBAZTEC) $(LIBIFPACK) \
           $(LIBLAPACK) $(LIBSPBLAS) $(LIBBLAS) $(LIBY12M)

#=======================================================================
# Petra test source files
#=======================================================================

TEST_CC = cc_main.cc basis.cc 
TEST_C = 
TEST_F = 

#=======================================================================
# TEST include files
#=======================================================================

TEST_INC = 

TEST_OBJ          = $(TEST_CC:.cc=.o) $(TEST_C:.c=.o) $(TEST_F:.f=.o)


TARGET_MPI = run.exe
TARGET_SERIAL = petra_nonlinear_serial.exe
TARGET = $(TARGET_$(TRILINOS_COMM))


$(TARGET): $(TEST_OBJ)
	@echo "Making " $(TARGET)
	$(LINKER) $(ARCH_FLAGS) $(TEST_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI) -o $(TARGET)

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  *~ $(TARGET_MPI) $(TARGET_SERIAL)
