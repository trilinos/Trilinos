include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

DEFINES=

ZOLTAN_CPP_INCLUDES = -I$(TRILINOS_HOME)/external/Zoltan_C++
INCLUDES = $(ARCH_INCLUDES) $(ZOLTAN_CPP_INCLUDES)

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
LDFLAGS=$(ARCH_LDFLAGS)

#=======================================================================
# ZOLTAN_CPP source files
#=======================================================================

LIBZOLTAN_CPP = $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libzoltan_cpp.a

ZOLTAN_CPP_CC = Zoltan_LoadBalance.cpp \
                Zoltan_MigrationContainer.cpp \
                Zoltan_MigrationFunctions.cpp \
                Zoltan_MigrationObject.cpp \
                Zoltan_QueryContainer.cpp \
                Zoltan_QueryFunctions.cpp \
                Zoltan_QueryObject.cpp

ZOLTAN_CPP_OBJ = $(ZOLTAN_CPP_CC:.cpp=.o) \

$(LIBZOLTAN_CPP): $(ZOLTAN_CPP_OBJ)
	mkdir -p $(TRILINOS_LIB_PATH)
	$(AR) $(ARFLAGS) $(LIBZOLTAN_CPP) $(ZOLTAN_CPP_OBJ)
	$(RANLIB) $(LIBZOLTAN_CPP)

clean:
	@rm -rf $(TRILINOS_TARGET)
	@rm -f *.o btf/*.o $(LIBZOLTAN_CPP) *~ 
