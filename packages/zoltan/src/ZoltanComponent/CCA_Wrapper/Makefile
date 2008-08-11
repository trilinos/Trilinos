## Makefile for ZoltanComponent. The user need not hack this. All the user-configurable
## stuff comes in via MakeIncl.ZoltanComponent, in the directory above.

# Suck in the user-configuarble stuff
include ../MakeIncl.ZoltanComponent

## suck in some CCA specific compile flags
include $(CCA_MAKEINCL_DIR)/MakeIncl.CCA_Component

.SUFFIXES: .c .cpp .o .so.o

SAPPOBJ = BaseLB.so.o PartitonerFactory.so.o RCB_prop.so.o BaseLB_prop.so.o RCB.so.o \
	  RIB.so.o RIB_prop.so.o HSFC.so.o HSFC_prop.so.o OctTree.so.o OctTree_prop.so.o \
	  ParMetis.so.o ParMetis_prop.so.o  PartitionerFactory_DataMig.so.o \
	  ZoltanApplQueryFunc.so.o ZoltanLB_wrapper.so.o 

libshared = ZoltanLB$(CCA_SHLIB_SUFFIX)

# Application specific Flags
C++APPFLAGS = -Wno-unused  -DNDEBUG

# Do we need to link in ParMetis and Jostle ?
ifdef PARMETIS_LIBPATH
PARMETIS_LIBS = -lparmetis -lmetis
endif

ifdef JOSTLE_LIBPATH
JOSTLE_LIBS = -ljostle
endif


# Application specific libraries like blas, lapack etc.
APPLIB = -L$(DCCAFE_LIB)/ -Wl,-rpath,$(DCCAFE_LIB) -l$(DCCAFE_LIBNAME) \
	 -L$(ZOLTAN_HOME)/Obj_$(ZOLTAN_ARCH) -lzoltan \
	 $(PARMETIS_LIBPATH) $(PARMETIS_LIBS) \
	 $(JOSTLE_LIBPATH) $(JOSTLE_LIBS)

LDLIBS = 
#-lstd++

# MPI and CCA specific includes
INCLUDE1 = -I$(MPI_HOME)/include
INCLUDE2 = $(CCA_INC)
INCLUDE3 = -I$(ZOLTAN_HOME) -I$(ZOLTAN_HOME)/include \
           -I$(ZOLTAN_HOME)/Utilities/Communication -I$(ZOLTAN_HOME)/Utilities/Memory \
	   -I$(ZOLTAN_HOME)/Utilities/DDirectory
INCLUDE4 = -I $(DCCAFE_INCLUDE)
INCLUDE5 = -I $(ZOLTAN_PORTS_DIR)/
INCLUDE = $(INCLUDE1) $(INCLUDE2) $(INCLUDE3) $(INCLUDE4) $(INCLUDE5) 

# Misc. definitions
ECHO = echo
RM = /bin/rm -f

default : component

component : $(SAPPOBJ)
	$(RM) $(libshared)
	$(CCA_SHLIB_LD)  -o $(libshared) $(SAPPOBJ) $(APPLIB) $(LDLIBS)

.cpp.so.o : 
	@$(ECHO)
	@$(ECHO) "Compiling Source File --" $<
	@$(ECHO) "---------------------"
	$(CCA_CXX) $(CCA_CXX_FLAGS) $(CCA_SHLIB_CFLAGS) $(C++APPFLAGS)  -c $(INCLUDE) -o $@ $< 
	@$(ECHO)

clean:
	$(RM) *.so.o core* *.cp *.log *~ *.so 

distclean: clean
	$(RM) *.cca

