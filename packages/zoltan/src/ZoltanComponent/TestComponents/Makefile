## Makefile for the Driver and Mesh components. These are wrappers around the 
## code in $(ZOLTAN_HOME)/driver and $(ZOLTAN_HOME)/ch. The user should not
## need to modify this file; all the user-configurable stuff need
## comes in from MakeIncl.ZoltanComponent

# Suck in the user-configuarble stuff
include ../MakeIncl.ZoltanComponent

## suck in some CCA specific compile flags
include $(CCA_MAKEINCL_DIR)/MakeIncl.CCA_Component

.SUFFIXES: .c .cpp .o .so.o .a

SAPPOBJ_Driver =  Driver.so.o  Driver_wrapper.so.o

SAPPOBJ_Mesh =  Mesh.so.o  Mesh_info.so.o  Mesh_wrapper.so.o Mesh_DataMig.so.o

ZOLTAN_BIN_DIR = $(ZOLTAN_HOME)/Obj_$(ZOLTAN_ARCH)

DR_OBJ = $(ZOLTAN_BIN_DIR)/dr_err.o \
         $(ZOLTAN_BIN_DIR)/dr_input.o \
         $(ZOLTAN_BIN_DIR)/dr_par_util.o \
         $(ZOLTAN_BIN_DIR)/dr_util.o \
         $(ZOLTAN_BIN_DIR)/dr_exoII_io.o \
         $(ZOLTAN_BIN_DIR)/dr_elem.o \
         $(ZOLTAN_BIN_DIR)/dr_loadbal.o \
         $(ZOLTAN_BIN_DIR)/dr_migrate.o \
         $(ZOLTAN_BIN_DIR)/dr_output.o \
         $(ZOLTAN_BIN_DIR)/dr_chaco_io.o \
         $(ZOLTAN_BIN_DIR)/dr_elem_util.o \
         $(ZOLTAN_BIN_DIR)/dr_maps.o \
         $(ZOLTAN_BIN_DIR)/dr_eval.o \
         $(ZOLTAN_BIN_DIR)/dr_gnuplot.o \
         $(ZOLTAN_BIN_DIR)/dr_random_io.o

CH_OBJ = $(ZOLTAN_BIN_DIR)/ch_dist_graph.o \
         $(ZOLTAN_BIN_DIR)/ch_input_assign.o \
         $(ZOLTAN_BIN_DIR)/ch_input_geom.o \
         $(ZOLTAN_BIN_DIR)/ch_input_graph.o \
         $(ZOLTAN_BIN_DIR)/ch_input_read.o \
         $(ZOLTAN_BIN_DIR)/ch_init_dist.o

libshared_driver = Driver$(CCA_SHLIB_SUFFIX)
libshared_mesh = Mesh$(CCA_SHLIB_SUFFIX)

# Application specific Flags
C++APPFLAGS = -Wno-unused  -DNDEBUG

# Do we need to link in ParMetis and Jostle ?
ifdef PARMETIS_LIBPATH
PARMETIS_LIBS = -lparmetis -lmetis
endif

ifdef JOSTLE_LIBPATH
JOSTLE_LIBS = -ljostle
endif

ifdef NEMESIS_LIBPATH
NEMESIS_LIBS = -lnemIc -lexoIIv2c -lnetcdf
endif

# Application specific libraries like blas, lapack etc.
APPLIB = -L. -lextras -L$(ZOLTAN_HOME)/Obj_$(ZOLTAN_ARCH) -lzoltan  \
	 -L$(DCCAFE_LIB)/ -Wl,-rpath,$(DCCAFE_LIB) -l$(DCCAFE_LIBNAME) \
	 $(PARMETIS_LIBPATH) $(PARMETIS_LIBS) $(JOSTLE_LIBPATH) $(JOSTLE_LIBS)

LDLIBS = 
#-lstdc++

# Zoltan specific include
INCLUDE1 = -I$(ZOLTAN_HOME) -I$(ZOLTAN_HOME)/include \
            -I$(ZOLTAN_HOME)/Utilities/Communication -I$(ZOLTAN_HOME)/Utilities/Memory 
INCLUDE2 = -I$(ZOLTAN_HOME)/Utilities/DDirectory -I$(ZOLTAN_HOME)/driver

# MPI and CCA specific includes
INCLUDE3 = $(CCA_INC)
INCLUDE4 = -I $(ZOLTAN_PORTS_DIR)
INCLUDE5 = -I$(DCCAFE_INCLUDE)
INCLUDE6 = -I$(MPI_HOME)/include 
INCLUDE = $(INCLUDE1) $(INCLUDE2) $(INCLUDE3) $(INCLUDE4) $(INCLUDE5) $(INCLUDE6)

# Misc. definition
ECHO = echo
RM = /bin/rm -f

default : driver mesh

driver : $(SAPPOBJ_Driver) libextras.a
	$(RM) $(libshared_driver)
	$(CCA_SHLIB_LD) $(CCA_CXX_FLAGS) $(CCA_SHLIB_LDFLAGS) -o $(libshared_driver)  \
	$(SAPPOBJ_Driver) $(APPLIB) $(LDLIBS)

mesh : $(SAPPOBJ_Mesh) libextras.a
	$(RM) $(libshared_mesh)
	$(CCA_SHLIB_LD) $(CCA_CXX_FLAGS) $(CCA_SHLIB_LDFLAGS) -o $(libshared_mesh)  \
	$(SAPPOBJ_Mesh) $(APPLIB) $(LDLIBS) $(NEMESIS_LIBPATH) $(NEMESIS_LIBS)

.cpp.so.o : 
	@$(ECHO)
	@$(ECHO) "Compiling Source File --" $<
	@$(ECHO) "---------------------"
	$(CCA_CXX) $(CCA_CXX_FLAGS) $(CCA_SHLIB_CFLAGS) $(C++APPFLAGS)  -c $(INCLUDE) -o $@ $< 
	@$(ECHO)

libextras.a : $(DR_OBJ) $(CH_OBJ)
	ar rv libextras.a $(DR_OBJ) $(CH_OBJ)
	ranlib libextras.a

clean:
	$(RM) *.so.o core* *.cp *.log *~ *.so *.a

distclean:
	$(RM) *.cca
