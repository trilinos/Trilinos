### NOTE: This makefile requires GNU make.
### Only the GNU COMPILER has been tested lately with this Makefile,
### However, the code has been compiled successfully on all normal Sierra platforms.

### EXTERNAL LIBRARIES (TPLS) NEEDED:
### exodusII: libexoIIv2c.a
### nemesis:  libnemIc.a
### netcdf:   libnetcdf.a
###
### By default, this Makefile looks for these libraries and associated include files
### in $(ACCESS)/lib and $(ACCESS)/inc.

EXOINC = -I${ACCESS}/inc 
#EXOLIB = -Wl,-rpath ${ACCESS}/lib -L${ACCESS}/lib -lexoIIv2c -lnemIc -lnetcdf 
EXOLIB = -Wl,-rpath ${ACCESS}/lib -L${ACCESS}/lib -lexoIIv2c -lnemIc -lnetcdf -lhdf5_hl -lhdf5 -lm -lz

### If compiling XDMF support, then need HDF5 libraries.
### Change the HDFINC, HDFLIB as needed. They can
### either be defined globally or in each systems block. 
### Note that this (HDF) is only needed if compiling in support for XDMF

#HDFINC = -I/usr/local/include
#HDFLIB = -L/usr/local/lib -lhdf5


COMPILER = GNU
DEBUG = -g 
#COMPILER = Intel
#COMPILER=Parallel

AR = ar
#ARFLAGS = cr
# Purple needs:  (also setenv OBJECT_MODE 64)
#ARFLAGS = -X64 srv

ifneq ($(COMPILER),Parallel)
DEBUG += -DNO_MPI
endif

MAKEDEPEND = makedepend
FLAGS = -DNO_XDMF_SUPPORT

SYSTEM = $(shell uname)

### Update to correct path for your platform.


### Update to correct path for your platform.
### Define here or within a system block...

#========================================================================
#       Linux
#========================================================================
ifeq ($(SYSTEM),Linux)

ifeq ($(COMPILER),Portland)
CC       = pgcc
CXX      = pgCC
CXXFLAGS  = -Mframe -DNO_CHEADERS --prelink_objects --instantiation_dir $(shell pwd)/bin/$(SYSTEM)/template.dir
endif

ifeq ($(COMPILER),Intel)
CC       = icc
CXX      = icpc 
CXXFLAGS  = -Xc 
endif

ifeq ($(COMPILER),GNU)
CC	= gcc
CXX	= g++
CXXFLAGS = -Wall -Wno-deprecated 
endif

ifeq ($(COMPILER),Parallel)
# These are for a 32-bit build using gcc
EXOINC = -I/usr/local/eng_sci/struct/i686/current-gcc/inc
EXOLIB = -L/usr/local/eng_sci/struct/i686/current-gcc/lib -lexoIIv2c -lnemIc -lnetcdf

# This path is for scico linux lan...
PATH := /usr/local/mpi/sierra/mpich/1.2.5.2/gcc-3.2.2/bin:/usr/bin:$(PATH)
CC	= mpicc
CXX	= mpiCC
CXXFLAGS = -Wno-deprecated
endif
endif

#========================================================================
#       SUN
#========================================================================
ifeq ($(SYSTEM),SunOS)
ifeq ($(COMPILER),GNU)
CC	= gcc
CXX	= g++
CXXFLAGS = -Wall
else
CC        = cc
CXX       = CC
CXXFLAGS  = -DNO_MEMBER_TEMPLATES -ptr$(shell pwd)/SunOS/SunWS_cache
endif
endif

#========================================================================
#       SGI
#========================================================================
ifeq ($(SYSTEM),IRIX64)
CC        = cc
CXX       = CC
CXXFLAGS  = -64 -LANG:std -DNO_CHEADERS
endif

#========================================================================
#       IBM
#========================================================================
ifeq ($(SYSTEM),AIX)
CC        = xlc
CXX       = /usr/vacpp/bin/xlC -q64
CXXFLAGS  = 
SYSLIB	  = -lm
endif

#========================================================================

#========================================================================
#       OSF (DEC Alpha)
#========================================================================
ifeq ($(SYSTEM),OSF1)
CC        = cc
CXX       = cxx
CXXFLAGS  = -std strict_ansi -ptr $(shell pwd)/OSF1/cxx_repository -DNO_CHEADERS -noimplicit_include -ttimestamp 
MAKEDEPEND = $(CXX) -Em
endif

### Shouldn't need to touch these...
vpath %.o ./bin/$(SYSTEM)
vpath %.C ./io_system/src
vpath %.C ./io_system/src/
vpath %.C ./io_system/src/init
vpath %.C ./io_system/src/exodusII
vpath %.C ./io_system/src/transform
vpath %.C ./io_system/src/heartbeat
vpath %.C ./io_system/src/generated
vpath %.C ./io_system/src/main

.SUFFIXES: .C .o

.C.o:
	$(CXX) $(DEBUG) $(FLAGS) ${CXXFLAGS} -c $(INCLUDES) $(EXOINC) $(HDFINC) -o ./bin/$(SYSTEM)/$@ $<


IOBJ  = io_info.o
MOBJ  = Main_io_shell.o
COBJ  = cth_pressure_map.o
SOBJ  = shell_to_hex.o
OBJS  =	Ionit_Initializer.o \
	Ioss_CommSet.o \
	Ioss_CompositeVariableType.o \
	Ioss_ConcreteVariableType.o \
	Ioss_ConstructedVariableType.o \
	Ioss_DatabaseIO.o \
	Ioss_EdgeBlock.o \
	Ioss_EdgeSet.o \
	Ioss_ElementBlock.o \
	Ioss_EntityBlock.o \
	Ioss_FaceBlock.o \
	Ioss_FaceSet.o \
	Ioss_Field.o \
	Ioss_FieldManager.o \
	Ioss_FileInfo.o \
	Ioss_GroupingEntity.o \
	Ioss_IOFactory.o \
	Ioss_Map.o \
	Ioss_NodeBlock.o \
	Ioss_NodeSet.o \
	Ioss_Property.o \
	Ioss_PropertyManager.o \
	Ioss_Region.o \
	Ioss_SerializeIO.o \
	Ioss_Utils.o \
	Ioss_ParallelUtils.o \
	Ioss_VariableType.o \
	Ioss_ElementTopology.o \
	Ioex_DatabaseIO.o \
	Ioex_IOFactory.o \
	Ioex_Internals.o \
	Ioex_SuperElement.o \
	Ioss_Bar2.o \
	Ioss_Bar3.o \
	Ioss_Edge2.o \
	Ioss_Edge2D2.o \
	Ioss_Edge2D3.o \
	Ioss_Edge3.o \
	Ioss_Hex20.o \
	Ioss_Hex27.o \
	Ioss_Hex8.o \
	Ioss_Initializer.o \
	Ioss_Node.o \
	Ioss_Pyramid13.o \
	Ioss_Pyramid5.o \
	Ioss_Quad4.o \
	Ioss_Quad8.o \
	Ioss_Quad9.o \
	Ioss_Shell4.o \
	Ioss_Shell8.o \
	Ioss_Shell9.o \
	Ioss_ShellLine2D2.o \
	Ioss_ShellLine2D3.o \
	Ioss_Sphere.o \
	Ioss_Super.o \
	Ioss_Tet10.o \
	Ioss_Tet4.o \
	Ioss_Tet8.o \
	Ioss_Tri3.o \
	Ioss_Tri4.o \
	Ioss_Tri6.o \
	Ioss_TriShell3.o \
	Ioss_TriShell6.o \
	Ioss_Unknown.o \
	Ioss_Wedge15.o \
	Ioss_Wedge6.o \
	Ioss_Transform.o \
	Iotr_Factory.o \
	Iotr_Initializer.o \
	Iotr_MinMax.o \
	Iotr_Offset.o \
	Iotr_Scale.o \
	Iotr_Offset3D.o \
	Iotr_Scale3D.o \
	Iotr_Tensor.o \
	Iotr_VectorMagnitude.o \
	Iogn_DatabaseIO.o \
	Iogn_GeneratedMesh.o \
	Iohb_DatabaseIO.o \
	Iohb_Layout.o

SRCS = $(OBJS:.o=.C) $(MOBJ:.o=.C)

# Directories to search for header files
INCLUDES = -I$(shell pwd)/io_system/include 

#========================================================================
LDFLAGS = 
LINKER = $(CXX)

DIR = ./bin/$(SYSTEM)

all: io_info io_shell cth_pressure_map shell_to_hex $(DIR)/libIoss.a

$(DIR)/libIoss.a: $(DIR) $(OBJS)
	cd $(DIR); $(AR) $(ARFLAGS) libIoss.a $(OBJS)

io_shell: $(DIR) $(MOBJ) $(DIR)/libIoss.a
	cd ./bin/$(SYSTEM); ${PREFIX} ${LINKER} ${CXXFLAGS} ${LDFLAGS} -o $@ $(MOBJ) libIoss.a $(EXOLIB) $(SYSLIB) 

io_info: $(DIR) $(IOBJ) $(DIR)/libIoss.a
	cd ./bin/$(SYSTEM); ${PREFIX} ${LINKER} ${CXXFLAGS} ${LDFLAGS} -o $@ $(IOBJ) libIoss.a $(EXOLIB) $(SYSLIB) 

cth_pressure_map: $(DIR) $(COBJ) $(DIR)/libIoss.a
	cd ./bin/$(SYSTEM); ${PREFIX} ${LINKER} ${CXXFLAGS} ${LDFLAGS} -o $@ $(COBJ) libIoss.a $(EXOLIB) $(SYSLIB) 

shell_to_hex: $(DIR) $(SOBJ) $(DIR)/libIoss.a
	cd ./bin/$(SYSTEM); ${PREFIX} ${LINKER} ${CXXFLAGS} ${LDFLAGS} -o $@ $(SOBJ) libIoss.a $(EXOLIB) $(SYSLIB) 

$(DIR):
	mkdir -p $(DIR)

clean:
	@${RM} core $(OBJS)  *~ io_shell *.ti *.ii
	@cd bin/$(SYSTEM); $(RM) libIoss.a $(MOBJ) $(OBJS) io_shell cth_pressure_map

really_clean:
	@$(RM) -rf bin

depend:
	touch depend.$(SYSTEM)
	$(MAKEDEPEND) -f depend.$(SYSTEM) $(CXXFLAGS) $(INCLUDES) $(SRCS)

# DO NOT DELETE
-include depend.$(SYSTEM)
