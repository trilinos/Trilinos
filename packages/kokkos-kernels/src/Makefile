#=======================================================================
#=================== Settings ==========================================
#=======================================================================

ifndef KOKKOSKERNELS_PATH
  MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
  KOKKOSKERNELS_PATH = $(subst Makefile.minimal,,$(MAKEFILE_PATH))..
endif

KOKKOS_PATH = ${KOKKOSKERNELS_PATH}/../kokkos


KOKKOS_DEVICES=OpenMP

KOKKOSKERNELS_SCALARS="'double,complex<double>'"
KOKKOSKERNELS_OPTIONS=eti-only

ifndef PREFIX
  PREFIX = ${CURDIR}/install
endif
KOKKOSKERNELS_INSTALL_PATH = ${PREFIX}

default: build

CXX = g++

CXXFLAGS = -O3 -g
LINK = ${CXX}
LINKFLAGS =

DEPFLAGS = -M

OBJ =
LIB = -lpthread


#=======================================================================
#=================== Include Kokkos and KokkosKernels ==================
#=======================================================================

include $(KOKKOSKERNELS_PATH)/Makefile.kokkos-kernels

kokkoskernels-install: kokkoskernels-build-lib
	mkdir -p $(KOKKOSKERNELS_INSTALL_PATH)
	mkdir -p $(KOKKOSKERNELS_INSTALL_PATH)/include
	mkdir -p $(KOKKOSKERNELS_INSTALL_PATH)/lib
	echo "# KokkosKernels Settings" > $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOSKERNELS_PATH = ${KOKKOSKERNELS_INSTALL_PATH}" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOSKERNELS_CPPFLAGS = -I$(KOKKOSKERNELS_INSTALL_PATH)/include $(KOKKOSKERNELS_CPPFLAGS)" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOSKERNELS_LDFLAGS = -L$(KOKKOSKERNELS_INSTALL_PATH)/lib $(KOKKOSKERNELS_LDFLAGS)" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOSKERNELS_LIBS = $(KOKKOSKERNELS_LIBS)" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOSKERNELS_LINK_DEPENDS = $(KOKKOSKERNELS_INSTALL_PATH)/lib/$(KOKKOSKERNELS_LINK_DEPENDS)" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "# KokkosKernels Internal Settings" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOSKERNELS_INTERNAL_ETI_ONLY = ${KOKKOSKERNELS_INTERNAL_ETI_ONLY}" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOSKERNELS_INTERNAL_INST_EXECSPACE_CUDA = ${KOKKOSKERNELS_INTERNAL_INST_EXECSPACE_CUDA}" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOSKERNELS_INTERNAL_INST_EXECSPACE_OPENMP = ${KOKKOSKERNELS_INTERNAL_INST_EXECSPACE_OPENMP}" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOSKERNELS_INTERNAL_INST_EXECSPACE_SERIAL = ${KOKKOSKERNELS_INTERNAL_INST_EXECSPACE_SERIAL}" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOSKERNELS_INTERNAL_INST_EXECSPACE_THREADS = ${KOKKOSKERNELS_INTERNAL_INST_EXECSPACE_THREADS}" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOSKERNELS_INTERNAL_INST_EXECSPACE_ROCM = ${KOKKOSKERNELS_INTERNAL_INST_EXECSPACE_ROCM}" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "# Kokkos Settings Used for KokkosKernels Compilation" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOS_PATH ?= ${KOKKOS_PATH}" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOS_OS ?= ${KOKKOS_OS}" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOS_CXXFLAGS ?= $(KOKKOS_CXXFLAGS)" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOS_CPPFLAGS ?= $(KOKKOS_CPPFLAGS)" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOS_LDFLAGS ?= $(KOKKOS_LDFLAGS)" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOS_LIBS ?= $(KOKKOS_LIBS)" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOS_LINK_DEPENDS ?= $(KOKKOS_LINK_DEPENDS)" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "# Kokkos Internal Settings" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOS_INTERNAL_USE_CUDA = ${KOKKOS_INTERNAL_USE_CUDA}" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOS_INTERNAL_USE_OPENMP = ${KOKKOS_INTERNAL_USE_OPENMP}" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOS_INTERNAL_USE_SERIAL = ${KOKKOS_INTERNAL_USE_SERIAL}" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOS_INTERNAL_USE_PTHREADS = ${KOKKOS_INTERNAL_USE_PTHREADS}" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "KOKKOS_INTERNAL_USE_ROCM = ${KOKKOS_INTERNAL_USE_ROCM}" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "# Fake target" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	echo "libkokkos_kernels.a:" >> $(KOKKOSKERNELS_INSTALL_PATH)/Makefile.kokkos-kernels
	$(CP) $(COPY_FLAG) KokkosKernels_config.h $(KOKKOSKERNELS_INSTALL_PATH)/include
	$(CP) $(COPY_FLAG) ${KOKKOSKERNELS_PATH}/src/*.hpp $(KOKKOSKERNELS_INSTALL_PATH)/include
	$(CP) $(COPY_FLAG) ${KOKKOSKERNELS_PATH}/src/impl/*.hpp $(KOKKOSKERNELS_INSTALL_PATH)/include
	$(CP) $(COPY_FLAG) ${KOKKOSKERNELS_PATH}/src/blas/*.hpp $(KOKKOSKERNELS_INSTALL_PATH)/include
	$(CP) $(COPY_FLAG) ${KOKKOSKERNELS_PATH}/src/blas/impl/*.hpp $(KOKKOSKERNELS_INSTALL_PATH)/include
	$(CP) $(COPY_FLAG) ${KOKKOSKERNELS_PATH}/src/batched/*.hpp $(KOKKOSKERNELS_INSTALL_PATH)/include
	$(CP) $(COPY_FLAG) ${KOKKOSKERNELS_PATH}/src/sparse/*.hpp $(KOKKOSKERNELS_INSTALL_PATH)/include
	$(CP) $(COPY_FLAG) ${KOKKOSKERNELS_PATH}/src/sparse/impl/*.hpp $(KOKKOSKERNELS_INSTALL_PATH)/include
	$(CP) $(COPY_FLAG) ${KOKKOSKERNELS_PATH}/src/graph/*.hpp $(KOKKOSKERNELS_INSTALL_PATH)/include
	$(CP) $(COPY_FLAG) ${KOKKOSKERNELS_PATH}/src/graph/impl/*.hpp $(KOKKOSKERNELS_INSTALL_PATH)/include
	$(CP) $(COPY_FLAG) ${KOKKOSKERNELS_PATH}/src/common/*.hpp $(KOKKOSKERNELS_INSTALL_PATH)/include
	$(CP) $(COPY_FLAG) -r ${KOKKOSKERNELS_PATH}/src/impl/generated_specializations_hpp $(KOKKOSKERNELS_INSTALL_PATH)/include
	$(CP) $(COPY_FLAG) ${KOKKOSKERNELS_PATH}/src/impl/tpls/*.hpp $(KOKKOSKERNELS_INSTALL_PATH)/include
	$(CP) $(COPY_FLAG) $(KOKKOSKERNELS_INTERNAL_LIBRARY) $(KOKKOSKERNELS_INSTALL_PATH)/lib

build: kokkoskernels-build-lib

install: kokkoskernels-install
 
clean: kokkoskernels-clean kokkos-clean
	rm -f *.o
#	rm -f *.o $(KOKKOSKERNELS_INTERNAL_LIBRARY) KokkosKernels_config.h KokkosKernels_config.tmp KokkosCore_config.h KokkosCore_config.tmp 

