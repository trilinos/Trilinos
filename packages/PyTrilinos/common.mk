####################################################
######   Common Makefile macros and targets   ######
####################################################
##                                                ##
##   The following macros should be set           ##
##   before including common.mk:                  ##
##                                                ##
##     ROOT        (root PyTrilinos directory)    ##
##     INTERFACES  (SWIG interface files, if any) ##
##     SRC         (C++ source files, if any)     ##
##     LIBRARY     (Library to be built, if any)  ##
##                                                ##
####################################################

# Get the operating system name
UNAME       := $(shell uname)

# Get the host name
HOSTNAME    := $(shell hostname)

# Get the location of the python executable
PYTHON_HOME := $(shell $(ROOT)pyLocate)

# Get the python name with version number
PYTHON_NAME := $(shell $(ROOT)pyLocate --name)

# Get the full present working directory name
PWD         := $(shell pwd)

# Get the base present working directory name
BASE_PWD    := $(shell basename $(PWD))

# Default Trilinos home
TRILINOS_HOME := /usr/local

# System-specific macros
ifeq ($(UNAME),Darwin)
  ifeq ($(HOSTNAME),samt5980.sandia.gov)
    AAL_HOME      := /Users/aalorbe
    TRILINOS_HOME := $(AAL_HOME)/local
  endif
  CXX               := c++
  CXXFLAGS          := -Wno-long-double
else
  ifeq ($(HOSTNAME),sadl12555)
    TRILINOS_HOME   := $(HOME)/scratch2/local
  endif
  ifeq ($(HOSTNAME),sahp4960)
    TRILINOS_HOME   := /usr/netpub/Trilinos-10_31_03
  endif
  ifeq ($(HOSTNAME),)
    TRILINOS_HOME   := /smallHD/scratch/install
  endif
  CXX               := g++ -g
endif

TRILINOS_INCLUDE1 := -I$(TRILINOS_HOME)/include

ifeq ($(TRILINOS_INCLUDE1),-I/usr/local/include)
  TRILINOS_INCLUDE2 :=
else
  TRILINOS_INCLUDE2 := $(TRILINOS_INCLUDE1)
endif

COMMON          := $(ROOT)src
COMMON_INCLUDE  := -I$(COMMON)
PYTHON_INCLUDE  := -I$(PYTHON_HOME)/include/$(PYTHON_NAME)
NUMERIC_INCLUDE := $(PYTHON_INCLUDE)/Numeric
AUTODEP         := $(ROOT)autodep
SWIG            := swig

# The wrapper files
WRAPPERS        := $(patsubst %.i, %_wrap.cxx, $(INTERFACES))

# Dynamic libraries for loading into python
DYNAMIC_LIBS    := $(patsubst %.i, _%.so, $(INTERFACES))

# Proxies and PYC files
PROXIES         := $(patsubst %.i, %.py,  $(INTERFACES))

# Dependency files
DEPDIR          := depend/
DEPEND          := $(patsubst %.i,   $(DEPDIR)%.d, $(INTERFACES)) \
	           $(patsubst %.cxx, $(DEPDIR)%.d, $(SRC)       )

# All the objects
OBJECTS         := $(patsubst %.cxx,      %.o,      $(SRC)     )
WRAP_OBJECTS    := $(patsubst %_wrap.cxx, %_wrap.o, $(WRAPPERS))

# Phony targets
.PHONY: sweep

# The user can add to this target in the calling Makefile,
# but this ensures it is the first target
all:

# Include the dependency files
include $(DEPEND)

# Generate a C++ wrapper and proxy file from a SWIG interface
%_wrap.cxx %.py: %.i
	$(SWIG) $(COMMON_INCLUDE) $(TRILINOS_INCLUDE1) -c -python -c++ -shadow $<

# Generate an object file from a C++ file and its header
%.o: %.cxx %.h
	$(CXX) $(CXXFLAGS) -DHAVE_CONFIG_H -I. $(PYTHON_INCLUDE) \
	$(NUMERIC_INCLUDE) $(TRILINOS_INCLUDE2) -c $<

# Generate a dependency file from a SWIG interface file
$(DEPDIR)%.d: %.i
	$(AUTODEP) -DHAVE_CONFIG_H -DSWIG $(TRILINOS_INCLUDE1) \
        $(PYTHON_INCLUDE) $(NUMERIC_INCLUDE) -S_wrap.cxx $< $@

# Generate a dependency file from a C++ file
$(DEPDIR)%.d: %.cxx
	$(AUTODEP) -DHAVE_CONFIG_H $(TRILINOS_INCLUDE1) $(PYTHON_INCLUDE) \
        $(NUMERIC_INCLUDE) $< $@

# Remove emacs backup files
sweep: clobber
	rm -f *~
