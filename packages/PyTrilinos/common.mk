# @HEADER
# ************************************************************************
#
#                  PyTrilinos: Rapid Prototyping Package
#                   Copyright (2005) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
# @HEADER

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
PYTHON_INCLUDE  := -I$(shell $(ROOT)pyLocate --include)
AUTODEP         := $(ROOT)autodep
SWIG            := swig
SWIG_FLAGS      := -Wall -noruntime -python -c++

# The wrapper files
WRAPPERS        := $(patsubst %.i, %_wrap.cxx, $(INTERFACES))

# Dynamic libraries for loading into python
DYNAMIC_LIBS    := $(patsubst %.i, _%.so, $(INTERFACES))

# Proxies and PYC files
PROXIES         := $(patsubst %.i, %.py,  $(INTERFACES))

# Dependency files
DEPDIR          := depend/
DEPEND          := $(patsubst %.i,   $(DEPDIR)%.d, $(INTERFACES)) #\
#	           $(patsubst %.cxx, $(DEPDIR)%.d, $(SRC)       )

# All the objects
#OBJECTS         := $(patsubst %.cxx,      %.o,      $(SRC)     )
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
	$(SWIG) $(COMMON_INCLUDE) $(TRILINOS_INCLUDE1) $(SWIG_FLAGS) $<

# # Generate an object file from a C++ file and its header
# %.o: %.cxx %.h
# 	$(CXX) $(CXXFLAGS) -DHAVE_CONFIG_H -I. $(PYTHON_INCLUDE) \
# 	$(TRILINOS_INCLUDE2) -c $<

# Generate a dependency file from a SWIG interface file
$(DEPDIR)%.d: %.i
	$(AUTODEP) -DHAVE_CONFIG_H -DSWIG $(TRILINOS_INCLUDE1) \
        $(PYTHON_INCLUDE) -S_wrap.cxx $< $@

# Generate a dependency file from a C++ file
$(DEPDIR)%.d: %.cxx
	$(AUTODEP) -DHAVE_CONFIG_H $(TRILINOS_INCLUDE1) $(PYTHON_INCLUDE) $< $@

# Remove emacs backup files
sweep: clobber
	rm -f *~
