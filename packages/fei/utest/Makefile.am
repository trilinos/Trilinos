#@HEADER
# ************************************************************************
# 
#                               FEI
#                 Copyright (2005) Sandia Corporation
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
# Questions? Contact Alan B. Williams (william@sandia.gov) 
# 
# ************************************************************************
#@HEADER

include $(top_builddir)/Makefile.export.fei

FETIDP_INCLUDE=@FETIDP_INCLUDE@
FETIDP_LIB_DIR=@FETIDP_LIB_DIR@

DEFAULT_INCLUDES = $(FEI_INCLUDES) \
	-I$(top_srcdir)/test_utils -I$(top_srcdir)/utest_cases \
	$(FETIDP_INCLUDE)

INSTRUMENT_CMD=@INSTRUMENT_CMD@

if HAVE_FEI_FETI
FEI_FETIDP_LIBS=-L$(FETIDP_LIB_DIR) -lfeti-mpi
else
FEI_FETIDP_LIBS=
endif

FEI_TEST=fei_unit_tests

SPEC_RUN_COMMAND=@SPECIFIED_RUN_COMMAND@

if RUN_COMMAND_SPECIFIED
RUN_COMMAND=$(SPEC_RUN_COMMAND)
else
if FEI_SERIAL
RUN_COMMAND=" "
else
RUN_COMMAND=mpirun -np
endif
endif


DEFS = 
CXXLD=$(INSTRUMENT_CMD) $(CXX)

EXEEXT = .exe

noinst_PROGRAMS = $(FEI_TEST)

# One SOURCES primary is needed for EACH exe.  (Notice for the following 
# primaries that the pattern is "(Exe_name)_(PRIMARY_NAME) = ").  List all source 
# files needed to compile the exe.
fei_unit_tests_SOURCES = $(srcdir)/fei_unit_tests.cpp

# One DEPENDENCIES primary is needed for each exe.  List any libraries that the 
# exe depends on.  (This primary is not absolutely essential.)
fei_unit_tests_DEPENDENCIES = \
	$(FEI_LIBS) \
	$(top_builddir)/utest_cases/libfei_utest_cases.a \
	$(top_builddir)/base/libfei_base.a \
	$(top_builddir)/support-Trilinos/libfei_trilinos.a \
	$(top_builddir)/test_utils/libfei_test_utils.a

# One CXXFLAGS primary can be defined for each exe.  Include any -I's necessary for 
# building the exe.
fei_unit_tests_CXXFLAGS = $(FEI_INCLUDES)

# One LDADD primary can be defined for each exe.
# Include any LDFLAGS necessary for building the exe.
fei_unit_tests_LDADD = \
	-L$(top_builddir)/utest_cases -lfei_utest_cases \
	-L$(top_builddir)/test_utils -lfei_test_utils \
	$(FEI_LIBS) \
	$(FEI_FETIDP_LIBS) \
	$(LIBS) $(FLIBS)

