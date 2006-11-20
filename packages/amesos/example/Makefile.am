#@HEADER
# ************************************************************************
# 
#                 Amesos: Direct Sparse Solver Package
#                 Copyright (2004) Sandia Corporation
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
#@HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package amesos subdirectory example built successfully."
	@echo ""

include $(top_builddir)/Makefile.export.amesos

if USING_GNUMAKE
EXAMPLE_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(AMESOS_LIBS) $(AMESOS_EXTRA_LIBS))
EXAMPLE_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(AMESOS_INCLUDES) $(AMESOS_EXTRA_INCLUDES))
else
EXAMPLE_LIBS     = $(AMESOS_LIBS)     $(AMESOS_EXTRA_LIBS)
EXAMPLE_INCLUDES = $(AMESOS_INCLUDES) $(AMESOS_EXTRA_INCLUDES)
endif

EXEEXT = .exe

if HAVE_AMESOS_PARAKLETE
# PARAKLETE_TEST = run_pk
PARAKLETE_TEST = 
else
PARAKLETE_TEST = 
endif


if HAVE_AMESOS_GALERI
EXAMPLES_GALERI = example_AmesosFactory example_AmesosFactory_HB compare_solvers $(PARAKLETE_TEST)
else
EXAMPLES_GALERI = 
endif

if HAVE_AMESOS_STRATIMIKOS
if HAVE_AMESOS_ML
if HAVE_AMESOS_THYRA
EXAMPLES_STRATIMIKOS = stratimikos_example
else
EXAMPLES_STRATIMIKOS = 
endif
else
EXAMPLES_STRATIMIKOS = 
endif
else
EXAMPLES_STRATIMIKOS = 
endif


if HAVE_AMESOS_MC64
EXAMPLES_MC64 = example_MC64
else
EXAMPLES_MC64 = 
endif

if HAVE_MPI
EXAMPLES_MPI = a_trivial_mpi_test
else
EXAMPLES_MPI = 
endif

noinst_PROGRAMS = $(EXAMPLES_MPI) example_AmesosFactory_Tridiag $(EXAMPLES_STRATIMIKOS)  $(EXAMPLES_GALERI) $(EXAMPLES_MC64)

EXTRA_DIST = In_662_bus.rsa _stratimikos_amesos_lapack.xml

stratimikos_example_SOURCES      = $(srcdir)/stratimikos_example.cpp  $(srcdir)/simpleStratimikosSolve.cpp $(srcdir)/Thyra_AmesosLinearOpWithSolveFactory.cpp
stratimikos_example_DEPENDENCIES = $(top_builddir)/src/libamesos.a $(top_builddir)/../stratimikos/src/libstratimikos.a 
stratimikos_example_CXXFLAGS     = $(EXAMPLE_INCLUDES)  -I$(top_builddir)/../stratimikos/src/  -I$(top_builddir)/../thyra/src/  -I$(top_srcdir)/../stratimikos/src/ -I$(top_srcdir)/../thyra/src/
stratimikos_example_LDADD        =  $(top_builddir)/../stratimikos/src/libstratimikos.a  $(top_builddir)/../aztecoo/thyra/src/libaztecoothyra.a  $(top_builddir)/../aztecoo/src/libaztecoo.a  $(top_builddir)/../ml/thyra/src/libmlthyra.a  $(top_builddir)/../ml/src/libml.a  $(STRATIMIKOS_LIBS) $(EXAMPLE_LIBS)

example_AmesosFactory_SOURCES      = $(srcdir)/example_AmesosFactory.cpp
example_AmesosFactory_DEPENDENCIES = $(top_builddir)/src/libamesos.a 
example_AmesosFactory_CXXFLAGS     = $(EXAMPLE_INCLUDES)
example_AmesosFactory_LDADD        = $(EXAMPLE_LIBS)

example_AmesosFactory_Tridiag_SOURCES      = $(srcdir)/example_AmesosFactory_Tridiag.cpp
example_AmesosFactory_Tridiag_DEPENDENCIES = $(top_builddir)/src/libamesos.a
example_AmesosFactory_Tridiag_CXXFLAGS     = $(EXAMPLE_INCLUDES)
example_AmesosFactory_Tridiag_LDADD        = $(EXAMPLE_LIBS)

example_AmesosFactory_HB_SOURCES      = $(srcdir)/example_AmesosFactory_HB.cpp 
example_AmesosFactory_HB_DEPENDENCIES = $(top_builddir)/src/libamesos.a 
example_AmesosFactory_HB_CXXFLAGS     = $(EXAMPLE_INCLUDES)
example_AmesosFactory_HB_LDADD        = $(EXAMPLE_LIBS)

compare_solvers_SOURCES      = $(srcdir)/compare_solvers.cpp
compare_solvers_DEPENDENCIES = $(top_builddir)/src/libamesos.a 
compare_solvers_CXXFLAGS     = $(EXAMPLE_INCLUDES)
compare_solvers_LDADD        = $(EXAMPLE_LIBS)

example_MC64_SOURCES      = $(srcdir)/example_MC64.cpp
example_MC64_DEPENDENCIES = $(top_builddir)/src/libamesos.a 
example_MC64_CXXFLAGS     = $(EXAMPLE_INCLUDES)
example_MC64_LDADD        = $(EXAMPLE_LIBS)

# run_p[k_SOURCES      = $(srcdir)/run_pk.c
# run_p[k_DEPENDENCIES = $(top_builddir)/src/libamesos.a 
# run_p[k_CXXFLAGS     = $(EXAMPLE_INCLUDES)
# run_p[k_LDADD        = $(EXAMPLE_LIBS)

a_trivial_mpi_test_SOURCES      = $(srcdir)/a_trivial_mpi_test.cpp
a_trivial_mpi_test_CXXFLAGS     = 
a_trivial_mpi_test_LDADD        = 

clean:  clean-am
	rm -f  662_bus.rsa stratimikos_amesos_lapack.xml

