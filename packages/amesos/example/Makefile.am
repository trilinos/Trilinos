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

include $(top_builddir)/Makefile.export.amesos

if USING_GNUMAKE
EXAMPLE_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(AMESOS_LIBS) $(AMESOS_EXTRA_LIBS))
EXAMPLE_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(AMESOS_INCLUDES) $(AMESOS_EXTRA_INCLUDES))
else
EXAMPLE_LIBS     = $(AMESOS_LIBS)     $(AMESOS_EXTRA_LIBS)
EXAMPLE_INCLUDES = $(AMESOS_INCLUDES) $(AMESOS_EXTRA_INCLUDES)
endif

AM_CPPFLAGS = $(EXAMPLE_INCLUDES)

EXEEXT = .exe

noinst_PROGRAMS = example_AmesosFactory example_AmesosFactory_HB example_AmesosFactory_Tridiag compare_solvers

EXTRA_DIST = In_662_bus.rsa

#copy these files over to the build area
GetNecessaryFiles: $(EXTRA_DIST)
	cp -f $(srcdir)/In_662_bus.rsa 662_bus.rsa
	chmod -x 662_bus.rsa

example_AmesosFactory_SOURCES      = $(srcdir)/example_AmesosFactory.cpp
example_AmesosFactory_DEPENDENCIES = $(top_builddir)/src/libamesos.a 
example_AmesosFactory_LDADD        = $(EXAMPLE_LIBS)

example_AmesosFactory_Tridiag_SOURCES      = $(srcdir)/example_AmesosFactory_Tridiag.cpp
example_AmesosFactory_Tridiag_DEPENDENCIES = $(top_builddir)/src/libamesos.a
example_AmesosFactory_Tridiag_LDADD        = $(EXAMPLE_LIBS)

example_AmesosFactory_HB_SOURCES      = $(srcdir)/example_AmesosFactory_HB.cpp 
example_AmesosFactory_HB_DEPENDENCIES = $(top_builddir)/src/libamesos.a GetNecessaryFiles
example_AmesosFactory_HB_LDADD        = $(EXAMPLE_LIBS)

compare_solvers_SOURCES      = $(srcdir)/compare_solvers.cpp
compare_solvers_DEPENDENCIES = $(top_builddir)/src/libamesos.a 
compare_solvers_LDADD        = $(EXAMPLE_LIBS)

clean:  clean-am
	rm -f  662_bus.rsa
