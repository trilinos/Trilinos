#@HEADER
# ************************************************************************
# 
#                     Stokhos Package
#                 Copyright (2008) Sandia Corporation
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
# Questions? Contact Eric T. Phipps (etphipp@sandia.gov) 
# 
# ************************************************************************
#@HEADER

# Include export Makefiles for correct paths
include $(top_builddir)/Makefile.export.stokhos

# The following line helps the test harness recover from build errors.

all-local: generate-sample-output.pl

generate-sample-output.pl : generate-sample-output.pl.stub
	cp generate-sample-output.pl.stub generate-sample-output.pl
	chmod a+x generate-sample-output.pl

EXEEXT = .exe

if ENABLE_STOKHOS_QUADRATURE

QUAD_EXAMPLES = product_basis_example exp_example

endif

noinst_PROGRAMS = hermite_example dakota_pce dakota_exp $(QUAD_EXAMPLES)

AM_CPPFLAGS = $(STOKHOS_INCLUDES)

hermite_example_SOURCES = $(srcdir)/hermite_example.cpp
hermite_example_DEPENDENCIES = $(top_builddir)/src/libstokhos.a 
hermite_example_LDADD = $(STOKHOS_LIBS)

product_basis_example_SOURCES = $(srcdir)/product_basis_example.cpp
product_basis_example_DEPENDENCIES = $(hermite_example_DEPENDENCIES)
product_basis_example_LDADD = $(hermite_example_LDADD)

dakota_pce_SOURCES = $(srcdir)/dakota_pce.cpp
dakota_pce_DEPENDENCIES = $(hermite_example_DEPENDENCIES)
dakota_pce_LDADD = $(hermite_example_LDADD)

exp_example_SOURCES = $(srcdir)/exp_example.cpp
exp_example_DEPENDENCIES = $(hermite_example_DEPENDENCIES)
exp_example_LDADD = $(hermite_example_LDADD)

dakota_exp_SOURCES = $(srcdir)/dakota_exp.cpp
dakota_exp_DEPENDENCIES = $(hermite_example_DEPENDENCIES)
dakota_exp_LDADD = $(hermite_example_LDADD)

SUBDIRS = 
