#@HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
# 
#            LOCA: Library of Continuation Algorithms Package
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
# 
# Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
# Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
# ************************************************************************
#@HEADER

include $(top_builddir)/Makefile.export.nox

if USING_GNUMAKE
TEST_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(LOCA_EPETRA_TEST_LIBS))
TEST_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(LOCA_EPETRA_TEST_INCLUDES))
else
TEST_LIBS     = $(LOCA_EPETRA_TEST_LIBS)
TEST_INCLUDES = $(LOCA_EPETRA_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(TEST_INCLUDES)

EXEEXT = .exe

# The XYZT tests require epetraext and mpi
if BUILD_NOX_EPETRAEXT_SUPPORT
BRUSSXYZT = brussXYZT_ifpack brussXYZT_global brussXYZT_sequential brussXYZT_sequential_IPS brussXYZT_sequential_OPS brussXYZT_parallel brussXYZT_blockdiagonal brussXYZT_bdsdt # brussXYZT_parareal
else
BRUSSXYZT = 
endif

if BUILD_LOCA_EPETRA_TESTS
noinst_PROGRAMS = $(BRUSSXYZT)
endif

brussXYZT_ifpack_SOURCES = \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Brusselator.C $(srcdir)/Brusselator.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H \
	$(srcdir)/brussXYZT_ifpack.C

brussXYZT_ifpack_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-epetra/libnoxepetra.a \
        $(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-epetra/liblocaepetra.a

brussXYZT_ifpack_LDADD = $(TEST_LIBS)


brussXYZT_global_SOURCES = \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Brusselator.C $(srcdir)/Brusselator.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H \
	$(srcdir)/brussXYZT_global.C

brussXYZT_global_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-epetra/libnoxepetra.a \
        $(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-epetra/liblocaepetra.a

brussXYZT_global_LDADD = $(TEST_LIBS)


brussXYZT_sequential_SOURCES = \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Brusselator.C $(srcdir)/Brusselator.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H \
	$(srcdir)/brussXYZT_sequential.C

brussXYZT_sequential_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-epetra/libnoxepetra.a \
        $(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-epetra/liblocaepetra.a

brussXYZT_sequential_LDADD = $(TEST_LIBS)


brussXYZT_sequential_IPS_SOURCES = \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Brusselator.C $(srcdir)/Brusselator.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H \
	$(srcdir)/brussXYZT_sequential_IPS.C

brussXYZT_sequential_IPS_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-epetra/libnoxepetra.a \
        $(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-epetra/liblocaepetra.a

brussXYZT_sequential_IPS_LDADD = $(TEST_LIBS)


brussXYZT_sequential_OPS_SOURCES = \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Brusselator.C $(srcdir)/Brusselator.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H \
	$(srcdir)/brussXYZT_sequential_OPS.C

brussXYZT_sequential_OPS_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-epetra/libnoxepetra.a \
        $(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-epetra/liblocaepetra.a

brussXYZT_sequential_OPS_LDADD = $(TEST_LIBS)


brussXYZT_parallel_SOURCES = \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Brusselator.C $(srcdir)/Brusselator.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H \
	$(srcdir)/brussXYZT_parallel.C

brussXYZT_parallel_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-epetra/libnoxepetra.a \
        $(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-epetra/liblocaepetra.a

brussXYZT_parallel_LDADD = $(TEST_LIBS)


# Currently breaks mpirun
#brussXYZT_parareal_SOURCES = \
#	$(srcdir)/Basis.C $(srcdir)/Basis.H \
#	$(srcdir)/Brusselator.C $(srcdir)/Brusselator.H \
#	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H \
#	$(srcdir)/brussXYZT_parareal.C
#
#brussXYZT_parareal_DEPENDENCIES = $(top_builddir)/src/libnox.a \
#        $(top_builddir)/src-epetra/libnoxepetra.a \
#        $(top_builddir)/src-loca/src/libloca.a \
#        $(top_builddir)/src-loca/src-epetra/liblocaepetra.a
#
#brussXYZT_parareal_LDADD = $(TEST_LIBS)


brussXYZT_blockdiagonal_SOURCES = \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Brusselator.C $(srcdir)/Brusselator.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H \
	$(srcdir)/brussXYZT_blockdiagonal.C

brussXYZT_blockdiagonal_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-epetra/libnoxepetra.a \
        $(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-epetra/liblocaepetra.a

brussXYZT_blockdiagonal_LDADD = $(TEST_LIBS)


brussXYZT_bdsdt_SOURCES = \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Brusselator.C $(srcdir)/Brusselator.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H \
	$(srcdir)/brussXYZT_bdsdt.C

brussXYZT_bdsdt_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-epetra/libnoxepetra.a \
        $(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-epetra/liblocaepetra.a

brussXYZT_bdsdt_LDADD = $(TEST_LIBS)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(CXXREPOS)

