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
TEST_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(NOX_TEST_LIBS))
TEST_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(NOX_TEST_INCLUDES))
else
TEST_LIBS     = $(NOX_TEST_LIBS)
TEST_INCLUDES = $(NOX_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(TEST_INCLUDES)

EXEEXT = .exe

if BUILD_BELOS_TESTS
1DFEM_BELOS = 1DfemBelos
endif

if BUILD_EPETRA_TESTS

if BUILD_NOX_EPETRAEXT_SUPPORT
1DFEM_WITHMODELEVALUATOR = 1DfemWithModelEvaluator
endif

noinst_PROGRAMS = 1Dfem 1DfemScaled 1DfemScaledInexactTrustRegion 1DfemTransient 1DfemMLPreconditioner 1DfemAztecSolverPrecReuse 1DfemAztecSolverPrecRebuild 1DfemAztecSolverPrecRecompute 1DfemAztecPrecRebuild 1DfemAztecPrecReuse 1DfemIfpackReuse 1DfemIfpackRecompute 1DfemIfpackRebuild 1DfemNonlinearCG $(1DFEM_BELOS) $(1DFEM_WITHMODELEVALUATOR) 

endif

1Dfem_SOURCES = \
	        $(srcdir)/1DfemInterface.H \
                $(srcdir)/1DfemInterface.C \
                $(srcdir)/1DfemPrePostOperator.H \
                $(srcdir)/1DfemPrePostOperator.C \
                $(srcdir)/1Dfem.C

1Dfem_DEPENDENCIES = \
                     $(top_builddir)/src/libnox.a \
                     $(top_builddir)/src-epetra/libnoxepetra.a

1Dfem_LDADD = $(TEST_LIBS)

1DfemScaled_SOURCES = \
	        $(srcdir)/1DfemInterface.H \
                $(srcdir)/1DfemInterface.C \
                $(srcdir)/1DfemPrePostOperator.H \
                $(srcdir)/1DfemPrePostOperator.C \
                $(srcdir)/1DfemScaled.C

1DfemScaled_DEPENDENCIES = \
                     $(top_builddir)/src/libnox.a \
                     $(top_builddir)/src-epetra/libnoxepetra.a

1DfemScaled_LDADD = $(TEST_LIBS)

1DfemScaledInexactTrustRegion_SOURCES = \
	        $(srcdir)/1DfemInterface.H \
                $(srcdir)/1DfemInterface.C \
                $(srcdir)/1DfemPrePostOperator.H \
                $(srcdir)/1DfemPrePostOperator.C \
                $(srcdir)/1DfemScaledInexactTrustRegion.C

1DfemScaledInexactTrustRegion_DEPENDENCIES = \
                     $(top_builddir)/src/libnox.a \
                     $(top_builddir)/src-epetra/libnoxepetra.a

1DfemScaledInexactTrustRegion_LDADD = $(TEST_LIBS)

1DfemTransient_SOURCES = \
	                 $(srcdir)/1DfemInterface.H \
                         $(srcdir)/1DfemInterface.C \
                         $(srcdir)/1DfemTransient.C

1DfemTransient_DEPENDENCIES = \
                              $(top_builddir)/src/libnox.a \
                              $(top_builddir)/src-epetra/libnoxepetra.a

1DfemTransient_LDADD = $(TEST_LIBS)

1DfemMLPreconditioner_SOURCES = \
	                 $(srcdir)/1DfemInterface.H \
                         $(srcdir)/1DfemInterface.C \
                         $(srcdir)/1DfemMLPreconditioner.C

1DfemMLPreconditioner_DEPENDENCIES = \
                              $(top_builddir)/src/libnox.a \
                              $(top_builddir)/src-epetra/libnoxepetra.a

1DfemMLPreconditioner_LDADD = $(TEST_LIBS)

1DfemBelos_SOURCES = \
	                 $(srcdir)/1DfemInterface.H \
                         $(srcdir)/1DfemInterface.C \
                         $(srcdir)/1DfemBelos.C

1DfemBelos_DEPENDENCIES = \
                              $(top_builddir)/src/libnox.a \
                              $(top_builddir)/src-epetra/libnoxepetra.a \
                              $(top_builddir)/src-belos/libnoxbelos.a

1DfemBelos_LDADD = $(TEST_LIBS)

1DfemWithModelEvaluator_SOURCES = \
                         $(srcdir)/1DfemPrePostOperator.H \
                         $(srcdir)/1DfemPrePostOperator.C \
	                 $(srcdir)/1DfemModelEvaluator.H \
                         $(srcdir)/1DfemModelEvaluator.C \
                         $(srcdir)/1DfemWithModelEvaluator.C

1DfemWithModelEvaluator_DEPENDENCIES = \
                              $(top_builddir)/src/libnox.a \
                              $(top_builddir)/src-epetra/libnoxepetra.a

1DfemWithModelEvaluator_LDADD = $(TEST_LIBS)

1DfemAztecSolverPrecReuse_SOURCES = \
                         $(srcdir)/1DfemInterface.H \
                         $(srcdir)/1DfemInterface.C \
                         $(srcdir)/1DfemAztecSolverPrecReuse.C

1DfemAztecSolverPrecReuse_DEPENDENCIES = \
                              $(top_builddir)/src/libnox.a \
                              $(top_builddir)/src-epetra/libnoxepetra.a

1DfemAztecSolverPrecReuse_LDADD = $(TEST_LIBS)

1DfemAztecSolverPrecRebuild_SOURCES = \
                         $(srcdir)/1DfemInterface.H \
                         $(srcdir)/1DfemInterface.C \
                         $(srcdir)/1DfemAztecSolverPrecRebuild.C

1DfemAztecSolverPrecRebuild_DEPENDENCIES = \
                              $(top_builddir)/src/libnox.a \
                              $(top_builddir)/src-epetra/libnoxepetra.a

1DfemAztecSolverPrecRebuild_LDADD = $(TEST_LIBS)

1DfemAztecSolverPrecRecompute_SOURCES = \
                         $(srcdir)/1DfemInterface.H \
                         $(srcdir)/1DfemInterface.C \
                         $(srcdir)/1DfemAztecSolverPrecRecompute.C

1DfemAztecSolverPrecRecompute_DEPENDENCIES = \
                              $(top_builddir)/src/libnox.a \
                              $(top_builddir)/src-epetra/libnoxepetra.a

1DfemAztecSolverPrecRecompute_LDADD = $(TEST_LIBS)

1DfemAztecPrecRebuild_SOURCES = \
                         $(srcdir)/1DfemInterface.H \
                         $(srcdir)/1DfemInterface.C \
                         $(srcdir)/1DfemAztecPrecRebuild.C

1DfemAztecPrecRebuild_DEPENDENCIES = \
                              $(top_builddir)/src/libnox.a \
                              $(top_builddir)/src-epetra/libnoxepetra.a

1DfemAztecPrecRebuild_LDADD = $(TEST_LIBS)

1DfemAztecPrecReuse_SOURCES = \
                         $(srcdir)/1DfemInterface.H \
                         $(srcdir)/1DfemInterface.C \
                         $(srcdir)/1DfemAztecPrecReuse.C

1DfemAztecPrecReuse_DEPENDENCIES = \
                              $(top_builddir)/src/libnox.a \
                              $(top_builddir)/src-epetra/libnoxepetra.a

1DfemAztecPrecReuse_LDADD = $(TEST_LIBS)

1DfemIfpackReuse_SOURCES = \
                         $(srcdir)/1DfemInterface.H \
                         $(srcdir)/1DfemInterface.C \
                         $(srcdir)/1DfemIfpackReuse.C

1DfemIfpackReuse_DEPENDENCIES = \
                              $(top_builddir)/src/libnox.a \
                              $(top_builddir)/src-epetra/libnoxepetra.a

1DfemIfpackReuse_LDADD = $(TEST_LIBS)

1DfemIfpackRecompute_SOURCES = \
                         $(srcdir)/1DfemInterface.H \
                         $(srcdir)/1DfemInterface.C \
                         $(srcdir)/1DfemIfpackRecompute.C

1DfemIfpackRecompute_DEPENDENCIES = \
                              $(top_builddir)/src/libnox.a \
                              $(top_builddir)/src-epetra/libnoxepetra.a

1DfemIfpackRecompute_LDADD = $(TEST_LIBS)

1DfemIfpackRebuild_SOURCES = \
                         $(srcdir)/1DfemInterface.H \
                         $(srcdir)/1DfemInterface.C \
                         $(srcdir)/1DfemIfpackRebuild.C

1DfemIfpackRebuild_DEPENDENCIES = \
                              $(top_builddir)/src/libnox.a \
                              $(top_builddir)/src-epetra/libnoxepetra.a

1DfemIfpackRebuild_LDADD = $(TEST_LIBS)

1DfemNonlinearCG_SOURCES = \
                         $(srcdir)/1DfemInterface.H \
                         $(srcdir)/1DfemInterface.C \
                         $(srcdir)/1DfemNonlinearCG.C

1DfemNonlinearCG_DEPENDENCIES = \
                              $(top_builddir)/src/libnox.a \
                              $(top_builddir)/src-epetra/libnoxepetra.a

1DfemNonlinearCG_LDADD = $(TEST_LIBS)
