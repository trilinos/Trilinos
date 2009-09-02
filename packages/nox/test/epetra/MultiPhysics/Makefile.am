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

if BUILD_EPETRA_TESTS
if BUILD_EPETRA_MULTIPHYSICS

## Disabling mpi
if HAVE_MPI
##
else
noinst_PROGRAMS = example_hmx example_brusselator example_yeckel Dummy
endif

endif
endif

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  -r $(srcdir)/gold_brusselator gold_brusselator
	cp  -r $(srcdir)/gold_hmx gold_hmx

Dummy_SOURCES = Dummy

example_hmx_SOURCES = $(srcdir)/Example_HMX.C \
	$(srcdir)/HMX_PDE.C $(srcdir)/HMX_PDE.H \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/GenericEpetraProblem.C $(srcdir)/GenericEpetraProblem.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H \
	$(srcdir)/NOX_Matlab_Interface.C $(srcdir)/NOX_Matlab_Interface.H \
	$(srcdir)/NOX_Multiphysics_Matlab_Interface.C $(srcdir)/NOX_Multiphysics_Matlab_Interface.H \
	$(srcdir)/Problem_Manager.C $(srcdir)/Problem_Manager.H \
	$(srcdir)/OffBlock_Manager.C $(srcdir)/OffBlock_Manager.H \
	$(srcdir)/Xfer_Operator.C $(srcdir)/Xfer_Operator.H \
	$(srcdir)/NOX_Epetra_SchurOperator.C $(srcdir)/NOX_Epetra_SchurOperator.H \
	$(srcdir)/NOX_Epetra_SchurInterface.H \
	$(srcdir)/NOX_Epetra_SchurCoupler.C $(srcdir)/NOX_Epetra_SchurCoupler.H

example_hmx_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-epetra/libnoxepetra.a

example_hmx_LDADD = $(TEST_LIBS)

example_brusselator_SOURCES = $(srcdir)/Example_Brusselator.C \
	$(srcdir)/Equation_A.C $(srcdir)/Equation_A.H \
	$(srcdir)/Equation_B.C $(srcdir)/Equation_B.H \
	$(srcdir)/Burgers.C $(srcdir)/Burgers.H \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/GenericEpetraProblem.C $(srcdir)/GenericEpetraProblem.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H \
	$(srcdir)/NOX_Matlab_Interface.C $(srcdir)/NOX_Matlab_Interface.H \
	$(srcdir)/NOX_Multiphysics_Matlab_Interface.C $(srcdir)/NOX_Multiphysics_Matlab_Interface.H \
	$(srcdir)/Problem_Manager.C $(srcdir)/Problem_Manager.H \
	$(srcdir)/OffBlock_Manager.C $(srcdir)/OffBlock_Manager.H \
	$(srcdir)/Xfer_Operator.C $(srcdir)/Xfer_Operator.H \
	$(srcdir)/NOX_Epetra_SchurOperator.C $(srcdir)/NOX_Epetra_SchurOperator.H \
	$(srcdir)/NOX_Epetra_SchurInterface.H \
	$(srcdir)/NOX_Epetra_SchurCoupler.C $(srcdir)/NOX_Epetra_SchurCoupler.H

example_brusselator_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-epetra/libnoxepetra.a

example_brusselator_LDADD = $(TEST_LIBS)

example_yeckel_SOURCES = $(srcdir)/Example_Yeckel2006.C \
	$(srcdir)/ConvDiff_PDE.C $(srcdir)/ConvDiff_PDE.H \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/GenericEpetraProblem.C $(srcdir)/GenericEpetraProblem.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H \
	$(srcdir)/NOX_Matlab_Interface.C $(srcdir)/NOX_Matlab_Interface.H \
	$(srcdir)/NOX_Multiphysics_Matlab_Interface.C $(srcdir)/NOX_Multiphysics_Matlab_Interface.H \
	$(srcdir)/Problem_Manager.C $(srcdir)/Problem_Manager.H \
	$(srcdir)/OffBlock_Manager.C $(srcdir)/OffBlock_Manager.H \
	$(srcdir)/Xfer_Operator.C $(srcdir)/Xfer_Operator.H \
	$(srcdir)/NOX_Epetra_SchurOperator.C $(srcdir)/NOX_Epetra_SchurOperator.H \
	$(srcdir)/NOX_Epetra_SchurInterface.H \
	$(srcdir)/NOX_Epetra_SchurCoupler.C $(srcdir)/NOX_Epetra_SchurCoupler.H

example_yeckel_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-epetra/libnoxepetra.a

example_yeckel_LDADD = $(TEST_LIBS)

EXTRA_DIST = \
	gold_brusselator/output.001.000_00000_map \
	gold_brusselator/output.001.000_00000_mesh_map \
	gold_brusselator/output.001.000_00000_mesh_vec \
	gold_brusselator/output.001.000_00000_vec \
	gold_brusselator/output.001.000_00001_map \
	gold_brusselator/output.001.000_00001_mesh_map \
	gold_brusselator/output.001.000_00001_mesh_vec \
	gold_brusselator/output.001.000_00001_vec \
	gold_brusselator/output.001.000_00002_map \
	gold_brusselator/output.001.000_00002_mesh_map \
	gold_brusselator/output.001.000_00002_mesh_vec \
	gold_brusselator/output.001.000_00002_vec \
	gold_brusselator/output.001.000_00003_map \
	gold_brusselator/output.001.000_00003_mesh_map \
	gold_brusselator/output.001.000_00003_mesh_vec \
	gold_brusselator/output.001.000_00003_vec \
	gold_brusselator/output.002.000_00000_map \
	gold_brusselator/output.002.000_00000_mesh_map \
	gold_brusselator/output.002.000_00000_mesh_vec \
	gold_brusselator/output.002.000_00000_vec \
	gold_brusselator/output.002.000_00001_map \
	gold_brusselator/output.002.000_00001_mesh_map \
	gold_brusselator/output.002.000_00001_mesh_vec \
	gold_brusselator/output.002.000_00001_vec \
	gold_brusselator/output.002.000_00002_map \
	gold_brusselator/output.002.000_00002_mesh_map \
	gold_brusselator/output.002.000_00002_mesh_vec \
	gold_brusselator/output.002.000_00002_vec \
	gold_brusselator/output.002.000_00003_map \
	gold_brusselator/output.002.000_00003_mesh_map \
	gold_brusselator/output.002.000_00003_mesh_vec \
	gold_brusselator/output.002.000_00003_vec \
	gold_brusselator/output.003.000_00000_map \
	gold_brusselator/output.003.000_00000_mesh_map \
	gold_brusselator/output.003.000_00000_mesh_vec \
	gold_brusselator/output.003.000_00000_vec \
	gold_brusselator/output.003.000_00001_map \
	gold_brusselator/output.003.000_00001_mesh_map \
	gold_brusselator/output.003.000_00001_mesh_vec \
	gold_brusselator/output.003.000_00001_vec \
	gold_brusselator/output.003.000_00002_map \
	gold_brusselator/output.003.000_00002_mesh_map \
	gold_brusselator/output.003.000_00002_mesh_vec \
	gold_brusselator/output.003.000_00002_vec \
	gold_brusselator/output.003.000_00003_map \
	gold_brusselator/output.003.000_00003_mesh_map \
	gold_brusselator/output.003.000_00003_mesh_vec \
	gold_brusselator/output.003.000_00003_vec \
	gold_hmx/output.001.000_00000_map \
	gold_hmx/output.001.000_00000_mesh_map \
	gold_hmx/output.001.000_00000_mesh_vec \
	gold_hmx/output.001.000_00000_vec \
	gold_hmx/output.001.000_00001_map \
	gold_hmx/output.001.000_00001_mesh_map \
	gold_hmx/output.001.000_00001_mesh_vec \
	gold_hmx/output.001.000_00001_vec \
	gold_hmx/output.001.000_00002_map \
	gold_hmx/output.001.000_00002_mesh_map \
	gold_hmx/output.001.000_00002_mesh_vec \
	gold_hmx/output.001.000_00002_vec \
	gold_hmx/output.001.000_00003_map \
	gold_hmx/output.001.000_00003_mesh_map \
	gold_hmx/output.001.000_00003_mesh_vec \
	gold_hmx/output.001.000_00003_vec \
	gold_hmx/output.002.000_00000_map \
	gold_hmx/output.002.000_00000_mesh_map \
	gold_hmx/output.002.000_00000_mesh_vec \
	gold_hmx/output.002.000_00000_vec \
	gold_hmx/output.002.000_00001_map \
	gold_hmx/output.002.000_00001_mesh_map \
	gold_hmx/output.002.000_00001_mesh_vec \
	gold_hmx/output.002.000_00001_vec \
	gold_hmx/output.002.000_00002_map \
	gold_hmx/output.002.000_00002_mesh_map \
	gold_hmx/output.002.000_00002_mesh_vec \
	gold_hmx/output.002.000_00002_vec \
	gold_hmx/output.002.000_00003_map \
	gold_hmx/output.002.000_00003_mesh_map \
	gold_hmx/output.002.000_00003_mesh_vec \
	gold_hmx/output.002.000_00003_vec \
	gold_hmx/output.003.000_00000_map \
	gold_hmx/output.003.000_00000_mesh_map \
	gold_hmx/output.003.000_00000_mesh_vec \
	gold_hmx/output.003.000_00000_vec \
	gold_hmx/output.003.000_00001_map \
	gold_hmx/output.003.000_00001_mesh_map \
	gold_hmx/output.003.000_00001_mesh_vec \
	gold_hmx/output.003.000_00001_vec \
	gold_hmx/output.003.000_00002_map \
	gold_hmx/output.003.000_00002_mesh_map \
	gold_hmx/output.003.000_00002_mesh_vec \
	gold_hmx/output.003.000_00002_vec \
	gold_hmx/output.003.000_00003_map \
	gold_hmx/output.003.000_00003_mesh_map \
	gold_hmx/output.003.000_00003_mesh_vec \
	gold_hmx/output.003.000_00003_vec \
	gold_hmx/output.004.000_00000_map \
	gold_hmx/output.004.000_00000_mesh_map \
	gold_hmx/output.004.000_00000_mesh_vec \
	gold_hmx/output.004.000_00000_vec \
	gold_hmx/output.004.000_00001_map \
	gold_hmx/output.004.000_00001_mesh_map \
	gold_hmx/output.004.000_00001_mesh_vec \
	gold_hmx/output.004.000_00001_vec \
	gold_hmx/output.004.000_00002_map \
	gold_hmx/output.004.000_00002_mesh_map \
	gold_hmx/output.004.000_00002_mesh_vec \
	gold_hmx/output.004.000_00002_vec \
	gold_hmx/output.004.000_00003_map \
	gold_hmx/output.004.000_00003_mesh_map \
	gold_hmx/output.004.000_00003_mesh_vec \
	gold_hmx/output.004.000_00003_vec \
	gold_hmx/output.004.000_00000_map \
	gold_hmx/output.004.000_00000_mesh_map \
	gold_hmx/output.004.000_00000_mesh_vec \
	gold_hmx/output.004.000_00000_vec \
	gold_hmx/output.004.000_00001_map \
	gold_hmx/output.004.000_00001_mesh_map \
	gold_hmx/output.004.000_00001_mesh_vec \
	gold_hmx/output.004.000_00001_vec \
	gold_hmx/output.004.000_00002_map \
	gold_hmx/output.004.000_00002_mesh_map \
	gold_hmx/output.004.000_00002_mesh_vec \
	gold_hmx/output.004.000_00002_vec \
	gold_hmx/output.004.000_00003_map \
	gold_hmx/output.004.000_00003_mesh_map \
	gold_hmx/output.004.000_00003_mesh_vec \
	gold_hmx/output.004.000_00003_vec
