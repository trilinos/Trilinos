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

EXEEXT = .exe
# Nothing should be done for a serial build
if  HAVE_MPI

if HAVE_AMESOS_SLUD
SLUDIST_TEST = TestOldSuperludist 
else
SLUDIST_TEST = 
endif

if HAVE_AMESOS_DSCPACK
DSCPACK_TEST = TestDscpack
else
DSCPACK_TEST = 
endif

if HAVE_AMESOS_UMFPACK
UMFPACK_TEST = TestUmfpack
else
UMFPACK_TEST = 
endif

if HAVE_AMESOS_SCALAPACK
SCALAPACK_TEST = TestScalapack
else
SCALAPACK_TEST = 
endif

if HAVE_AMESOS_SUPERLU
SUPERLU_TEST = Test_Superlu
else
SUPERLU_TEST = 
endif

if HAVE_AMESOS_KLU
KLU_TEST = TestKlu
else
KLU_TEST = 
endif

if HAVE_AMESOS_SUPERLUDIST
SUPERLUDIST_TEST = TestSuperlud
else
SUPERLUDIST_TEST = 
endif

if HAVE_AMESOS_MUMPS
MUMPS_TEST = TestMumps
else
MUMPS_TEST = 
endif


noinst_PROGRAMS = $(SLUDIST_TEST) $(DSCPACK_TEST) $(UMFPACK_TEST) \
	$(MUMPS_TEST) $(KLU_TEST) $(SUPERLUDIST_TEST) $(SCALAPACK_TEST)

if HAVE_AMESOS_SLUD
TestOldSuperludist_SOURCES = $(srcdir)/TestOldSuperludist

TestOldSuperludist$(EXEEXT): $(srcdir)/TestOldSuperludist
	cp $(srcdir)/TestOldSuperludist TestOldSuperludist$(EXEEXT) 
endif

if HAVE_AMESOS_DSCPACK
TestDscpack_SOURCES = $(srcdir)/TestDscpack

TestDscpack$(EXEEXT): $(srcdir)/TestDscpack
	cp $(srcdir)/TestDscpack TestDscpack$(EXEEXT) 
endif

if HAVE_AMESOS_UMFPACK
TestUmfpack_SOURCES = $(srcdir)/TestUmfpack

TestUmfpack$(EXEEXT): $(srcdir)/TestUmfpack
	cp $(srcdir)/TestUmfpack TestUmfpack$(EXEEXT) 
endif

if HAVE_AMESOS_SCALAPACK
TestScalapack_SOURCES = $(srcdir)/TestScalapack

TestScalapack$(EXEEXT): $(srcdir)/TestScalapack
	cp $(srcdir)/TestScalapack TestScalapack$(EXEEXT) 
endif

if HAVE_AMESOS_KLU
TestKlu_SOURCES = $(srcdir)/TestKlu

TestKlu$(EXEEXT): $(srcdir)/TestKlu
	cp $(srcdir)/TestKlu TestKlu$(EXEEXT) 
endif

if HAVE_AMESOS_SUPERLUDIST
TestSuperlud_SOURCES = $(srcdir)/TestSuperlud

TestSuperlud$(EXEEXT): $(srcdir)/TestSuperlud
	cp $(srcdir)/TestSuperlud TestSuperlud$(EXEEXT) 
endif

if HAVE_AMESOS_MUMPS
TestMumps_SOURCES = $(srcdir)/TestMumps

TestMumps$(EXEEXT): $(srcdir)/TestMumps
	cp $(srcdir)/TestMumps TestMumps$(EXEEXT) 
endif

endif

