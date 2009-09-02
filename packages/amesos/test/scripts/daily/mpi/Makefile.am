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

#
#  This is a bit confusing because some of these tests were written by Ken Stanley,
#  while others were written by Marzio Sala.  
#
#  TestAmesos runs cxx_main.cpp tests in the following directories:
#     Test_Epetra_RowMatrix 
#     Test_Epetra_CrsMatrix 
#     Test_Epetra_VbrMatrix 
#     Test_Detailed 
#     Test_LAPACK 
#     Test_KLU 
#     Test_UMFPACK 
#     Test_SuperLU 
#     Test_SuperLU_DIST 
#     Test_MUMPS 
#     Test_DSCPACK 
#  TestBasic runs TestBasic.csh in Test_Basic
#  TestExamples runs executables in the example directories 
#  TestKlu.exe runs Test_Basic/AmesosKlu  
#  TestDscpack.exe runs Test_Basic/AmesosDscpack  
#  TestSuperludist.exe runs Test_Basic/AmesosSuperludist  
#  TestSuperlu.exe runs Test_Basic/AmesosSuperlu
#  TestUmfpack.exe runs Test_Basic/AmesosUmfpack  
#  TestMumps.exe runs Test_Basic/AmesosMumps  
#  TestTaucs.exe runs Test_Basic/AmesosTaucs  
#  TestParaklete.exe runs Test_Basic/AmesosParaklete  
#  TestPardiso.exe runs Test_Basic/AmesosPardiso  
#  TestKluSerial.exe runs Test_Basic/AmesosKluSerial
#  TestSuperluSerial.exe runs Test_Basic/AmesosSuperluSerial
#  TestUmpfpackSerial.exe runs Test_Basic/AmesosUmpfpackSerial
#
#  TestValgrind.exe runs TestOptions/TestOptions.exe (under valgrind) 
#
#

EXEEXT = .exe
# Nothing should be done for a serial build
if  HAVE_MPI


TestAmesos_SOURCES = $(srcdir)/TestAmesos.sh
TestAmesos$(EXEEXT): ${TestAmesos_SOURCES}
	cp $(srcdir)/TestAmesos.sh TestAmesos$(EXEEXT) 

TestQuietAmesos_SOURCES = $(srcdir)/TestQuietAmesos.sh
TestQuietAmesos$(EXEEXT): ${TestQuietAmesos_SOURCES}
	cp $(srcdir)/TestQuietAmesos.sh TestQuietAmesos$(EXEEXT) 

TestExamples_SOURCES = $(srcdir)/TestExamples.sh
TestExamples$(EXEEXT): ${TestExamples_SOURCES}
	cp $(srcdir)/TestExamples.sh TestExamples$(EXEEXT) 

endif

if HAVE_AMESOS_DSCPACK
DSCPACK_TEST = TestDscpack
else
DSCPACK_TEST = 
endif

if HAVE_VALGRIND
VALGRIND_TEST = TestValgrind
TestValgrind_SOURCES = $(srcdir)/TestValgrind

TestValgrind$(EXEEXT): $(srcdir)/TestValgrind
	cp $(srcdir)/TestValgrind TestValgrind$(EXEEXT) 
else
VALGRIND_TEST = 
endif

if HAVE_COVERAGE
COVERAGE_TEST = TestCoverage
TestCoverage_SOURCES = $(srcdir)/TestCoverage

TestCoverage$(EXEEXT): $(srcdir)/TestCoverage
	cp $(srcdir)/TestCoverage TestCoverage$(EXEEXT) 
else
COVERAGE_TEST = 
endif


if HAVE_AMESOS_UMFPACK
UMFPACK_TEST = TestUmfpack
else
UMFPACK_TEST = 
endif

if HAVE_AMESOS_SUPERLU
SUPERLU_TEST = TestSuperlu
else
SUPERLU_TEST = 
endif

if HAVE_AMESOS_KLU
KLU_TEST = TestKlu
TestKlu_SOURCES = $(srcdir)/TestKlu

TestKlu$(EXEEXT): $(srcdir)/TestKlu
	cp $(srcdir)/TestKlu TestKlu$(EXEEXT) 
else
KLU_TEST = 
endif

if HAVE_AMESOS_SUPERLUDIST
SUPERLUDIST_TEST = TestSuperludist
else
SUPERLUDIST_TEST = 
endif

if HAVE_AMESOS_MUMPS
MUMPS_TEST = TestMumps
else
MUMPS_TEST = 
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

if HAVE_AMESOS_LAPACK
LAPACK_TEST = TestLapack
TestLapack_SOURCES = $(srcdir)/TestLapack

TestLapack$(EXEEXT): $(srcdir)/TestLapack
	cp $(srcdir)/TestLapack TestLapack$(EXEEXT) 
else
LAPACK_TEST = 
endif


if HAVE_AMESOS_PARDISO
PARDISO_TEST = TestPardiso
TestPardiso_SOURCES = $(srcdir)/TestPardiso

TestPardiso$(EXEEXT): $(srcdir)/TestPardiso
	cp $(srcdir)/TestPardiso TestPardiso$(EXEEXT) 
else
PARDISO_TEST = 
endif


if HAVE_AMESOS_PARAKLETE
PARAKLETE_TEST = TestParaklete
TestParaklete_SOURCES = $(srcdir)/TestParaklete

TestParaklete$(EXEEXT): $(srcdir)/TestParaklete
	cp $(srcdir)/TestParaklete TestParaklete$(EXEEXT) 
else
PARAKLETE_TEST = 
endif

if HAVE_AMESOS_SCALAPACK
SCALAPACK_TEST = TestScalapack
TestScalapack_SOURCES = $(srcdir)/TestScalapack

TestScalapack$(EXEEXT): $(srcdir)/TestScalapack
	cp $(srcdir)/TestScalapack TestScalapack$(EXEEXT) 
else
SCALAPACK_TEST = 
endif

if HAVE_AMESOS_TAUCS
TAUCS_TEST = TestTaucs
TestTaucs_SOURCES = $(srcdir)/TestTaucs

TestTaucs$(EXEEXT): $(srcdir)/TestTaucs
	cp $(srcdir)/TestTaucs TestTaucs$(EXEEXT) 
else
TAUCS_TEST = 
endif

if HAVE_AMESOS_KLU
endif

if HAVE_AMESOS_SUPERLU
TestSuperlu_SOURCES = $(srcdir)/TestSuperlu

TestSuperlu$(EXEEXT): $(srcdir)/TestSuperlu
	cp $(srcdir)/TestSuperlu TestSuperlu$(EXEEXT) 
endif


if HAVE_AMESOS_SUPERLUDIST
TestSuperludist_SOURCES = $(srcdir)/TestSuperludist

TestSuperludist$(EXEEXT): $(srcdir)/TestSuperludist
	cp $(srcdir)/TestSuperludist TestSuperludist$(EXEEXT) 
endif

if HAVE_AMESOS_MUMPS
TestMumps_SOURCES = $(srcdir)/TestMumps

TestMumps$(EXEEXT): $(srcdir)/TestMumps
	cp $(srcdir)/TestMumps TestMumps$(EXEEXT) 
endif

#
#  Note:  COVERAGE_TEST must be listed last so that it will execute last
#
noinst_PROGRAMS = TestAmesos TestQuietAmesos TestExamples $(DSCPACK_TEST) \
	$(SCALAPACK_TEST) $(LAPACK_TEST) $(UMFPACK_TEST) \
	$(TAUCS_TEST) $(PARDISO_TEST) \
	$(MUMPS_TEST) $(PARAKLETE_TEST) $(KLU_TEST)  $(SUPERLU_TEST) \
	$(SUPERLUDIST_TEST) \
	$(VALGRIND_TEST)  $(COVERAGE_TEST)

