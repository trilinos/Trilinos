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

#
#  All of the AmeosShellScript and TestAmesosScript stuff is designed to 
#  copy some files form the source tree to the build tree
#
#  The idea behind the mechanism is to fool autotools into thinking that
#  AmesosShellScript is an executable and that the way that you create
#  that executable is the "cp" command.
#
#  The files are renamed as they are copied so that configure can work 
#  within the source directory (even though we recommend that users 
#  run configure from a separate build directory.
#

noinst_PROGRAMS = amesos_test AmesosShellScript AmesosDscpack AmesosDscpackShort \
	AmesosShortSuperLU AmesosShortKundert AmesosDistScript AmesosShortScript \
	AmesosJanusScript SpoolesJanusScript SpoolesScript AmesosSuperLUserial \
	AmesosUmfpack AmesosUmfpackShort AmesosUmfpackSerial AmesosMumps \
	AmesosKlu AmesosKluShort AmesosKluSerial \
	AmesosScalapack \
	Amesos_Superlu Amesos_SuperluSerial \
	AmesosSuperludist AmesosSuperludistShort \
	AmesosSuperludistJanus AmesosUmfpackJanus AmesosKluJanus \
	A_Superludist  B_Superludist  C_Superludist Timings DscpackTiming 

AmesosShellScript_SOURCES = $(srcdir)/AmesosShellScript
A_Superludist_SOURCES = $(srcdir)/A_Superludist
B_Superludist_SOURCES = $(srcdir)/B_Superludist
C_Superludist_SOURCES = $(srcdir)/C_Superludist
Timings_SOURCES = $(srcdir)/Timings
AmesosDscpack_SOURCES = $(srcdir)/AmesosDscpack
AmesosDscpackShort_SOURCES = $(srcdir)/AmesosDscpackShort
AmesosMumps_SOURCES = $(srcdir)/AmesosMumps
AmesosUmfpack_SOURCES = $(srcdir)/AmesosUmfpack
AmesosUmfpackShort_SOURCES = $(srcdir)/AmesosUmfpackShort
AmesosUmfpackSerial_SOURCES = $(srcdir)/AmesosUmfpackSerial
AmesosScalapack_SOURCES = $(srcdir)/AmesosScalapack
AmesosKlu_SOURCES = $(srcdir)/AmesosKlu
AmesosKluShort_SOURCES = $(srcdir)/AmesosKluShort
AmesosKluSerial_SOURCES = $(srcdir)/AmesosKluSerial
Amesos_Superlu_SOURCES = $(srcdir)/Amesos_Superlu
Amesos_SuperluSerial_SOURCES = $(srcdir)/Amesos_SuperluSerial
AmesosSuperludist_SOURCES = $(srcdir)/AmesosSuperludist
AmesosSuperludistShort_SOURCES = $(srcdir)/AmesosSuperludistShort
AmesosShortSuperLU_SOURCES = $(srcdir)/AmesosShortSuperLU
AmesosShortKundert_SOURCES = $(srcdir)/AmesosShortKundert
AmesosDistScript_SOURCES = $(srcdir)/AmesosDistScript
AmesosShortScript_SOURCES = $(srcdir)/AmesosShortScript
AmesosSuperLUserial_SOURCES = $(srcdir)/AmesosSuperLUserial
AmesosJanusScript_SOURCES = $(srcdir)/AmesosJanusScript
SpoolesJanusScript_SOURCES = $(srcdir)/SpoolesJanusScript
SpoolesScript_SOURCES = $(srcdir)/SpoolesScript

AmesosShellScript$(EXEEXT): $(srcdir)/AmesosShellScript
	cp  $(srcdir)/AmesosShellScript  AmesosShellScript$(EXEEXT)

A_Superludist$(EXEEXT): $(srcdir)/A_Superludist
	cp  $(srcdir)/A_Superludist  A_Superludist$(EXEEXT)

B_Superludist$(EXEEXT): $(srcdir)/B_Superludist
	cp  $(srcdir)/B_Superludist  B_Superludist$(EXEEXT)

C_Superludist$(EXEEXT): $(srcdir)/C_Superludist
	cp  $(srcdir)/C_Superludist  C_Superludist$(EXEEXT)

Timings$(EXEEXT): $(srcdir)/Timings
	cp  $(srcdir)/Timings  Timings$(EXEEXT)

DscpackTiming$(EXEEXT): $(srcdir)/DscpackTiming
	cp  $(srcdir)/DscpackTiming  DscpackTiming$(EXEEXT)

AmesosDscpack$(EXEEXT): $(srcdir)/AmesosDscpack
	cp  $(srcdir)/AmesosDscpack  AmesosDscpack$(EXEEXT)

AmesosDscpackShort$(EXEEXT): $(srcdir)/AmesosDscpackShort
	cp  $(srcdir)/AmesosDscpackShort  AmesosDscpackShort$(EXEEXT)

AmesosMumps$(EXEEXT): $(srcdir)/AmesosMumps
	cp  $(srcdir)/AmesosMumps  AmesosMumps$(EXEEXT)

AmesosUmfpack$(EXEEXT): $(srcdir)/AmesosUmfpack
	cp  $(srcdir)/AmesosUmfpack  AmesosUmfpack$(EXEEXT)

AmesosUmfpackShort$(EXEEXT): $(srcdir)/AmesosUmfpackShort
	cp  $(srcdir)/AmesosUmfpackShort  AmesosUmfpackShort$(EXEEXT)

AmesosUmfpackSerial$(EXEEXT): $(srcdir)/AmesosUmfpackSerial
	cp  $(srcdir)/AmesosUmfpackSerial  AmesosUmfpackSerial$(EXEEXT)

AmesosScalapack$(EXEEXT): $(srcdir)/AmesosScalapack
	cp  $(srcdir)/AmesosScalapack  AmesosScalapack$(EXEEXT)

AmesosKlu$(EXEEXT): $(srcdir)/AmesosKlu
	cp  $(srcdir)/AmesosKlu  AmesosKlu$(EXEEXT)

AmesosKluShort$(EXEEXT): $(srcdir)/AmesosKluShort
	cp  $(srcdir)/AmesosKluShort  AmesosKluShort$(EXEEXT)

AmesosKluSerial$(EXEEXT): $(srcdir)/AmesosKluSerial
	cp  $(srcdir)/AmesosKluSerial  AmesosKluSerial$(EXEEXT)

Amesos_Superlu$(EXEEXT): $(srcdir)/Amesos_Superlu
	cp  $(srcdir)/Amesos_Superlu  Amesos_Superlu$(EXEEXT)

Amesos_SuperluSerial$(EXEEXT): $(srcdir)/Amesos_SuperluSerial
	cp  $(srcdir)/Amesos_SuperluSerial  Amesos_SuperluSerial$(EXEEXT)

AmesosSuperludist$(EXEEXT): $(srcdir)/AmesosSuperludist
	cp  $(srcdir)/AmesosSuperludist  AmesosSuperludist$(EXEEXT)

AmesosSuperludistShort$(EXEEXT): $(srcdir)/AmesosSuperludistShort
	cp  $(srcdir)/AmesosSuperludistShort  AmesosSuperludistShort$(EXEEXT)

AmesosShortSuperLU$(EXEEXT): $(srcdir)/AmesosShortSuperLU
	cp  $(srcdir)/AmesosShortSuperLU  AmesosShortSuperLU$(EXEEXT)

AmesosShortKundert$(EXEEXT): $(srcdir)/AmesosShortKundert
	cp  $(srcdir)/AmesosShortKundert  AmesosShortKundert$(EXEEXT)

AmesosDistScript$(EXEEXT): $(srcdir)/AmesosDistScript
	cp  $(srcdir)/AmesosDistScript  AmesosDistScript$(EXEEXT)

AmesosShortScript$(EXEEXT): $(srcdir)/AmesosShortScript
	cp  $(srcdir)/AmesosShortScript  AmesosShortScript$(EXEEXT)

AmesosSuperLUserial$(EXEEXT): $(srcdir)/AmesosSuperLUserial
	cp  $(srcdir)/AmesosSuperLUserial  AmesosSuperLUserial$(EXEEXT)

AmesosJanusScript$(EXEEXT): $(srcdir)/AmesosJanusScript
	cp  $(srcdir)/AmesosJanusScript  AmesosJanusScript$(EXEEXT)

SpoolesJanusScript$(EXEEXT): $(srcdir)/SpoolesJanusScript
	cp  $(srcdir)/SpoolesJanusScript  SpoolesJanusScript$(EXEEXT)

AmesosSuperludistJanus$(EXEEXT): $(srcdir)/AmesosSuperludistJanus
	cp  $(srcdir)/AmesosSuperludistJanus  AmesosSuperludistJanus$(EXEEXT)

AmesosKluJanus$(EXEEXT): $(srcdir)/AmesosKluJanus
	cp  $(srcdir)/AmesosKluJanus  AmesosKluJanus$(EXEEXT)

AmesosUmfpackJanus$(EXEEXT): $(srcdir)/AmesosUmfpackJanus
	cp  $(srcdir)/AmesosUmfpackJanus  AmesosUmfpackJanus$(EXEEXT)

SpoolesScript$(EXEEXT): $(srcdir)/SpoolesScript
	cp  $(srcdir)/SpoolesScript  SpoolesScript$(EXEEXT)
	cp  $(srcdir)/662_bus.rsa  662_bus_out.rsa
	cp  $(srcdir)/meg1.rua  meg1_out.rua
	cp  $(srcdir)/In_ImpcolA.rua  ImpcolA.rua
	cp  $(srcdir)/In_ImpcolB.rua  ImpcolB.rua
	cp  $(srcdir)/In_ImpcolC.rua  ImpcolC.rua
	cp  $(srcdir)/In_ImpcolD.rua  ImpcolD.rua
	cp  $(srcdir)/In_ImpcolE.rua  ImpcolE.rua
	cp  $(srcdir)/In_fidapm05.rua  fidapm05.rua
	cp  $(srcdir)/In_SuperLU.rua  SuperLU.rua
	cp  $(srcdir)/In_bcsstk18.rsa  bcsstk18.rsa
	cp  $(srcdir)/In_bcsstk24.rsa  bcsstk24.rsa
	cp  $(srcdir)/In_K4989.triS  K4989.triS
	cp  $(srcdir)/In_K5000.triS  K5000.triS
	cp  $(srcdir)/In_SuperLU.triU  SuperLU.triU
	cp  $(srcdir)/In_Khead.triS  Khead.triS
	cp  $(srcdir)/In_bcsstk01.mtx bcsstk01.mtx
	cp  $(srcdir)/In_bcsstk02.mtx bcsstk02.mtx 
	cp  $(srcdir)/In_bcsstk04.mtx bcsstk04.mtx 
	cp  $(srcdir)/In_bcsstk08.mtx bcsstk08.mtx 
	cp  $(srcdir)/In_bcsstk13.mtx bcsstk13.mtx 
	cp  $(srcdir)/In_nos1.mtx nos1.mtx 
	cp  $(srcdir)/In_nos5.mtx nos5.mtx 
	cp  $(srcdir)/In_FourByFour.mtx FourByFour.mtx 
	cp  $(srcdir)/In_Diagonal.mtx Diagonal.mtx 
	cp  $(srcdir)/In_KheadK.mtx KheadK.mtx 
	cp  $(srcdir)/In_KheadSorted.mtx KheadSorted.mtx 
	cp  $(srcdir)/In_ZD2.triU ZD2.triU 
	cp  $(srcdir)/In_ZD5.triU ZD5.triU 
	cp  $(srcdir)/In_ZD16.triU ZD16.triU 




SRC_CXX =  $(srcdir)/SparseDirectTimingVars.cpp $(srcdir)/SparseSolverResult.cpp \
	$(srcdir)/Comm_assert_equal.cpp $(srcdir)/Trilinos_Util_CountTriples.cpp \
	$(srcdir)/Trilinos_Util_ReadTriples2Epetra.cpp \
	$(srcdir)/Trilinos_Util_CountMatrixMarket.cpp \
	$(srcdir)/Trilinos_Util_ReadMatrixMarket2Epetra.cpp

amesos_test_SOURCES = $(srcdir)/Amesos_TestDriver.cpp  $(srcdir)/Amesos_TestSolver.cpp \
	 $(srcdir)/Amesos_TestMultiSolver.cpp  $(srcdir)/Amesos_TestMrhsSolver.cpp \
	 $(SRC_CXX)

#Specify locations of include files within Trilinos, but outside amesos
AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src  -I$(top_srcdir)/../amd/src \
	-I$(top_builddir)/../epetra/src \
	-I$(top_srcdir)/../teuchos/src -I$(top_builddir)/../teuchos/src \
	-I$(top_srcdir)/../triutils/src -I$(top_builddir)/../triutils/src 

amesos_test_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libamesos.a \
	$(top_builddir)/amd/src/libamd.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../triutils/src/libtriutils.a 

amesos_test_CXXFLAGS = -I$(top_srcdir)/src

amesos_test_LDADD = $(top_builddir)/src/libamesos.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../triutils/src/libtriutils.a \
	$(AMESOS_LIBS) $(LAPACK_LIBS) \
	$(BLAS_LIBS) $(FLIBS) 

#np# Set "SUBDIRS= (list all subdirectories that are autotool'ed here)"
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.
#np# For packages whose examples are located in this directory, rather than
#np# in a subdirectory of this directory, use the Makefile.am located in
#np# example1/ rather than this file.
# SUBDIRS = example1
# SUBDIRS = scripts 
SUBDIRS = scripts TestOptions TestWithTriutils

