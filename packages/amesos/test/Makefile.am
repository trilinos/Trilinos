#@HEADER
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
	AmesosUmfpack AmesosUmfpackShort AmesosUmfpackSerial AmesosMumps

AmesosShellScript_SOURCES = $(srcdir)/AmesosShellScript
AmesosDscpack_SOURCES = $(srcdir)/AmesosDscpack
AmesosDscpackShort_SOURCES = $(srcdir)/AmesosDscpackShort
AmesosMumps_SOURCES = $(srcdir)/AmesosMumps
AmesosUmfpack_SOURCES = $(srcdir)/AmesosUmfpack
AmesosUmfpackShort_SOURCES = $(srcdir)/AmesosUmfpackShort
AmesosUmfpackSerial_SOURCES = $(srcdir)/AmesosUmfpackSerial
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

SpoolesScript$(EXEEXT): $(srcdir)/SpoolesScript
	cp  $(srcdir)/SpoolesScript  SpoolesScript$(EXEEXT)
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




SRC_CXX =  $(srcdir)/SparseDirectTimingVars.cpp $(srcdir)/SparseSolverResult.cpp \
	$(srcdir)/Comm_assert_equal.cpp $(srcdir)/Trilinos_Util_CountTriples.cpp \
	$(srcdir)/Trilinos_Util_ReadTriples2Epetra.cpp \
	$(srcdir)/Trilinos_Util_CountMatrixMarket.cpp \
	$(srcdir)/Trilinos_Util_ReadMatrixMarket2Epetra.cpp

amesos_test_SOURCES = $(srcdir)/Amesos_TestDriver.cpp  $(srcdir)/Amesos_TestSolver.cpp \
	 $(srcdir)/Amesos_TestMultiSolver.cpp  $(srcdir)/Amesos_TestMrhsSolver.cpp \
	 $(SRC_CXX)

#Specify locations of include files within Trilinos, but outside amesos
AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src \
	-I$(top_srcdir)/../triutils/src -I$(top_builddir)/../triutils/src 

amesos_test_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libamesos.a \
	$(top_builddir)/../triutils/src/libtriutils.a 

amesos_test_CXXFLAGS = -I$(top_srcdir)/src

amesos_test_LDADD = $(top_builddir)/src/libamesos.a \
	$(top_builddir)/../epetra/src/libepetra.a \
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
SUBDIRS = scripts 

