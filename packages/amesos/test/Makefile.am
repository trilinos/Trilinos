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

noinst_PROGRAMS = amesos_test AmesosShellScript

AmesosShellScript_SOURCES = $(srcdir)/AmesosShellScript

AmesosShellScript$(EXEEXT): $(srcdir)/AmesosShellScript
	cp  $(srcdir)/AmesosShellScript  AmesosShellScript$(EXEEXT)
	cp  $(srcdir)/AmesosDistScript  AmesosDistScript$(EXEEXT)
	cp  $(srcdir)/AmesosShortScript  AmesosShortScript$(EXEEXT)
	cp  $(srcdir)/AmesosJanusScript  AmesosJanusScript$(EXEEXT)
	cp  $(srcdir)/SpoolesJanusScript  SpoolesJanusScript$(EXEEXT)
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




SRC_CXX =  $(srcdir)/SparseDirectTimingVars.cpp $(srcdir)/SparseSolverResult.cpp \
	$(srcdir)/Comm_assert_equal.cpp $(srcdir)/Trilinos_Util_CountTriples.cpp \
	$(srcdir)/Trilinos_Util_ReadTriples2Epetra.cpp

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

amesos_test_LDADD = -L$(top_builddir)/src -lamesos \
	-L$(top_builddir)/../epetra/src -lepetra \
	-L$(top_builddir)/../triutils/src -ltriutils \
	$(AMESOS_LIBS) $(LAPACK_LIBS) \
	$(BLAS_LIBS) $(FLIBS)

#np# Set "SUBDIRS= (list all subdirectories that are autotool'ed here)"
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.
#np# For packages whose examples are located in this directory, rather than
#np# in a subdirectory of this directory, use the Makefile.am located in
#np# example1/ rather than this file.
# SUBDIRS = example1
SUBDIRS = 
