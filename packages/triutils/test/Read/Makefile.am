#@HEADER
#@HEADER

EXEEXT = .exe

#
#

noinst_PROGRAMS = read_test Dummy


Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(srcdir)/In_bcsstk01.rsa  bcsstk01.rsa
	cp  $(srcdir)/In_bcsstk01.mtx  bcsstk01.mtx
	cp  $(srcdir)/In_bcsstk01.triS  bcsstk01.triS

	cp  $(srcdir)/In_fs_183_4.rua  fs_183_4.rua
	cp  $(srcdir)/In_fs_183_4.mtx  fs_183_4.mtx
	cp  $(srcdir)/In_fs_183_4.triU  fs_183_4.triU

	cp  $(srcdir)/In_impcol_a.rua  impcol_a.rua
	cp  $(srcdir)/In_impcol_a.mtx  impcol_a.mtx
	cp  $(srcdir)/In_impcol_a.triU  impcol_a.triU




SRC_CXX =  $(srcdir)/TestTriutilReads.cpp

read_test_SOURCES = 	 $(SRC_CXX)

#Specify locations of include files within Trilinos, but outside amesos
AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src  -I$(top_srcdir)/../amd/src \
	-I$(top_builddir)/../epetra/src \
	-I$(top_srcdir)/../triutils/src -I$(top_builddir)/../triutils/src 

read_test_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libtriutils.a 

read_test_CXXFLAGS = -I$(top_srcdir)/src

read_test_LDADD = $(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libtriutils.a \
	$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) 

# SUBDIRS = example1
SUBDIRS = 

