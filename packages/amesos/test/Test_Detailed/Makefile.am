#@HEADER
#@HEADER

EXEEXT = .exe
TEST = Test_Detailed

noinst_PROGRAMS = $(TEST)

if HAVE_AMESOS_TRIUTILS
TRIUTILS_INCDIR = -I$(top_srcdir)/../triutils/src -I$(top_builddir)/../triutils/src
TRIUTILS_LIBS = $(top_builddir)/../triutils/src/libtriutils.a
else
TRIUTILS_INCDIR = 
TRIUTILS_LIBS = 
endif

#Specify locations of include files within Trilinos, but outside amesos
AM_CPPFLAGS = $(TRIUTILS_INCDIR) \
	-I$(top_builddir)/../epetra/src  -I$(top_srcdir)/../epetra/src \
	-I$(top_builddir)/../teuchos/src -I$(top_srcdir)/../teuchos/src

Test_Detailed_SOURCES = $(srcdir)/cxx_main.cpp

Test_Detailed_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libamesos.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(TRIUTILS_LIBS)

Test_Detailed_CXXFLAGS = -I$(top_srcdir)/src

#
#  libamesos.a is loaded twice.  The second time allows umfpack to pick up AMD.
#
Test_Detailed_LDADD = $(top_builddir)/src/libamesos.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(TRIUTILS_LIBS) \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(AMESOS_LIBS) \
	$(top_builddir)/src/libamesos.a \
        $(LAPACK_LIBS) \
	$(BLAS_LIBS) $(FLIBS)  $(FLIBS)
