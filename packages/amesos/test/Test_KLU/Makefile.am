#@HEADER
#@HEADER

EXEEXT = .exe
TEST = Test_KLU

noinst_PROGRAMS = $(TEST)

#Specify locations of include files within Trilinos, but outside amesos
AM_CPPFLAGS = \
	-I$(top_builddir)/../epetra/src  -I$(top_srcdir)/../epetra/src \
	-I$(top_builddir)/../teuchos/src -I$(top_srcdir)/../teuchos/src

Test_KLU_SOURCES = $(srcdir)/cxx_main.cpp

Test_KLU_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libamesos.a \
	$(top_builddir)/../teuchos/src/libteuchos.a

Test_KLU_CXXFLAGS = -I$(top_srcdir)/src

#
#  libamesos.a is loaded twice.  The second time allows umfpack to pick up AMD.
#
Test_KLU_LDADD = $(top_builddir)/src/libamesos.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(AMESOS_LIBS) \
	$(top_builddir)/src/libamesos.a \
        $(LAPACK_LIBS) \
	$(BLAS_LIBS) $(FLIBS)  $(FLIBS)
