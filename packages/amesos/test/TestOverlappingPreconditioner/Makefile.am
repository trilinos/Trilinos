#@HEADER
#@HEADER

TEST = TestPreconditioner

EXEEXT = .exe

noinst_PROGRAMS = $(TEST)

#Specify locations of include files within Trilinos, but outside amesos
AM_CPPFLAGS = \
	-I$(top_srcdir)/../triutils/src   \
	-I$(top_srcdir)/../epetra/src     \
	-I$(top_builddir)/../epetra/src   \
	-I$(top_builddir)/../triutils/src \
	-I$(top_builddir)/../aztecoo/src  \
	-I$(top_srcdir)/../aztecoo/src    \
	-I$(top_builddir)/../teuchos/src  \
	-I$(top_srcdir)/../teuchos/src

TestPreconditioner_SOURCES = $(srcdir)/TestPreconditioner.cpp

TestPreconditioner_DEPENDENCIES = \
	$(top_builddir)/../triutils/src/libtriutils.a \
	$(top_builddir)/../epetra/src/libepetra.a     \
	$(top_builddir)/src/libamesos.a               \
	$(top_builddir)/../teuchos/src/libteuchos.a   \
	$(top_builddir)/../aztecoo/src/libaztecoo.a   \
	$(top_builddir)/../triutils/src/libtriutils.a 

TestPreconditioner_CXXFLAGS = -I$(top_srcdir)/src

#
#  libamesos.a is loaded twice.  The second time allows umfpack to pick up AMD.
#
TestPreconditioner_LDADD = \
	$(top_builddir)/src/libamesos.a               \
	$(top_builddir)/../aztecoo/src/libaztecoo.a   \
	$(top_builddir)/../teuchos/src/libteuchos.a   \
	$(top_builddir)/../triutils/src/libtriutils.a \
	$(top_builddir)/../epetra/src/libepetra.a     \
	$(AMESOS_LIBS)                                \
	$(top_builddir)/src/libamesos.a               \
        $(LAPACK_LIBS) \
	$(BLAS_LIBS) $(FLIBS)  $(FLIBS)

