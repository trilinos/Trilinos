#@HEADER
#@HEADER

EXEEXT = .exe

TEST = TestOptions 

noinst_PROGRAMS = $(TEST)  Dummy 

#Specify locations of include files within Trilinos, but outside amesos
AM_CPPFLAGS = -I$(top_srcdir)/../triutils/src -I$(top_srcdir)/../epetra/src \
  -I$(top_builddir)/../epetra/src -I$(top_builddir)/../triutils/src \
   -I$(top_builddir)/../teuchos/src -I$(top_srcdir)/../teuchos/src

TestOptions_SOURCES = $(srcdir)/TestOptions.cpp  \
	$(srcdir)/TestOptions.h \
	$(srcdir)/CreateTridi.cpp \
	$(srcdir)/CreateTridi.h \
	$(srcdir)/TestSuperludist.cpp \
	$(srcdir)/TestSuperludist.h \
	$(srcdir)/TestAllClasses.cpp \
	$(srcdir)/TestAllClasses.h \
	$(srcdir)/PerformOneSolveAndTest.cpp \
	$(srcdir)/PerformOneSolveAndTest.h \
	$(srcdir)/TestOtherClasses.cpp \
	$(srcdir)/TestOtherClasses.h

TestOptions_DEPENDENCIES = \
	$(top_builddir)/../triutils/src/libtriutils.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libamesos.a \
	$(top_builddir)/amd/src/libamd.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../triutils/src/libtriutils.a 

TestOptions_CXXFLAGS = -I$(top_srcdir)/src

TestOptions_LDADD = $(top_builddir)/src/libamesos.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../triutils/src/libtriutils.a \
	$(AMESOS_LIBS) \
        $(LAPACK_LIBS) \
	$(BLAS_LIBS) $(FLIBS)  $(FLIBS)

Dummy$(EXEEXT):  $(srcdir)/Dummy
	cp $(srcdir)/Dummy Dummy$(EXEEXT)
	cp $(srcdir)/Suppressions Suppressions.in 
