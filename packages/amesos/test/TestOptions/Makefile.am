#@HEADER
#@HEADER

EXEEXT = .exe

TEST = TestOptions 

noinst_PROGRAMS = $(TEST) Suppressions

if HAVE_AMESOS_EPETRAEXT
   EPETRA_EXT_LIB = 	$(top_builddir)/../epetraext/src/libepetraext.a 
   EPETRA_EXT_INCLUDE =	-I$(top_srcdir)/../epetraext/src/transform -I$(top_builddir)/../epetraext/src/transform -I$(top_srcdir)/../epetraext/src -I$(top_builddir)/../epetraext/src 
else
   EPETRA_EXT_LIB = 
   EPETRA_EXT_INCLUDE = 
endif 

#Specify locations of include files within Trilinos, but outside amesos
AM_CPPFLAGS = -I$(top_srcdir)/../triutils/src \
-I$(top_srcdir)/../epetra/src  -I$(top_builddir)/../epetra/src \
	${EPETRA_EXT_INCLUDE} \
-I$(top_builddir)/../triutils/src \
   -I$(top_builddir)/../teuchos/src -I$(top_srcdir)/../teuchos/src

Suppressions_SOURCES =  $(srcdir)/Suppressions

TestOptions_SOURCES = $(srcdir)/TestOptions.cpp  \
	$(srcdir)/CreateTridi.cpp \
	$(srcdir)/CreateTridi.h \
	$(srcdir)/TestSuperludist.cpp \
	$(srcdir)/TestScalapack.cpp \
	$(srcdir)/TestScalapack.h \
	$(srcdir)/TestKlu.cpp \
	$(srcdir)/TestKlu.h \
	$(srcdir)/TestSuperludist.h \
	$(srcdir)/TestAllClasses.cpp \
	$(srcdir)/TestAllClasses.h \
	$(srcdir)/PerformOneSolveAndTest.cpp \
	$(srcdir)/PerformOneSolveAndTest.h \
	$(srcdir)/TestOtherClasses.cpp \
	$(srcdir)/TestOtherClasses.h

TestOptions_DEPENDENCIES = \
	$(top_builddir)/../triutils/src/libtriutils.a \
	${EPETRA_EXT_LIB} \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libamesos.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../triutils/src/libtriutils.a 

TestOptions_CXXFLAGS = -I$(top_srcdir)/src

#
#  libamesos.a is loaded twice.  The second time allows umfpack to pick up AMD.
#
TestOptions_LDADD = $(top_builddir)/src/libamesos.a \
	${EPETRA_EXT_LIB} \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../triutils/src/libtriutils.a \
	$(AMESOS_LIBS) \
        $(top_builddir)/src/libamesos.a \
        $(LAPACK_LIBS) \
	$(BLAS_LIBS) $(FLIBS)  $(FLIBS)

Suppressions$(EXEEXT):  $(srcdir)/Suppressions
	cp $(srcdir)/Suppressions Suppressions$(EXEEXT)
