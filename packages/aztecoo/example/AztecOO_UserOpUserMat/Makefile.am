#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = Poisson2d

Poisson2d_SOURCES = $(srcdir)/cxx_main.cpp $(srcdir)/Poisson2dOperator.cpp \
		$(srcdir)/Poisson2dOperator.h \
		$(top_srcdir)/../triutils/src/Trilinos_Util.h
#Poisson2d_DEPENDENCIES = $(top_builddir)/src/libaztecoo.a 

Poisson2d_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/../epetra/src \
	-I$(top_srcdir)/../triutils/src -I$(top_builddir)/../epetra/src \
	-I$(top_builddir)/../triutils/src

Poisson2d_LDADD = $(top_builddir)/src/libaztecoo.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)


