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
	-I$(top_builddir)/../triutils/src -I$(top_builddir)/../ifpack/src
	-I$(top_srcdir)/../ifpack/src

if BUILD_Y12M
REFERENCE_Y12M_LIBS=-L$(top_builddir)/../../external/y12m -ly12m
endif

if BUILD_BLAS
REFERENCE_BLAS_LIBS=-L$(top_builddir)/../../external/blas -lblas
endif
if BUILD_LAPACK
REFERENCE_LAPACK_LIBS=-L$(top_builddir)/../../external/lapack -llapack
endif

if BUILD_OLDIFPACK
OLDIFPACK_LIBS=-L$(top_builddir)/../ifpack/src -lifpack
endif 

Poisson2d_LDADD = -L$(top_builddir)/src -laztecoo $(OLDIFPACK_LIBS) \
	-L$(top_builddir)/../epetra/src -lepetra $(REFERENCE_Y12M_LIBS) \
	$(LAPACK_LIBS) $(REFERENCE_LAPACK_LIBS) $(BLAS_LIBS) \
	$(REFERENCE_BLAS_LIBS) $(LIBS) $(FLIBS)


