#@HEADER
#@HEADER

EXAMPLE = UG_Ex1

noinst_PROGRAMS = $(EXAMPLE)

UG_Ex1_SOURCES = $(srcdir)/cxx_main.cpp

UG_Ex1_DEPENDENCIES = $(top_builddir)/src/libepetra.a

UG_Ex1_CXXFLAGS = -I$(top_srcdir)/src

UG_Ex1_LDADD = $(top_builddir)/src/libepetra.a $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)
