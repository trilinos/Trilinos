#@HEADER
#@HEADER

EXAMPLE = UG_ExMaps

noinst_PROGRAMS = $(EXAMPLE)

UG_ExMaps_SOURCES = $(srcdir)/cxx_main.cpp

UG_ExMaps_DEPENDENCIES = $(top_builddir)/src/libepetra.a

UG_ExMaps_CXXFLAGS = -I$(top_srcdir)/src

UG_ExMaps_LDADD = $(top_builddir)/src/libepetra.a $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)
