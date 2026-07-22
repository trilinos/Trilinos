
include $(top_builddir)/Makefile.export.zoltan

F77=$(FC)
AM_FFLAGS=$(FCFLAGS)

#if USING_GNUMAKE
#EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(ZOLTAN_LIBS))
#EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(ZOLTAN_INCLUDES))
#else
EXPORT_LIBS = $(ZOLTAN_LIBS)
EXPORT_INC_PATH = $(ZOLTAN_INCLUDES)
#endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)

EXEEXT = .exe

noinst_PROGRAMS = zfdrive

if NAG_F90_COMPILER
FARG_SRC = $(srcdir)/farg_nagf95.f
else
FARG_SRC = $(srcdir)/farg_typical.f
endif

zfdrive_SOURCES = \
	$(srcdir)/mpi_h.f \
	$(FARG_SRC) \
	$(srcdir)/mmio.f \
	$(srcdir)/fdr_sort.f90 \
	$(srcdir)/fdr_const.f90 \
	$(srcdir)/fdr_input.f90 \
	$(srcdir)/fdr_chaco_io.f90 \
	$(srcdir)/fdr_param_file.f90 \
	$(srcdir)/fdr_mm_io.f90 \
	$(srcdir)/fdr_migrate.f90 \
	$(srcdir)/fdr_loadbal.f90 \
	$(srcdir)/fdr_main.f90 

zfdrive_DEPENDENCIES = \
			$(ZOLTAN_DEPS)

zfdrive_CPPFLAGS = $(EXPORT_INC_PATH) -I$(srcdir)/ -I$(top_srcdir)/src/ch/

zfdrive_LDADD = $(EXPORT_LIBS)

EXTRA_DIST = $(srcdir)/CMakeLists.txt \
	$(srcdir)/Makefile.am \
	$(srcdir)/Makefile.in \
	$(srcdir)/README.mpich \
	$(srcdir)/farg_nagf95.f \
	$(srcdir)/farg_typical.f \
	$(srcdir)/makefile \
	$(srcdir)/zoltan_user_data.f90.old

MOSTLYCLEANFILES = *.mod

.NOTPARALLEL:
