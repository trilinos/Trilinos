


check_PROGRAMS=stressTestRIB stressTestRCB stressTestPHG stressTestGRAPH

include $(top_builddir)/Makefile.export.zoltan

LDADD=../../src/libzoltan.a -lm

EXPORT_INCLUDES = $(ZOLTAN_LIBRARY_INCLUDES)

AM_CPPFLAGS = $(EXPORT_INCLUDES) @MPI_RECV_LIMIT_FLAG@

