
include $(top_builddir)/Makefile.export.zoltan

EXPORT_LIBS = $(ZOLTAN_LIBS)
EXPORT_INC_PATH = $(ZOLTAN_INCLUDES)

AM_CPPFLAGS = $(EXPORT_INC_PATH)

EXEEXT = .exe

#ZOLTAN_C_TEST_DRIVER = zdrive 
if HAVE_ZOLTAN_CPPDRIVER
  ZOLTAN_CPP_TEST_DRIVER = zCPPdrive 
endif

noinst_PROGRAMS = zdrive $(ZOLTAN_CPP_TEST_DRIVER)

COMMON_DRIVER_SRC = \
	$(srcdir)/dr_err.c $(srcdir)/dr_input.c $(srcdir)/dr_par_util.c \
	$(srcdir)/dr_util.c $(srcdir)/dr_elem.c $(srcdir)/dr_setfixed.c \
	$(srcdir)/dr_output.c $(srcdir)/dr_chaco_io.c \
	$(srcdir)/dr_random_io.c $(srcdir)/dr_elem_util.c $(srcdir)/dr_eval.c \
	$(srcdir)/dr_gnuplot.c $(srcdir)/dr_hg_io.c \
	$(srcdir)/dr_hg_readfile.c $(srcdir)/dr_mm_readfile.c \
	$(srcdir)/dr_mmio.c \
	$(top_srcdir)/src/ch/ch_dist_graph.c \
	$(top_srcdir)/src/ch/ch_input_assign.c  \
	$(top_srcdir)/src/ch/ch_input_graph.c \
	$(top_srcdir)/src/ch/ch_init_dist.c \
	$(top_srcdir)/src/ch/ch_input_geom.c \
	$(top_srcdir)/src/ch/ch_input_read.c

if BUILD_GZIP
COMMON_DRIVER_SRC += $(srcdir)/dr_compress.c
endif

# List headers here that are used for both the C and C++ drivers.
COMMON_DRIVER_HDRS = \
	$(srcdir)/dr_compress_const.h $(srcdir)/dr_hg_readfile.h \
	$(srcdir)/dr_externs.h $(srcdir)/dr_const.h $(srcdir)/dr_input_const.h \
	$(srcdir)/dr_util_const.h $(srcdir)/dr_err_const.h \
	$(srcdir)/dr_par_util_const.h $(srcdir)/dr_elem_const.h \
	$(srcdir)/dr_elem_util_const.h $(srcdir)/dr_maps_const.h \
	$(srcdir)/dr_eval_const.h $(srcdir)/dr_loadbal_const.h \
	$(srcdir)/dr_mmio.h $(srcdir)/dr_output_const.h \
	$(srcdir)/dr_param_file.h $(srcdir)/dr_dd.h \
	$(top_srcdir)/src/ch/ch_input_const.h \
	$(top_srcdir)/src/ch/ch_init_dist_const.h

zdrive_SOURCES = \
	$(COMMON_DRIVER_SRC) $(srcdir)/dr_main.c $(srcdir)/dr_exoII_io.c \
	$(srcdir)/dr_loadbal.c $(srcdir)/dr_migrate.c $(srcdir)/dr_maps.c \
	$(srcdir)/dr_dd.c $(srcdir)/dr_param_file.c $(srcdir)/order_test.c

# List headers here used for the C driver, but not the C++ driver.
# Headers for both are listed under COMMON_DRIVER_HDRS.
#zdrive_HDRS =
zCPPdrive_SOURCES = \
	$(COMMON_DRIVER_SRC) $(srcdir)/dr_mainCPP.cpp \
	$(srcdir)/dr_exoII_ioCPP.cpp $(srcdir)/dr_loadbalCPP.cpp \
	$(srcdir)/dr_migrateCPP.cpp $(srcdir)/dr_mapsCPP.cpp \
	$(srcdir)/dr_ddCPP.cpp $(srcdir)/dr_param_fileCPP.cpp

# List headers here used for the C++ driver, but not the C driver.
# Headers for both are listed under COMMON_DRIVER_HDRS.
zCPPdrive_HDRS = \
	$(srcdir)/dr_param_fileCPP.h

zdrive_DEPENDENCIES = \
			$(ZOLTAN_DEPS)

zCPPdrive_DEPENDENCIES = \
			$(ZOLTAN_DEPS)

zdrive_CPPFLAGS = $(EXPORT_INC_PATH) -I$(srcdir)/ -I$(top_srcdir)/src/ch/

zCPPdrive_CPPFLAGS = $(EXPORT_INC_PATH) -I$(srcdir)/ -I$(top_srcdir)/src/ch/

zdrive_LDADD = $(EXPORT_LIBS)

zCPPdrive_LDADD = $(EXPORT_LIBS)

EXTRA_DIST = $(COMMON_DRIVER_HDRS) $(zCPPdrive_HDRS) $(zdrive_HDRS) \
             $(srcdir)/dr_compress.c

# Dummy for new export makefiles
#LIBS = 

