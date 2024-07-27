
include $(top_builddir)/Makefile.export.pamgen

# The following line helps the test harness recover from build errors.

all-local: 

EXEEXT = .exe

noinst_PROGRAMS = pamgen_test_driver 

pamgen_test_driver_SOURCES = $(srcdir)/../example/pamgen_test_driver.c $(srcdir)/../example/getopts.C

#np# Specify locations of include files within Trilinos, but outside np
#np# Add packages that np is dependent on; remove epetra if np does not need it

AM_CPPFLAGS = $(PAMGEN_INCLUDES)

#np# - List libraries that np is dependent on.  The DEPENDENCIES primary is
#np# optional.
pamgen_test_driver_DEPENDENCIES = \
	$(top_builddir)/src/libpamgen.a 

#np# List libraries that np is dependent on.
pamgen_test_driver_LDADD = $(PAMGEN_LIBS)

#
#  This is a crude but effective way to move files from the source tree to the build tree
#  We recommend renaming files because, although we recommend against the practice, some
#  users will configure and build in the same tree.
#
#  This requires at least one file to end in .exe (or whatever EXEEXT is set to).
#  The other files can have arbitrary names. (In this case Suppressions)
#

pamgen_test_driver$(EXEEXT): $(srcdir)/test_files/run_test_script_01
	cp ../example/pamgen_test_driver.exe .
	cp $(srcdir)/test_files/run_test_script_01 run_test_script_01$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_02 run_test_script_02$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_03 run_test_script_03$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_04 run_test_script_04$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_05 run_test_script_05$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_06 run_test_script_06$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_07 run_test_script_07$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_08 run_test_script_08$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_09 run_test_script_09$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_10 run_test_script_10$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_11 run_test_script_11$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_12 run_test_script_12$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_13 run_test_script_13$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_14 run_test_script_14$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_15 run_test_script_15$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_16 run_test_script_16$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_17 run_test_script_17$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_18 run_test_script_18$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_19 run_test_script_19$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_20 run_test_script_20$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_21 run_test_script_21$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_22 run_test_script_22$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_23 run_test_script_23$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_24 run_test_script_24$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_25 run_test_script_25$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_26 run_test_script_26$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_27 run_test_script_27$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_28 run_test_script_28$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_29 run_test_script_29$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_30 run_test_script_30$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_31 run_test_script_31$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_32 run_test_script_32$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_33 run_test_script_33$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_34 run_test_script_34$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_35 run_test_script_35$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_36 run_test_script_36$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_37 run_test_script_37$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_38 run_test_script_38$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_39 run_test_script_39$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_40 run_test_script_40$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_41 run_test_script_41$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_42 run_test_script_42$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_43 run_test_script_43$(EXEEXT)
	cp $(srcdir)/test_files/run_test_script_44 run_test_script_44$(EXEEXT)
	cp $(srcdir)/test_files/*.txt .
	cp $(srcdir)/test_files/*.res .

#np# Set "SUBDIRS= (list all subdirectories that are autotool'ed here)"
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.
#np# For packages whose examples are located in this directory, rather than
#np# in a subdirectory of this directory, use the Makefile.am located in
#np# example1/ rather than this file.
# SUBDIRS = example1
#SUBDIRS = 

EXTRA_DIST = \
$(srcdir)/test_files/run_test_script_01 \
$(srcdir)/test_files/run_test_script_02 \
$(srcdir)/test_files/run_test_script_03 \
$(srcdir)/test_files/run_test_script_04 \
$(srcdir)/test_files/run_test_script_05 \
$(srcdir)/test_files/run_test_script_06 \
$(srcdir)/test_files/run_test_script_07 \
$(srcdir)/test_files/run_test_script_08 \
$(srcdir)/test_files/run_test_script_09 \
$(srcdir)/test_files/run_test_script_10 \
$(srcdir)/test_files/run_test_script_11 \
$(srcdir)/test_files/run_test_script_12 \
$(srcdir)/test_files/run_test_script_13 \
$(srcdir)/test_files/run_test_script_14 \
$(srcdir)/test_files/run_test_script_15 \
$(srcdir)/test_files/run_test_script_16 \
$(srcdir)/test_files/run_test_script_17 \
$(srcdir)/test_files/run_test_script_18 \
$(srcdir)/test_files/run_test_script_19 \
$(srcdir)/test_files/run_test_script_20 \
$(srcdir)/test_files/run_test_script_21 \
$(srcdir)/test_files/run_test_script_22 \
$(srcdir)/test_files/run_test_script_23 \
$(srcdir)/test_files/run_test_script_24 \
$(srcdir)/test_files/run_test_script_25 \
$(srcdir)/test_files/run_test_script_26 \
$(srcdir)/test_files/run_test_script_27 \
$(srcdir)/test_files/run_test_script_28 \
$(srcdir)/test_files/run_test_script_29 \
$(srcdir)/test_files/run_test_script_30 \
$(srcdir)/test_files/run_test_script_31 \
$(srcdir)/test_files/run_test_script_32 \
$(srcdir)/test_files/run_test_script_33 \
$(srcdir)/test_files/run_test_script_34 \
$(srcdir)/test_files/run_test_script_35 \
$(srcdir)/test_files/run_test_script_36 \
$(srcdir)/test_files/run_test_script_37 \
$(srcdir)/test_files/run_test_script_38 \
$(srcdir)/test_files/run_test_script_39 \
$(srcdir)/test_files/run_test_script_40 \
$(srcdir)/test_files/run_test_script_41 \
$(srcdir)/test_files/run_test_script_42 \
$(srcdir)/test_files/run_test_script_43 \
$(srcdir)/test_files/run_test_script_44 \
$(srcdir)/test_files/ur.txt \
$(srcdir)/test_files/uded.txt \
$(srcdir)/test_files/report3.txt \
$(srcdir)/test_files/trisect_multi_theta.txt \
$(srcdir)/test_files/ui2d.txt \
$(srcdir)/test_files/2dtribis.txt \
$(srcdir)/test_files/all_diffs.txt \
$(srcdir)/test_files/tr.txt \
$(srcdir)/test_files/report9.txt \
$(srcdir)/test_files/mesh.txt \
$(srcdir)/test_files/ok.txt \
$(srcdir)/test_files/report17.txt \
$(srcdir)/test_files/btri.txt \
$(srcdir)/test_files/report12.txt \
$(srcdir)/test_files/tri.txt \
$(srcdir)/test_files/foo.txt \
$(srcdir)/test_files/3dbj.txt \
$(srcdir)/test_files/m2.txt \
$(srcdir)/test_files/trisect_multi_theta_2d.txt \
$(srcdir)/test_files/p.txt \
$(srcdir)/test_files/report5.txt \
$(srcdir)/test_files/tp.txt \
$(srcdir)/test_files/report2.txt \
$(srcdir)/test_files/ff.txt \
$(srcdir)/test_files/tri_proc.txt \
$(srcdir)/test_files/report13.txt \
$(srcdir)/test_files/2dtri.txt \
$(srcdir)/test_files/m3.txt \
$(srcdir)/test_files/messed.txt \
$(srcdir)/test_files/report15.txt \
$(srcdir)/test_files/geom_warp2d.txt \
$(srcdir)/test_files/report14.txt \
$(srcdir)/test_files/brick.txt \
$(srcdir)/test_files/report10.txt \
$(srcdir)/test_files/hbrick.txt \
$(srcdir)/test_files/mesh7.txt \
$(srcdir)/test_files/report6.txt \
$(srcdir)/test_files/geom_warp.txt \
$(srcdir)/test_files/report1.txt \
$(srcdir)/test_files/report11.txt \
$(srcdir)/test_files/ui.txt \
$(srcdir)/test_files/report7.txt \
$(srcdir)/test_files/report16.txt \
$(srcdir)/test_files/report8.txt \
$(srcdir)/test_files/report4.txt \
$(srcdir)/test_files/run_test_script_29.res \
$(srcdir)/test_files/run_test_script_34.res \
$(srcdir)/test_files/run_test_script_19.res \
$(srcdir)/test_files/run_test_script_31.res \
$(srcdir)/test_files/run_test_script_04.res \
$(srcdir)/test_files/run_test_script_33.res \
$(srcdir)/test_files/run_test_script_18.res \
$(srcdir)/test_files/run_test_script_37.res \
$(srcdir)/test_files/run_test_script_13.res \
$(srcdir)/test_files/run_test_script_08.res \
$(srcdir)/test_files/run_test_script_41.res \
$(srcdir)/test_files/run_test_script_16.res \
$(srcdir)/test_files/run_test_script_36.res \
$(srcdir)/test_files/run_test_script_12.res \
$(srcdir)/test_files/run_test_script_11.res \
$(srcdir)/test_files/run_test_script_05.res \
$(srcdir)/test_files/run_test_script_32.res \
$(srcdir)/test_files/run_test_script_15.res \
$(srcdir)/test_files/run_test_script_40.res \
$(srcdir)/test_files/run_test_script_27.res \
$(srcdir)/test_files/run_test_script_06.res \
$(srcdir)/test_files/run_test_script_38.res \
$(srcdir)/test_files/run_test_script_30.res \
$(srcdir)/test_files/run_test_script_21.res \
$(srcdir)/test_files/run_test_script_28.res \
$(srcdir)/test_files/run_test_script_25.res \
$(srcdir)/test_files/run_test_script_42.res \
$(srcdir)/test_files/run_test_script_43.res \
$(srcdir)/test_files/run_test_script_44.res \
$(srcdir)/test_files/run_test_script_07.res \
$(srcdir)/test_files/run_test_script_20.res \
$(srcdir)/test_files/run_test_script_10.res \
$(srcdir)/test_files/run_test_script_09.res \
$(srcdir)/test_files/run_test_script_26.res \
$(srcdir)/test_files/run_test_script_03.res \
$(srcdir)/test_files/run_test_script_35.res \
$(srcdir)/test_files/run_test_script_17.res \
$(srcdir)/test_files/run_test_script_14.res \
$(srcdir)/test_files/run_test_script_23.res \
$(srcdir)/test_files/run_test_script_02.res \
$(srcdir)/test_files/run_test_script_24.res \
$(srcdir)/test_files/run_test_script_39.res \
$(srcdir)/test_files/run_test_script_22.res \
$(srcdir)/test_files/run_test_script_01.res

