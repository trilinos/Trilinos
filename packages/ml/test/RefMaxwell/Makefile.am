# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

EXTRA_DIST = M0.dat.gz M1.dat.gz coord_node.dat.gz Ms.dat.gz S.dat.gz Tclean.dat.gz

noinst_PROGRAMS = RefMaxwell GetNecessaryFiles

RefMaxwell_SOURCES      = $(srcdir)/cxx_main.cpp
RefMaxwell_DEPENDENCIES = $(top_builddir)/src/libml.a
RefMaxwell_CXXFLAGS     = ${TEST_INCLUDES}
RefMaxwell_LDADD        = ${TEST_LIBS}

GetNecessaryFiles_SOURCES= $(top_srcdir)/test/RefMaxwell

#copy these files over to the build area
GetNecessaryFiles$(EXEEXT): $(GetNecessaryFiles_SOURCES)
	@if test ! -e ./M0.dat; then cp -f $(srcdir)/M0.dat.gz .; gunzip -f ./M0.dat.gz ; fi
	@if test ! -e ./M1.dat; then cp -f $(srcdir)/M1.dat.gz .; gunzip -f ./M1.dat.gz ; fi
	@if test ! -e ./Ms.dat; then cp -f $(srcdir)/Ms.dat.gz .; gunzip -f ./Ms.dat.gz ; fi
	@if test ! -e ./S.dat ; then cp -f $(srcdir)/S.dat.gz  .; gunzip -f ./S.dat.gz  ; fi
	@if test ! -e ./Tclean.dat    ; then cp -f $(srcdir)/Tclean.dat.gz     .; gunzip -f Tclean.dat.gz     ; fi
	@if test ! -e ./coord_node.dat; then cp -f $(srcdir)/coord_node.dat.gz .; gunzip -f coord_node.dat.gz ; fi

