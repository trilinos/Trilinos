include $(top_builddir)/Makefile.export.stratimikos

_EXPORT_INC_PATH = $(STRATIMIKOS_TEST_INCLUDES)
_EXPORT_LIBS = $(STRATIMIKOS_TEST_LIBS)

if USING_GNUMAKE
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(_EXPORT_INC_PATH))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(_EXPORT_LIBS))
else
EXPORT_INC_PATH = $(_EXPORT_INC_PATH)
EXPORT_LIBS = $(_EXPORT_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)

EXEEXT = .exe

#
# Actual programs
#

if ENABLE_STRATIMIKOS_THYRA_EPETRAEXT
TEST_SINGLE_BELOS_THYRA_SOLVER_DRIVER = test_single_belos_thyra_solver_driver
#TEST_BELOS_THYRA_DRIVER = test_belos_thyra_driver
endif

noinst_PROGRAMS = $(TEST_SINGLE_BELOS_THYRA_SOLVER_DRIVER) $(TEST_BELOS_THYRA_DRIVER)

test_single_belos_thyra_solver_driver_SOURCES = $(srcdir)/test_single_belos_thyra_solver_driver.cpp $(common_sources)
test_single_belos_thyra_solver_driver_DEPENDENCIES = $(common_dependencies)
test_single_belos_thyra_solver_driver_LDADD = $(common_ldadd)

#test_belos_thyra_driver_SOURCES = $(srcdir)/test_belos_thyra_driver.cpp $(common_sources)
#test_belos_thyra_driver_DEPENDENCIES = $(common_dependencies) matrices
#test_belos_thyra_driver_LDADD = $(common_ldadd)

FourByFour.mtx:
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_bcsstk01.mtx bcsstk01.mtx                              
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_bcsstk02.mtx bcsstk02.mtx                              
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_bcsstk04.mtx bcsstk04.mtx                              
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_nos1.mtx nos1.mtx                                      
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_nos5.mtx nos5.mtx                                      
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_FourByFour.mtx FourByFour.mtx                          
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_Diagonal.mtx Diagonal.mtx                              
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_MissingADiagonal.mtx MissingADiagonal.mtx                              
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_KheadK.mtx KheadK.mtx                                  
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_KheadSorted.mtx  KheadSorted.mtx
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_Khead.triS Khead.triS  

#
# Common stuff
#

common_sources = $(srcdir)/test_single_belos_thyra_solver.cpp

common_dependencies = $(top_builddir)/../belos/src/libbelos.a FourByFour.mtx

common_ldadd = $(EXPORT_LIBS)

EXTRA_DIST = $(srcdir)/test_single_belos_thyra_solver.hpp 

# This is already added part of BELOS_LIBS and therefore automake does not
# need to add this again!
LIBS =
