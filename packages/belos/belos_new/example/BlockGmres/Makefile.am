#@HEADER
#@HEADER

include $(top_builddir)/Makefile.export.belos

if USING_GNUMAKE
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(BELOS_TEST_LIBS))
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(BELOS_TEST_INCLUDES))
else
EXPORT_LIBS = $(BELOS_TEST_LIBS)
EXPORT_INC_PATH = $(BELOS_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)

EXEEXT = .exe

if BUILD_TRIUTILS
#if BUILD_IFPACK
#BELOS_BGMRES_IFPACK_EX = bl_pgmres_hb bl_fgmres_hb
#endif
if BUILD_EPETRAEXT
BELOS_BGMRES_NOPREC_EX = BlockGmresEpetraExFile Dummy 
endif
endif

if BUILD_TRIUTILS
if BUILD_TEUCHOS_COMPLEX
BELOS_BGMRES_COMPLEX_EX = bl_gmres_complex_hb ComplexDummy
endif
endif


#noinst_PROGRAMS = $(BELOS_BGMRES_EX) $(BELOS_BGMRES_IFPACK_EX) $(BELOS_BGMRES_NOPREC_EX) $(BELOS_BGMRES_COMPLEX_EX)
noinst_PROGRAMS = $(BELOS_BGMRES_NOPREC_EX)

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(srcdir)/orsirr1.hb orsirr1.hb
	cp  $(srcdir)/fidap036.hb fidap036.hb

Dummy_SOURCES = Dummy

#ComplexDummy$(EXEEXT): $(srcdir)/Dummy
#	cp  $(top_srcdir)/test/MVOPTester/mhd1280b.cua mhd1280b.cua
#	cp  $(top_srcdir)/test/MVOPTester/mhd1280a.cua mhd1280a.cua

#ComplexDummy_SOURCES = Dummy

#bl_pgmres_hb_SOURCES = $(srcdir)/test_bl_pgmres_hb.cpp $(srcdir)/createEpetraProblem.cpp

#bl_fgmres_hb_SOURCES = $(srcdir)/test_bl_fgmres_hb.cpp $(srcdir)/createEpetraProblem.cpp

BlockGmresEpetraExFile_SOURCES = $(srcdir)/BlockGmresEpetraExFile.cpp 

#bl_gmres_complex_hb_SOURCES = $(srcdir)/test_bl_gmres_complex_hb.cpp

#bl_pgmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a			

#bl_fgmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a			

BlockGmresEpetraExFile_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

#bl_gmres_complex_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
#			$(top_builddir)/../epetra/src/libepetra.a \
#			$(top_builddir)/../teuchos/src/libteuchos.a

#bl_pgmres_hb_CXXFLAGS = -I$(top_srcdir)/src

#bl_fgmres_hb_CXXFLAGS = -I$(top_srcdir)/src 

BlockGmresEpetraExFile_CXXFLAGS = -I$(top_srcdir)/src 

#bl_gmres_complex_hb_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/test/MVOPTester

#bl_pgmres_hb_LDADD = $(EXPORT_LIBS)

#bl_fgmres_hb_LDADD = $(EXPORT_LIBS)

BlockGmresEpetraExFile_LDADD = $(EXPORT_LIBS)

#bl_gmres_complex_hb_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 
#
