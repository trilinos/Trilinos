# @HEADER
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package belos subdirectory src built successfully."
	@echo ""

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside Ifpack. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

CORE_H = \
	$(srcdir)/BelosBlockCGSolMgr.hpp \
	$(srcdir)/BelosBlockCGIter.hpp \
	$(srcdir)/BelosBlockGmresIter.hpp \
	$(srcdir)/BelosBlockFGmresIter.hpp \
	$(srcdir)/BelosBlockGmresSolMgr.hpp \
	$(srcdir)/BelosCGIter.hpp \
	$(srcdir)/BelosCGIteration.hpp \
	$(srcdir)/BelosConfigDefs.hpp \
	$(srcdir)/BelosDGKSOrthoManager.hpp \
	$(srcdir)/BelosEpetraAdapter.hpp \
	$(srcdir)/BelosEpetraOperator.h \
	$(srcdir)/BelosGCRODRIter.hpp \
	$(srcdir)/BelosGCRODRSolMgr.hpp \
        $(srcdir)/BelosGmresIteration.hpp \
	$(srcdir)/BelosGmresPolySolMgr.hpp \
	$(srcdir)/BelosICGSOrthoManager.hpp \
	$(srcdir)/BelosIMGSOrthoManager.hpp \
	$(srcdir)/BelosIteration.hpp \
	$(srcdir)/BelosLinearProblem.hpp \
	$(srcdir)/BelosMatOrthoManager.hpp \
	$(srcdir)/BelosMultiVec.hpp \
	$(srcdir)/BelosMultiVecTraits.hpp \
	$(srcdir)/BelosMVOPTester.hpp \
	$(srcdir)/BelosOperator.hpp \
	$(srcdir)/BelosOperatorTraits.hpp \
	$(srcdir)/BelosOrthoManager.hpp \
	$(srcdir)/BelosOutputManager.hpp \
	$(srcdir)/BelosPseudoBlockGmresIter.hpp \
	$(srcdir)/BelosPseudoBlockGmresSolMgr.hpp \
	$(srcdir)/BelosRCGIter.hpp \
	$(srcdir)/BelosRCGSolMgr.hpp \
	$(srcdir)/BelosSolverManager.hpp \
	$(srcdir)/BelosStatusTest.hpp \
	$(srcdir)/BelosStatusTestCombo.hpp \
	$(srcdir)/BelosStatusTestMaxIters.hpp \
	$(srcdir)/BelosStatusTestOutput.hpp \
	$(srcdir)/BelosStatusTestResNorm.hpp \
	$(srcdir)/BelosStatusTestGenResNorm.hpp \
	$(srcdir)/BelosStatusTestImpResNorm.hpp \
	$(srcdir)/BelosTypes.hpp

CORE = \
	$(srcdir)/BelosEpetraAdapter.cpp \
	$(srcdir)/BelosEpetraOperator.cpp

UTIL = \
	$(srcdir)/BelosVersion.cpp

### add in paths as needed

# ------------------------------------------------------------------------
# Includes for required and optional library dependencies
# ------------------------------------------------------------------------
include $(top_builddir)/Makefile.export.belos

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(BELOS_INCLUDES))
else
EXPORT_INCLUDES = $(BELOS_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)
# ------------------------------------------------------------------------

SERIAL_H = 

SERIAL = 

if HAVE_MPI

MPI_H = 

MPI = 

else

MPI_H =
MPI =

endif

#if HAVE_THREADS

MPISMP_H = 

MPISMP = 

#else

#MPISMP_H =
#MPISMP =

#endif

LANG_WRAPS_H = 

LANG_WRAPS = 

DIST_H = 

DIST = 

COMP_H = \
	$(srcdir)/BelosBlockCGSolMgr.hpp \
	$(srcdir)/BelosBlockGmresSolMgr.hpp \
        $(srcdir)/BelosPseudoBlockGmresSolMgr.hpp \
        $(srcdir)/BelosCGIter.hpp \
        $(srcdir)/BelosBlockCGIter.hpp \
        $(srcdir)/BelosBlockGmresIter.hpp \
        $(srcdir)/BelosPseudoBlockGmresIter.hpp

COMP = 



# ------------------------------------------------------------------------
# Belos library specifications
# ------------------------------------------------------------------------

BELOS_LIB = libbelos.a


BELOS_H = \
	$(CORE_H) \
	$(SERIAL_H) \
	$(LANG_WRAPS_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) \
	$(MPISMP_H)

libbelos_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(SERIAL) \
	$(LANG_WRAPS) \
	$(DIST) \
	$(COMP) \
	$(MPI) \
	$(MPISMP)

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libbelos_a_AR = $(ALTERNATE_AR)
else

libbelos_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libbelos_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(BELOS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(BELOS_H)

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Belos_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(includedir)
nodist_execinclude_HEADERS = Belos_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 
