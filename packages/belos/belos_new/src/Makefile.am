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
	$(srcdir)/BelosBlockCG.hpp \
	$(srcdir)/BelosBlockGmres.hpp \
	$(srcdir)/BelosBlockGmresIter.hpp \
	$(srcdir)/BelosBlockGmresSolMgr.hpp \
	$(srcdir)/BelosCG.hpp \
	$(srcdir)/BelosConfigDefs.hpp \
	$(srcdir)/BelosDGKSOrthoManager.hpp \
	$(srcdir)/BelosEpetraAdapter.hpp \
	$(srcdir)/BelosEpetraOperator.h \
	$(srcdir)/BelosICGSOrthoManager.hpp \
	$(srcdir)/BelosIteration.hpp \
	$(srcdir)/BelosIterativeSolver.hpp \
	$(srcdir)/BelosLinearProblem.hpp \
	$(srcdir)/BelosMatOrthoManager.hpp \
	$(srcdir)/BelosMultiVec.hpp \
	$(srcdir)/BelosMultiVecTraits.hpp \
	$(srcdir)/BelosOperator.hpp \
	$(srcdir)/BelosOperatorTraits.hpp \
	$(srcdir)/BelosOrthoManager.hpp \
	$(srcdir)/BelosOutputManager.hpp \
	$(srcdir)/BelosPseudoBlockGmres.hpp \
	$(srcdir)/BelosStatusTest.hpp \
	$(srcdir)/BelosStatusTestCombo.hpp \
	$(srcdir)/BelosStatusTestMaxIters.hpp \
	$(srcdir)/BelosStatusTestMaxRestarts.hpp \
	$(srcdir)/BelosStatusTestOutputter.hpp \
	$(srcdir)/BelosStatusTestResNorm.hpp \
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
	$(srcdir)/BelosBlockCG.hpp \
	$(srcdir)/BelosBlockGmres.hpp \
        $(srcdir)/BelosCG.hpp \
        $(srcdir)/BelosPseudoBlockGmres.hpp

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

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = Belos_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 


# ------------------------------------------------------------------------
# Library Target for Circular Dependencies
# ------------------------------------------------------------------------
if BUILD_IFPACK
BUILD_IFPACK_LIBS = cd $(top_builddir)/../ifpack ; $(MAKE) all-libs
else
BUILD_IFPACK_LIBS =
endif

all-libs:
	$(MAKE) libbelos.a
	cd $(top_builddir)/../epetra ; $(MAKE) all-libs
	cd $(top_builddir)/../triutils ; $(MAKE) all-libs
	cd $(top_builddir)/../teuchos ; $(MAKE) all-libs
	$(BUILD_IFPACK_LIBS)
#	$(BUILD_TSFCORE_EPETRA_LIBS)
