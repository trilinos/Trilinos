# @HEADER
# ************************************************************************
# 
#  For copyright information, ask xiaoye@nersc.gov 
# 
# ************************************************************************
# @HEADER

#np# Note the following paragraph.
# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
#np# There is, at present, no distinction between most of the categories, 
#np# other than those which start with EXTRA_.  All categories which do not
#np# start with EXTRA_ are included in all builds.  Categories which start
#np# with EXTRA_ are included in some builds (typically controlled by 
#np# --with-xyz or --without-xyz)
#
#np# Make sure to list all source files in one of the following categories.
CORE_H = \
	$(srcdir)/Cnames.h \
	$(srcdir)/f2c.h \
	$(srcdir)/scomplex.h \
	$(srcdir)/superlu_f2c.h

CORE = \
	$(srcdir)/caxpy.c \
	$(srcdir)/ccopy.c \
	$(srcdir)/cdotc.c \
	$(srcdir)/cgemv.c \
	$(srcdir)/cgerc.c \
	$(srcdir)/cgeru.c \
	$(srcdir)/chemv.c \
	$(srcdir)/cher2.c \
	$(srcdir)/cmyblas2.c \
	$(srcdir)/cscal.c \
	$(srcdir)/ctrsv.c \
	$(srcdir)/dasum.c \
	$(srcdir)/daxpy.c \
	$(srcdir)/dcabs1.c \
	$(srcdir)/dcopy.c \
	$(srcdir)/ddot.c \
	$(srcdir)/dgemm.c \
	$(srcdir)/dgemv.c \
	$(srcdir)/dger.c \
	$(srcdir)/dmyblas2.c \
	$(srcdir)/dnrm2.c \
	$(srcdir)/drot.c \
	$(srcdir)/dscal.c \
	$(srcdir)/dsymv.c \
	$(srcdir)/dsyr2.c \
	$(srcdir)/dtrsm.c \
	$(srcdir)/dtrsv.c \
	$(srcdir)/dzasum.c \
	$(srcdir)/dznrm2.c \
	$(srcdir)/icamax.c \
	$(srcdir)/idamax.c \
	$(srcdir)/isamax.c \
	$(srcdir)/izamax.c \
	$(srcdir)/sasum.c \
	$(srcdir)/saxpy.c \
	$(srcdir)/scasum.c \
	$(srcdir)/scnrm2.c \
	$(srcdir)/scopy.c \
	$(srcdir)/sdot.c \
	$(srcdir)/sgemv.c \
	$(srcdir)/sger.c \
	$(srcdir)/smyblas2.c \
	$(srcdir)/snrm2.c \
	$(srcdir)/srot.c \
	$(srcdir)/sscal.c \
	$(srcdir)/ssymv.c \
	$(srcdir)/ssyr2.c \
	$(srcdir)/strsv.c \
	$(srcdir)/zaxpy.c \
	$(srcdir)/zcopy.c \
	$(srcdir)/zdotc.c \
	$(srcdir)/zgemm.c \
	$(srcdir)/zgemv.c \
	$(srcdir)/zgerc.c \
	$(srcdir)/zgeru.c \
	$(srcdir)/zhemv.c \
	$(srcdir)/zher2.c \
	$(srcdir)/zmyblas2.c \
	$(srcdir)/zscal.c \
	$(srcdir)/ztrsm.c \
	$(srcdir)/ztrsv.c

UTIL_H = 

UTIL = 

ENUM_H = 

UNIPROC_H = 

UNIPROC = 

SERIAL_H = 

SERIAL = 

EXTRA_MPI_H =

EXTRA_MPI = 

if HAVE_MPI 
MPI_H = \
	$(EXTRA_MPI_H)
MPI = \
	$(EXTRA_MPI)

else

MPI_H =
MPI =

endif


USERBASE_H = 


LANG_WRAPS_H = 

LANG_WRAPS = 

DIST_H = 

DIST = 

COMP_H =

COMP = 


AM_CPPFLAGS = 

# ------------------------------------------------------------------------
# Cblas library specifications
# ------------------------------------------------------------------------
#np# replace cblas with the name of the package being autotool'ed here
CBLAS_LIB = libcblas.a

#np# replace cblas with the name of the package being autotool'ed here
CBLAS_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(ENUM_H) \
	$(UNIPROC_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) 

#np# replace cblas with the name of the package being autotool'ed here
libcblas_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(UNIPROC) \
	$(SERIAL) \
	$(DIST) \
	$(COMP) \
	$(MPI) 

#np# replace cblas with the name of the package being autotool'ed here
EXTRA_libcblas_a_SOURCES = \
	$(EXTRA_MPI) \
	$(EXTRA_MPI_H) 

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

#np# replace cblas with the name of the package being autotool'ed here
libcblas_a_AR = $(ALTERNATE_AR)
else

#np# replace cblas with the name of the package being autotool'ed here
libcblas_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

#np# replace cblas with the name of the package being autotool'ed here
libcblas_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#np# replace cblas with the name of the package being autotool'ed here
lib_LIBRARIES = $(CBLAS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#np# replace cblas with the name of the package being autotool'ed here
include_HEADERS = $(CBLAS_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

#np# replace cblas with the name of the package being autotool'ed here
# SPECIAL NOTE: Cblas_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
#
#  The cblas are not dependent on any compiler features, so 
#  we don't have to include Superludist_config.h
#np# replace cblas with the name of the package being autotool'ed here
nodist_execinclude_HEADERS = 

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






