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
	$(srcdir)/Cnames.h  \
	$(srcdir)/dcomplex.h \
	$(srcdir)/superlu_ddefs.h \
	$(srcdir)/superlu_zdefs.h \
	$(srcdir)/util.h \
	$(srcdir)/colamd.h  \
	$(srcdir)/machines.h  \
	$(srcdir)/superlu_defs.h   \
	$(srcdir)/supermatrix.h

CORE = \
	$(srcdir)/colamd.c \
	$(srcdir)/comm.c \
	$(srcdir)/dcomplex.c \
	$(srcdir)/ddistribute.c \
	$(srcdir)/dgsequ.c \
	$(srcdir)/dlamch.c \
	$(srcdir)/dlangs.c \
	$(srcdir)/dlaqgs.c \
	$(srcdir)/dldperm.c \
	$(srcdir)/dmemory.c \
	$(srcdir)/dreadhb.c \
	$(srcdir)/dsp_blas2.c \
	$(srcdir)/dsp_blas3.c \
	$(srcdir)/dutil.c \
	$(srcdir)/etree.c \
	$(srcdir)/get_perm_c.c \
	$(srcdir)/lsame.c \
	$(srcdir)/memory.c \
	$(srcdir)/mmd.c \
	$(srcdir)/pddistribute.c \
	$(srcdir)/pdgsequ.c \
	$(srcdir)/pdgsmv_AXglobal.c \
	$(srcdir)/pdgsmv.c \
	$(srcdir)/pdgsrfs_ABXglobal.c \
	$(srcdir)/pdgsrfs.c \
	$(srcdir)/pdgssvx_ABglobal.c \
	$(srcdir)/pdgssvx.c \
	$(srcdir)/pdgstrf_irecv.c \
	$(srcdir)/pdgstrs1.c \
	$(srcdir)/pdgstrs_Bglobal.c \
	$(srcdir)/pdgstrs.c \
	$(srcdir)/pdgstrs_lsum.c \
	$(srcdir)/pdlangs.c \
	$(srcdir)/pdlaqgs.c \
	$(srcdir)/pdutil.c \
	$(srcdir)/pxerbla.c \
	$(srcdir)/pzdistribute.c \
	$(srcdir)/pzgsequ.c \
	$(srcdir)/pzgsmv_AXglobal.c \
	$(srcdir)/pzgsmv.c \
	$(srcdir)/pzgsrfs_ABXglobal.c \
	$(srcdir)/pzgsrfs.c \
	$(srcdir)/pzgssvx_ABglobal.c \
	$(srcdir)/pzgssvx.c \
	$(srcdir)/pzgstrf_irecv.c \
	$(srcdir)/pzgstrs1.c \
	$(srcdir)/pzgstrs_Bglobal.c \
	$(srcdir)/pzgstrs.c \
	$(srcdir)/pzgstrs_lsum.c \
	$(srcdir)/pzlangs.c \
	$(srcdir)/pzlaqgs.c \
	$(srcdir)/pzutil.c \
	$(srcdir)/slamch.c \
	$(srcdir)/sp_colorder.c \
	$(srcdir)/sp_ienv.c \
	$(srcdir)/superlu_grid.c \
	$(srcdir)/superlu_timer.c \
	$(srcdir)/symbfact.c \
	$(srcdir)/util.c \
	$(srcdir)/xerbla.c \
	$(srcdir)/zdistribute.c \
	$(srcdir)/zdistribute_mark.c \
	$(srcdir)/zgsequ.c \
	$(srcdir)/zlangs.c \
	$(srcdir)/zlaqgs.c \
	$(srcdir)/zldperm.c \
	$(srcdir)/zmemory.c \
	$(srcdir)/zreadhb.c \
	$(srcdir)/zsp_blas2.c \
	$(srcdir)/zsp_blas3.c \
	$(srcdir)/mc21.f \
	$(srcdir)/mc64ad.f \
	$(srcdir)/zutil.c 

UTIL_H = 

UTIL = 

ENUM_H = 

UNIPROC_H = 

UNIPROC = 

SERIAL_H = 

SERIAL = 

#
#np#  Files which are included in select builds are typically have
#np#  names that start with EXTRA
#
#np#  EXTRA_MPI is control by --with-mpi and --without-mpi which is
#np#  set up by TAC_ARC_CONFIG_MPI 
#
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


#np# List locations of include files within Trilinos but outside of this 
#np# package
#np# - Add packages that this package is dependent on
#np#   Remove epetra if this package does not need it
#  
AM_CPPFLAGS = 

#np# The following section can be duplicated for building multiple libraries -
#np# make sure to list the correct categories with each library.
#np# replace new_package with the name of the package being autotool'ed here
# ------------------------------------------------------------------------
# Superludist library specifications
# ------------------------------------------------------------------------
#np# replace superludist with the name of the package being autotool'ed here
SUPERLUDIST_LIB = libsuperludist.a

#np# replace superludist with the name of the package being autotool'ed here
SUPERLUDIST_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(ENUM_H) \
	$(UNIPROC_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) 

#np# replace superludist with the name of the package being autotool'ed here
libsuperludist_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(UNIPROC) \
	$(SERIAL) \
	$(DIST) \
	$(COMP) \
	$(MPI) 

#np# replace superludist with the name of the package being autotool'ed here
EXTRA_libsuperludist_a_SOURCES = \
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

#np# replace superludist with the name of the package being autotool'ed here
libsuperludist_a_AR = $(ALTERNATE_AR)
else

#np# replace superludist with the name of the package being autotool'ed here
libsuperludist_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

#np# replace superludist with the name of the package being autotool'ed here
libsuperludist_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#np# replace superludist with the name of the package being autotool'ed here
lib_LIBRARIES = $(SUPERLUDIST_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#np# replace superludist with the name of the package being autotool'ed here
include_HEADERS = $(SUPERLUDIST_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

#np# replace superludist with the name of the package being autotool'ed here
# SPECIAL NOTE: Superludist_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
#np# replace superludist with the name of the package being autotool'ed here
nodist_execinclude_HEADERS = Superludist_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






