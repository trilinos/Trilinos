# @HEADER
# ************************************************************************
# 
#	      Finite Element Interface to Linear Solvers
#		 Copyright (2001) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# Questions? Contact Alan Williams (william@sandia.gov) 
# 
# ************************************************************************
# @HEADER

include $(top_builddir)/Makefile.export.fei

if FEI_HAVE_HYPRE
HYPRE_INCLUDES=@HYPRE_INCLUDE@
HYPRE_CPPFLAGS=$(HYPRE_INCLUDES)
else
HYPRE_INCLUDES=
HYPRE_CPPFLAGS=
endif

if FEI_HAVE_PETSC
PETSC_INCLUDES=@PETSC_INCLUDE@
PETSC_CPPFLAGS=$(PETSC_INCLUDES)
else
PETSC_INCLUDES=
PETSC_CPPFLAGS=
endif

if FEI_HAVE_PROMETHEUS
PROMETHEUS_INCLUDES=@PROMETHEUS_INCLUDE@
PROMETHEUS_CPPFLAGS=$(PROMETHEUS_INCLUDES)
else
PROMETHEUS_INCLUDES=
PROMETHEUS_CPPFLAGS=
endif

if FEI_HAVE_FETI
FETIDP_INCLUDES=@FETIDP_INCLUDE@
FETIDP_CPPFLAGS=$(FETIDP_INCLUDES)
else
FETIDP_INCLUDES=
FETIDP_CPPFLAGS=
endif

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside FEI. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

DEFS =
DEFAULT_INCLUDES = -I$(top_srcdir)/utest_cases $(FEI_INCLUDES) $(HYPRE_CPPFLAGS) $(PETSC_CPPFLAGS) $(PROMETHEUS_CPPFLAGS) $(FETIDP_CPPFLAGS)

CORE_H = \
	$(top_srcdir)/utest_cases/fei_unit_CommUtils.hpp \
	$(top_srcdir)/utest_cases/fei_unit_CSRMat_CSVec.hpp \
	$(top_srcdir)/utest_cases/fei_unit_DirBC.hpp \
	$(top_srcdir)/utest_cases/fei_unit_Reducer.hpp \
	$(top_srcdir)/utest_cases/fei_unit_testcase.hpp \
	$(top_srcdir)/utest_cases/fei_unit_testcontainer.hpp \
	$(top_srcdir)/utest_cases/fei_unit_register_tests.hpp \
	$(top_srcdir)/utest_cases/fei_unit_test_runner.hpp \
	$(top_srcdir)/utest_cases/fei_unit_utils.hpp

CORE = \
	$(srcdir)/fei_unit_testcase.cpp \
	$(srcdir)/fei_unit_testcontainer.cpp \
	$(srcdir)/fei_unit_register_tests.cpp \
	$(srcdir)/fei_unit_test_runner.cpp \
	$(srcdir)/fei_unit_CSRMat_CSVec.cpp \
	$(srcdir)/fei_unit_utils.cpp \
	$(srcdir)/fei_unit_CommUtils.cpp \
	$(srcdir)/fei_unit_DirBC.cpp \
	$(srcdir)/fei_unit_Reducer.cpp

# ------------------------------------------------------------------------
# FEI library specifications
# ------------------------------------------------------------------------

FEI_LIB = libfei_utest_cases.a

noinst_LIBRARIES = $(FEI_LIB)


libfei_utest_cases_a_SOURCES = \
	$(CORE_H) $(CORE) 

# EXTRA_libfei_utest_cases_a_SOURCES = 

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR if USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libfei_utest_cases_a_AR = $(ALTERNATE_AR)
else

libfei_utest_cases_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libfei_utest_cases_a_LIBADD = $(FEI_TEMPLATE_OBJS) $(XTRALDADD)

includedir = ${prefix}/utest_cases

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

