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

if HAVE_FEI_PETSC
PETSC_INCLUDES=@PETSC_INCLUDE@
PETSC_CPPFLAGS=$(PETSC_INCLUDES)
else
PETSC_INCLUDES=
PETSC_CPPFLAGS=
endif

if HAVE_FEI_FETI
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
DEFAULT_INCLUDES = $(FEI_INCLUDES) -I@abs_top_srcdir@/test_utils -I@abs_top_srcdir@/support-Trilinos $(PETSC_CPPFLAGS) $(FETIDP_CPPFLAGS)

CORE_H = \
	$(top_srcdir)/test_utils/AccessPattern.hpp \
	$(top_srcdir)/test_utils/BCNodeSet.hpp \
	$(top_srcdir)/test_utils/CoefAccess.hpp \
	$(top_srcdir)/test_utils/CommNodeSet.hpp \
	$(top_srcdir)/test_utils/CRSet.hpp \
	$(top_srcdir)/test_utils/DataReader.hpp \
	$(top_srcdir)/test_utils/driverData.hpp \
	$(top_srcdir)/test_utils/ElemBlock.hpp \
	$(top_srcdir)/test_utils/FEData.hpp \
	$(top_srcdir)/test_utils/FEI_tester.hpp \
	$(top_srcdir)/test_utils/feitester.hpp \
	$(top_srcdir)/test_utils/fei_test_utils.hpp \
	$(top_srcdir)/test_utils/HexBeamCR.hpp \
	$(top_srcdir)/test_utils/HexBeam.hpp \
	$(top_srcdir)/test_utils/InputData.hpp \
	$(top_srcdir)/test_utils/InputFileReader.hpp \
	$(top_srcdir)/test_utils/LibraryFactory.hpp \
	$(top_srcdir)/test_utils/LinSysCore.hpp \
	$(top_srcdir)/test_utils/PoissonData.hpp \
	$(top_srcdir)/test_utils/Poisson_Elem.hpp \
	$(top_srcdir)/test_utils/poisson_beam_mains.hpp \
	$(top_srcdir)/test_utils/snl_fei_tester.hpp \
	$(top_srcdir)/test_utils/SolnCheck.hpp \
	$(top_srcdir)/test_utils/test_Algebraic.hpp \
	$(top_srcdir)/test_utils/test_AztecWrappers.hpp \
	$(top_srcdir)/test_utils/test_Database.hpp \
	$(top_srcdir)/test_utils/testData.hpp \
	$(top_srcdir)/test_utils/test_EqnBuffer.hpp \
	$(top_srcdir)/test_utils/test_EqnCommMgr.hpp \
	$(top_srcdir)/test_utils/tester.hpp \
	$(top_srcdir)/test_utils/test_Factory.hpp \
	$(top_srcdir)/test_utils/test_Factory_helper.hpp \
	$(top_srcdir)/test_utils/test_FEI.hpp \
	$(top_srcdir)/test_utils/test_FEI_Implementation.hpp \
	$(top_srcdir)/test_utils/test_LinearSystem.hpp \
	$(top_srcdir)/test_utils/test_MatrixGraph.hpp \
	$(top_srcdir)/test_utils/test_Matrix.hpp \
	$(top_srcdir)/test_utils/test_misc.hpp \
	$(top_srcdir)/test_utils/test_PointBlockMap.hpp \
	$(top_srcdir)/test_utils/test_Set.hpp \
	$(top_srcdir)/test_utils/test_SNL_FEI_Structure.hpp \
	$(top_srcdir)/test_utils/test_FEI_Impl.hpp \
	$(top_srcdir)/test_utils/test_Tables.hpp \
	$(top_srcdir)/test_utils/test_Utils.hpp \
	$(top_srcdir)/test_utils/test_Vector.hpp \
	$(top_srcdir)/test_utils/test_VectorSpace.hpp \
	$(top_srcdir)/test_utils/test_benchmarks.hpp

CORE = \
	$(srcdir)/BCNodeSet.cpp \
	$(srcdir)/CommNodeSet.cpp \
	$(srcdir)/CRSet.cpp \
	$(srcdir)/cFeiTester_main.cpp \
	$(srcdir)/beam_main.cpp \
	$(srcdir)/beam_oldfei_main.cpp \
	$(srcdir)/DataReader.cpp \
	$(srcdir)/driverData.cpp \
	$(srcdir)/ElemBlock.cpp \
	$(srcdir)/feiDriver_main.cpp \
	$(srcdir)/fei_test_utils.cpp \
	$(srcdir)/FEData.cpp \
	$(srcdir)/FEI_tester.cpp \
	$(srcdir)/HexBeam.cpp \
	$(srcdir)/HexBeamCR.cpp \
	$(srcdir)/InputFileReader.cpp \
	$(srcdir)/LibraryFactory.cpp \
	$(srcdir)/PoissonData.cpp \
	$(srcdir)/Poisson_Elem.cpp \
	$(srcdir)/poisson3_main.cpp \
	$(srcdir)/poisson_main.cpp \
	$(srcdir)/snl_fei_tester.cpp \
	$(srcdir)/SolnCheck.cpp \
	$(srcdir)/tester.cpp \
	$(srcdir)/test_Algebraic.cpp \
	$(srcdir)/test_AztecWrappers.cpp \
	$(srcdir)/test_Database.cpp \
	$(srcdir)/test_EqnBuffer.cpp \
	$(srcdir)/test_EqnCommMgr.cpp \
	$(srcdir)/test_Factory.cpp \
	$(srcdir)/test_Factory_helper.cpp \
	$(srcdir)/test_FEI.cpp \
	$(srcdir)/test_FEI_Implementation.cpp \
	$(srcdir)/test_LinearSystem.cpp \
	$(srcdir)/test_Matrix.cpp \
	$(srcdir)/test_MatrixGraph.cpp \
	$(srcdir)/test_misc.cpp \
	$(srcdir)/test_PointBlockMap.cpp \
	$(srcdir)/test_Set.cpp \
	$(srcdir)/test_SNL_FEI_Structure.cpp \
	$(srcdir)/test_FEI_Impl.cpp \
	$(srcdir)/test_Tables.cpp \
	$(srcdir)/test_Utils.cpp \
	$(srcdir)/test_Vector.cpp \
	$(srcdir)/test_VectorSpace.cpp \
	$(srcdir)/test_benchmarks.cpp


# ------------------------------------------------------------------------
# FEI library specifications
# ------------------------------------------------------------------------

FEI_LIB = libfei_test_utils.a

noinst_LIBRARIES = $(FEI_LIB)


libfei_test_utils_a_SOURCES = \
	$(CORE_H) $(CORE) 

# EXTRA_libfei_test_utils_a_SOURCES = 

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR if USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libfei_test_utils_a_AR = $(ALTERNATE_AR)
else

libfei_test_utils_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libfei_test_utils_a_LIBADD = $(FEI_TEMPLATE_OBJS) $(XTRALDADD)

includedir = ${prefix}/test_utils

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

