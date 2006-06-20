# @HEADER
# ************************************************************************
# 
#                    Teuchos: Common Tools Package
#                 Copyright (2004) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package teuchos subdirectory src built successfully."
	@echo ""

# Note the following paragraph.
# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
# There is, at present, no distinction between most of the categories, 
# other than those which start with EXTRA_.  All categories which do not
# start with EXTRA_ are included in all builds.  Categories which start
# with EXTRA_ are included in some builds (typically controlled by 
# --with-xyz or --without-xyz)
#
# Make sure to list all source files in one of the following categories.
CORE_H = \
	$(srcdir)/Teuchos_AbstractFactory.hpp \
	$(srcdir)/Teuchos_AbstractFactoryStd.hpp \
	$(srcdir)/Teuchos_any.hpp \
	$(srcdir)/Teuchos_Array.hpp \
	$(srcdir)/Teuchos_arrayArg.hpp \
	$(srcdir)/Teuchos_basic_oblackholestream.hpp \
	$(srcdir)/Teuchos_BLAS.hpp \
	$(srcdir)/Teuchos_BLAS_wrappers.hpp \
	$(srcdir)/Teuchos_BLAS_types.hpp \
	$(srcdir)/Teuchos_Comm.hpp \
  $(srcdir)/Teuchos_CommHelpers.hpp \
  $(srcdir)/Teuchos_CommUtilities.hpp \
	$(srcdir)/Teuchos_CommandLineProcessor.hpp \
	$(srcdir)/Teuchos_CompileTimeAssert.hpp \
	$(srcdir)/Teuchos_CompObject.hpp \
	$(srcdir)/Teuchos_ConfigDefs.hpp \
  $(srcdir)/Teuchos_ConstNonconstObjectContainer.hpp \
  $(srcdir)/Teuchos_DefaultComm.hpp \
  $(srcdir)/Teuchos_DefaultMpiComm.hpp \
	$(srcdir)/Teuchos_DefaultSerialComm.hpp \
	$(srcdir)/Teuchos_DataAccess.hpp \
	$(srcdir)/Teuchos_Describable.hpp \
	$(srcdir)/Teuchos_dyn_cast.hpp \
	$(srcdir)/Teuchos_ErrorPolling.hpp \
	$(srcdir)/Teuchos_exit.h \
	$(srcdir)/Teuchos_FancyOStream.hpp \
	$(srcdir)/Teuchos_Flops.hpp \
	$(srcdir)/Teuchos_F77_wrappers.h \
	$(srcdir)/Teuchos_GlobalMPISession.hpp \
	$(srcdir)/Teuchos_getConst.hpp \
	$(srcdir)/Teuchos_LAPACK.hpp \
	$(srcdir)/Teuchos_LAPACK_wrappers.hpp \
	$(srcdir)/Teuchos_map.hpp \
	$(srcdir)/Teuchos_oblackholestream.hpp \
	$(srcdir)/Teuchos_MPIComm.hpp \
	$(srcdir)/Teuchos_MPIContainerComm.hpp \
	$(srcdir)/Teuchos_MPISession.hpp \
	$(srcdir)/Teuchos_MPITraits.hpp \
	$(srcdir)/Teuchos_Object.hpp \
  $(srcdir)/Teuchos_OpaqueWrapper.hpp \
	$(srcdir)/Teuchos_OrdinalTraits.hpp \
	$(srcdir)/Teuchos_ParameterEntry.hpp \
	$(srcdir)/Teuchos_ParameterList.hpp \
	$(srcdir)/Teuchos_ParameterListAcceptor.hpp \
	$(srcdir)/Teuchos_PrimitiveTypeTraits.hpp \
  $(srcdir)/Teuchos_ReductionOp.hpp \
  $(srcdir)/Teuchos_ReductionOpHelpers.hpp \
	$(srcdir)/Teuchos_PerformanceMonitorBase.hpp \
	$(srcdir)/Teuchos_PerformanceMonitorUtils.hpp \
	$(srcdir)/Teuchos_PolynomialTraits.hpp \
	$(srcdir)/Teuchos_PolynomialDecl.hpp \
	$(srcdir)/Teuchos_Polynomial.hpp \
	$(srcdir)/Teuchos_PromotionTraits.hpp \
	$(srcdir)/Teuchos_Range1D.hpp \
	$(srcdir)/Teuchos_RawMPITraits.hpp \
	$(srcdir)/Teuchos_ReductionOp.hpp \
	$(srcdir)/Teuchos_RefCountPtrDecl.hpp \
	$(srcdir)/Teuchos_RefCountPtr.hpp \
	$(srcdir)/Teuchos_RefCountPtrBoostSharedPtrConversionsDecl.hpp \
	$(srcdir)/Teuchos_RefCountPtrBoostSharedPtrConversions.hpp \
	$(srcdir)/Teuchos_ScalarTraits.hpp \
  $(srcdir)/Teuchos_SerializationTraits.hpp \
  $(srcdir)/Teuchos_SerializationTraitsHelpers.hpp \
  $(srcdir)/Teuchos_Serializer.hpp \
  $(srcdir)/Teuchos_SerializerHelpers.hpp \
	$(srcdir)/Teuchos_SerialDenseMatrix.hpp \
	$(srcdir)/Teuchos_SerialDenseVector.hpp \
	$(srcdir)/Teuchos_StandardCatchMacros.hpp \
	$(srcdir)/Teuchos_StandardCompositionMacros.hpp \
	$(srcdir)/Teuchos_StandardMemberCompositionMacros.hpp \
	$(srcdir)/Teuchos_StringToIntMap.hpp \
	$(srcdir)/Teuchos_TableColumn.hpp \
	$(srcdir)/Teuchos_TableEntry.hpp \
	$(srcdir)/Teuchos_TableFormat.hpp \
	$(srcdir)/Teuchos_TestForException.hpp \
	$(srcdir)/Teuchos_Time.hpp \
	$(srcdir)/Teuchos_TimeMonitor.hpp \
	$(srcdir)/Teuchos_Utils.hpp \
	$(srcdir)/Teuchos_VerboseObject.hpp \
	$(srcdir)/Teuchos_VerbosityLevel.hpp \
	$(srcdir)/Teuchos_Workspace.hpp

CORE = \
	$(srcdir)/Teuchos_BLAS.cpp \
	$(srcdir)/Teuchos_CommandLineProcessor.cpp \
	$(srcdir)/Teuchos_CompObject.cpp \
	$(srcdir)/Teuchos_Describable.cpp \
	$(srcdir)/Teuchos_dyn_cast.cpp \
	$(srcdir)/Teuchos_ErrorPolling.cpp \
	$(srcdir)/Teuchos_exit.cpp \
	$(srcdir)/Teuchos_Flops.cpp \
	$(srcdir)/Teuchos_GlobalMPISession.cpp \
	$(srcdir)/Teuchos_MPIComm.cpp \
	$(srcdir)/Teuchos_MPISession.cpp \
	$(srcdir)/Teuchos_Object.cpp \
	$(srcdir)/Teuchos_ParameterEntry.cpp \
	$(srcdir)/Teuchos_ParameterList.cpp \
	$(srcdir)/Teuchos_ParameterListAcceptor.cpp \
	$(srcdir)/Teuchos_PerformanceMonitorUtils.cpp \
	$(srcdir)/Teuchos_Range1D.cpp \
	$(srcdir)/Teuchos_RefCountPtr.cpp \
	$(srcdir)/Teuchos_ScalarTraits.cpp \
	$(srcdir)/Teuchos_StringToIntMap.cpp \
	$(srcdir)/Teuchos_TableColumn.cpp \
	$(srcdir)/Teuchos_TableEntry.cpp \
	$(srcdir)/Teuchos_TableFormat.cpp \
	$(srcdir)/Teuchos_TestForException.cpp \
	$(srcdir)/Teuchos_Time.cpp \
	$(srcdir)/Teuchos_TimeMonitor.cpp \
	$(srcdir)/Teuchos_Utils.cpp \
	$(srcdir)/Teuchos_VerboseObject.cpp \
	$(srcdir)/Teuchos_Workspace.cpp

EXTRA_EXTENDED_H = \
	$(srcdir)/Teuchos_ExpatHandlerAdapter.hpp \
	$(srcdir)/Teuchos_FileInputSource.hpp \
	$(srcdir)/Teuchos_FileInputStream.hpp \
	$(srcdir)/Teuchos_HashSet.hpp \
	$(srcdir)/Teuchos_Hashtable.hpp \
	$(srcdir)/Teuchos_HashUtils.hpp \
	$(srcdir)/Teuchos_Libxml2HandlerAdapter.hpp \
	$(srcdir)/Teuchos_ParameterXMLFileReader.hpp \
	$(srcdir)/Teuchos_StringInputSource.hpp \
	$(srcdir)/Teuchos_StringInputStream.hpp \
	$(srcdir)/Teuchos_StrUtils.hpp \
	$(srcdir)/Teuchos_TreeBuildingXMLHandler.hpp\
	$(srcdir)/Teuchos_XMLInputSource.hpp\
	$(srcdir)/Teuchos_XMLInputStream.hpp\
	$(srcdir)/Teuchos_XMLObject.hpp\
	$(srcdir)/Teuchos_XMLObjectImplem.hpp \
	$(srcdir)/Teuchos_XMLParameterListHelpers.hpp \
	$(srcdir)/Teuchos_XMLParameterListReader.hpp \
	$(srcdir)/Teuchos_XMLParameterListWriter.hpp \
	$(srcdir)/Teuchos_XMLParser.hpp

EXTRA_EXTENDED = \
	$(srcdir)/Teuchos_ExpatHandlerAdapter.cpp \
	$(srcdir)/Teuchos_FileInputSource.cpp \
	$(srcdir)/Teuchos_FileInputStream.cpp \
	$(srcdir)/Teuchos_HashUtils.cpp \
	$(srcdir)/Teuchos_Libxml2HandlerAdapter.cpp \
	$(srcdir)/Teuchos_ParameterXMLFileReader.cpp \
	$(srcdir)/Teuchos_StrUtils.cpp \
	$(srcdir)/Teuchos_StringInputSource.cpp \
	$(srcdir)/Teuchos_StringInputStream.cpp \
	$(srcdir)/Teuchos_TreeBuildingXMLHandler.cpp \
	$(srcdir)/Teuchos_XMLInputSource.cpp \
	$(srcdir)/Teuchos_XMLInputStream.cpp \
	$(srcdir)/Teuchos_XMLObject.cpp \
	$(srcdir)/Teuchos_XMLObjectImplem.cpp \
	$(srcdir)/Teuchos_XMLParameterListHelpers.cpp \
	$(srcdir)/Teuchos_XMLParameterListReader.cpp \
	$(srcdir)/Teuchos_XMLParameterListWriter.cpp \
	$(srcdir)/Teuchos_XMLParser.cpp 

UTIL_H = \
	$(srcdir)/Teuchos_Version.hpp

UTIL =

ENUM_H = 

UNIPROC_H = 

UNIPROC = 

SERIAL_H = 

SERIAL = 

#
#  Files which are included in select builds are typically have
#  names that start with EXTRA
#
#  EXTRA_MPI is control by --with-mpi and --without-mpi which is
#  set up by TAC_ARC_CONFIG_MPI 
#
EXTRA_MPI_H = \
	$(srcdir)/Teuchos_MpiComm.hpp \
	$(srcdir)/Teuchos_MpiReductionOpSetter.hpp 

EXTRA_MPI = \
	$(srcdir)/Teuchos_MpiReductionOpSetter.cpp 

if HAVE_MPI 
MPI_H = \
	$(EXTRA_MPI_H)
MPI = \
	$(EXTRA_MPI)

else

MPI_H =
MPI =

endif

# -------------------------------------------------------------
# Enable compilation of experimental code if feature is enabled
# -------------------------------------------------------------
if HAVE_TEUCHOS_EXTENDED

EXTENDED_H = $(EXTRA_EXTENDED_H)
EXTENDED = $(EXTRA_EXTENDED)

else

EXTENDED_H =
EXTENDED =

endif

USERBASE_H = 

LANG_WRAPS_H = 

LANG_WRAPS = 

DIST_H = 

DIST = 

COMP_H =

COMP = 

# The following section can be duplicated for building multiple libraries -
# make sure to list the correct categories with each library.
# replace teuchos with the name of the package being autotool'ed here
# ------------------------------------------------------------------------
# Teuchos library specifications
# ------------------------------------------------------------------------
# replace teuchos with the name of the package being autotool'ed here
TEUCHOS_LIB = libteuchos.a

# replace new_package with the name of the package being autotool'ed here
TEUCHOS_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(ENUM_H) \
	$(UNIPROC_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) \
	$(EXTENDED_H)


# replace new_package with the name of the package being autotool'ed here
libteuchos_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(UNIPROC) \
	$(SERIAL) \
	$(DIST) \
	$(COMP) \
	$(MPI) \
	$(EXTENDED)

# replace new_package with the name of the package being autotool'ed here
EXTRA_libteuchos_a_SOURCES = \
	$(EXTRA_MPI) \
	$(EXTRA_MPI_H) \
	$(EXTRA_ARPREC) \
	$(EXTRA_ARPREC_H) 

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libteuchos_a_AR = $(ALTERNATE_AR)
else

libteuchos_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libteuchos_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD) 

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(TEUCHOS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(TEUCHOS_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: New_Package_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = Teuchos_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Library Target for Circular Dependencies
# ------------------------------------------------------------------------

all-libs:
	$(MAKE) libteuchos.a
