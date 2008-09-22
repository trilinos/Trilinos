# @HEADER
# ************************************************************************
# 
#        Phalanx: A Partial Differential Equation Field Evaluation 
#       Kernel for Flexible Management of Complex Dependency Chains
#                  Copyright (2008) Sandia Corporation
# 
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
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
# 
# Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
# National Laboratories.
# 
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package phalanx subdirectory src built successfully."
	@echo ""


# ------------------------------------------------------------------------
# List of all files
# ------------------------------------------------------------------------

PHALANX_HPP = \
	$(srcdir)/Phalanx.hpp \
	$(srcdir)/Phalanx_Allocator_Contiguous.hpp \
	$(srcdir)/Phalanx_Allocator_New.hpp \
	$(srcdir)/Phalanx_ArrayPrivate.hpp \
	$(srcdir)/Phalanx_Array.hpp \
	$(srcdir)/Phalanx_ConfigDefs.hpp \
	$(srcdir)/Phalanx_DataContainer_Base_Def.hpp \
	$(srcdir)/Phalanx_DataContainer_Base.hpp \
	$(srcdir)/Phalanx_DataContainer_Def.hpp \
	$(srcdir)/Phalanx_DataContainer.hpp \
	$(srcdir)/Phalanx_DataContainer_TemplateManager.hpp \
	$(srcdir)/Phalanx_DataLayout_FlatLayout.hpp \
	$(srcdir)/Phalanx_DataLayout_MDALayout_Def.hpp \
	$(srcdir)/Phalanx_DataLayout_MDALayout.hpp \
	$(srcdir)/Phalanx_DataLayout.hpp \
	$(srcdir)/Phalanx_DoxygenDocumentation.hpp \
	$(srcdir)/Phalanx_EvaluationContainer_Base_Def.hpp \
	$(srcdir)/Phalanx_EvaluationContainer_Base.hpp \
	$(srcdir)/Phalanx_EvaluationContainer_Def.hpp \
	$(srcdir)/Phalanx_EvaluationContainer.hpp \
	$(srcdir)/Phalanx_EvaluationContainer_TemplateBuilder.hpp \
	$(srcdir)/Phalanx_EvaluationContainer_TemplateManager.hpp \
	$(srcdir)/Phalanx_Evaluator_Base.hpp \
	$(srcdir)/Phalanx_Evaluator_Derived.hpp \
	$(srcdir)/Phalanx_Evaluator_Factory_Def.hpp \
	$(srcdir)/Phalanx_Evaluator_Factory.hpp \
	$(srcdir)/Phalanx_Evaluator_Factory_UFO.hpp \
	$(srcdir)/Phalanx_Evaluator.hpp \
	$(srcdir)/Phalanx_Evaluator_Macros.hpp \
	$(srcdir)/Phalanx_Evaluator_Manager_Def.hpp \
	$(srcdir)/Phalanx_Evaluator_Manager.hpp \
	$(srcdir)/Phalanx_Evaluator_TemplateBuilder.hpp \
	$(srcdir)/Phalanx_Evaluator_TemplateManager.hpp \
	$(srcdir)/Phalanx_Evaluator_Utilities.hpp \
	$(srcdir)/Phalanx_Evaluator_WithBaseImpl_Def.hpp \
	$(srcdir)/Phalanx_Evaluator_WithBaseImpl.hpp \
	$(srcdir)/Phalanx_Field_Def.hpp \
	$(srcdir)/Phalanx_Field.hpp \
	$(srcdir)/Phalanx_FieldManager_Def.hpp \
	$(srcdir)/Phalanx_FieldManager.hpp \
	$(srcdir)/Phalanx_FieldTag.hpp \
	$(srcdir)/Phalanx_FieldTag_STL_Functors.hpp \
	$(srcdir)/Phalanx_FieldTag_Tag_Def.hpp \
	$(srcdir)/Phalanx_FieldTag_Tag.hpp \
	$(srcdir)/Phalanx_MDField_Def.hpp \
	$(srcdir)/Phalanx_MDField.hpp \
	$(srcdir)/Phalanx_Print_Utilities.hpp \
	$(srcdir)/Phalanx_SimpleArrayOps.hpp \
	$(srcdir)/Phalanx_TemplateIterator.hpp \
	$(srcdir)/Phalanx_TemplateManager_Def.hpp \
	$(srcdir)/Phalanx_TemplateManager.hpp \
	$(srcdir)/Phalanx_Traits_Base.hpp \
	$(srcdir)/Phalanx_TypeStrings.hpp \
	$(srcdir)/Phalanx_Version.hpp

PHALANX_CPP = \
	$(srcdir)/Phalanx_ArrayPrivate.cpp \
	$(srcdir)/Phalanx_DataLayout_FlatLayout.cpp \
	$(srcdir)/Phalanx_DataLayout.cpp \
	$(srcdir)/Phalanx_FieldTag.cpp

# ------------------------------------------------------------------------
# Includes for required and optional library dependencies
# ------------------------------------------------------------------------
include $(top_builddir)/Makefile.export.phalanx

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(PHALANX_INCLUDES))
else
EXPORT_INCLUDES = $(PHALANX_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

#np# The following section can be duplicated for building multiple libraries -
#np# make sure to list the correct categories with each library.
#np# replace new_package with the name of the package being autotool'ed here
# ------------------------------------------------------------------------
# Phalanx library specifications
# ------------------------------------------------------------------------
#np# replace new_package with the name of the package being autotool'ed here
PHALANX_LIB = libphalanx.a

#np# replace new_package with the name of the package being autotool'ed here
libphalanx_a_SOURCES = $(PHALANX_CPP)

#np# replace new_package with the name of the package being autotool'ed here
EXTRA_libphalanx_a_SOURCES =

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

#np# replace new_package with the name of the package being autotool'ed here
libphalanx_a_AR = $(ALTERNATE_AR)
else

#np# replace new_package with the name of the package being autotool'ed here
libphalanx_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
libphalanx_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
lib_LIBRARIES = libphalanx.a

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
include_HEADERS = $(PHALANX_HPP) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
# SPECIAL NOTE: New_Package_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
#np# replace new_package with the name of the package being autotool'ed here
nodist_execinclude_HEADERS = Phalanx_config.hpp

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in
