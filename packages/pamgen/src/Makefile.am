# @HEADER
# ************************************************************************
# 
#                         Pamgen Package
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
	@echo "Trilinos package pamgen subdirectory pamgen_src built successfully."
	@echo ""

#np# Note the following paragraph.
# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside Pamgen. The first is
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
CORE = \
	$(srcdir)/Random.h \
	$(srcdir)/StrLoopLimits.h \
	$(srcdir)/Vector.h \
	$(srcdir)/bc_specification.h \
	$(srcdir)/brick_inline_mesh_desc.C \
	$(srcdir)/brick_inline_mesh_desc.h \
	$(srcdir)/calc_decomp_cuts.C \
	$(srcdir)/calc_decomp_cuts.h \
	$(srcdir)/create_inline_mesh.C \
	$(srcdir)/element_density_function.C \
	$(srcdir)/element_density_function.h \
	$(srcdir)/geometry_transform.C \
	$(srcdir)/geometry_transform.h \
	$(srcdir)/inline_geometries.h \
	$(srcdir)/inline_mes_driver_lt.h \
	$(srcdir)/inline_mesh_desc.C \
	$(srcdir)/inline_mesh_desc.h \
	$(srcdir)/inline_mesh_driver_lt.C \
	$(srcdir)/inline_mesh_driver_lt.h \
	$(srcdir)/legacy_inline_mesh_desc.C \
	$(srcdir)/legacy_inline_mesh_desc.h \
	$(srcdir)/parse_routines.C \
	$(srcdir)/parse_routines.h \
	$(srcdir)/radial_inline_mesh_desc.C \
	$(srcdir)/radial_inline_mesh_desc.h \
	$(srcdir)/radial_trisection_inline_mesh_desc.C \
	$(srcdir)/radial_trisection_inline_mesh_desc.h \
	$(srcdir)/topology_enum.h \
	$(srcdir)/uns_inline_decomp.C \
	$(srcdir)/uns_inline_decomp.h

PARSER = $(srcdir)/../parser/InputBlock.C \
	$(srcdir)/../parser/InputBlock.h \
	$(srcdir)/../parser/keyword.h \
	$(srcdir)/../parser/lexer.C \
	$(srcdir)/../parser/parse_table.C \
	$(srcdir)/../parser/parse_table.h \
	$(srcdir)/../parser/parser.C \
	$(srcdir)/../parser/parser.h \
	$(srcdir)/../parser/token.C \
	$(srcdir)/../parser/token.h \
	$(srcdir)/../parser/token_enum.h \
	$(srcdir)/../parser/token_stream.C \
	$(srcdir)/../parser/token_stream.h \
	$(srcdir)/../parser/token_value.h 

MS_LT = $(srcdir)/../mesh_spec_lt/im_ex_c_interface.C \
	$(srcdir)/../mesh_spec_lt/im_ex_c_interface_l.C \
	$(srcdir)/../mesh_spec_lt/im_exodusII.h \
	$(srcdir)/../mesh_spec_lt/im_exodusII_ext.h \
	$(srcdir)/../mesh_spec_lt/im_ne_c_interface.C \
	$(srcdir)/../mesh_spec_lt/im_ne_c_interface_l.C \
	$(srcdir)/../mesh_spec_lt/im_ne_nemesisI.h \
	$(srcdir)/../mesh_spec_lt/pamgen_mesh_specification.C \
	$(srcdir)/../mesh_spec_lt/pamgen_mesh_specification.h

RTCOMPILER = $(srcdir)/../rtcompiler/ArrayIndexRTC.hh \
	$(srcdir)/../rtcompiler/ArrayNumberRTC.hh \
	$(srcdir)/../rtcompiler/ArrayVarRTC.hh \
	$(srcdir)/../rtcompiler/Bessel_I.C \
	$(srcdir)/../rtcompiler/Bessel_I.hh \
	$(srcdir)/../rtcompiler/Block.C \
	$(srcdir)/../rtcompiler/BlockRTC.hh \
	$(srcdir)/../rtcompiler/ConditionalBlock.C \
	$(srcdir)/../rtcompiler/ConditionalBlockRTC.hh  \
	$(srcdir)/../rtcompiler/ExecutableRTC.hh \
	$(srcdir)/../rtcompiler/ForBlock.C \
	$(srcdir)/../rtcompiler/ForBlockRTC.hh \
	$(srcdir)/../rtcompiler/Function.C \
	$(srcdir)/../rtcompiler/FunctionRTC.hh \
	$(srcdir)/../rtcompiler/IfElseifElseBlock.C \
	$(srcdir)/../rtcompiler/IfElseifElseBlockRTC.hh \
	$(srcdir)/../rtcompiler/Line.C \
	$(srcdir)/../rtcompiler/LineRTC.hh \
	$(srcdir)/../rtcompiler/NormalBlock.C \
	$(srcdir)/../rtcompiler/NormalBlockRTC.hh \
	$(srcdir)/../rtcompiler/ObjectRTC.hh \
	$(srcdir)/../rtcompiler/Operator.C \
	$(srcdir)/../rtcompiler/OperatorRTC.hh \
	$(srcdir)/../rtcompiler/Registrar.C \
	$(srcdir)/../rtcompiler/RegistrarRTC.hh \
	$(srcdir)/../rtcompiler/ScalarNumberRTC.hh \
	$(srcdir)/../rtcompiler/ScalarVarRTC.hh \
	$(srcdir)/../rtcompiler/Tokenizer.C \
	$(srcdir)/../rtcompiler/TokenizerRTC.hh \
	$(srcdir)/../rtcompiler/ValueRTC.hh \
	$(srcdir)/../rtcompiler/VariableRTC.hh \
	$(srcdir)/../rtcompiler/WhileBlock.C \
	$(srcdir)/../rtcompiler/WhileBlockRTC.hh \
	$(srcdir)/../rtcompiler/common.C \
	$(srcdir)/../rtcompiler/commonRTC.hh

ASRC = $(srcdir)/../asrc/code_types.h \
	$(srcdir)/../asrc/element_dictionary.C \
	$(srcdir)/../asrc/element_dictionary.h \
	$(srcdir)/../asrc/fudges.C \
	$(srcdir)/../asrc/fudges.h 


CORE_H =  $(srcdir)/create_inline_mesh.h \
	$(srcdir)/../mesh_spec_lt/im_exodusII.h \
	$(srcdir)/../mesh_spec_lt/im_exodusII_ext.h \
	$(srcdir)/../mesh_spec_lt/im_ne_nemesisI.h 

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
AM_CPPFLAGS = -I../rtcompiler -I../mesh_spec_lt

#np# The following section can be duplicated for building multiple libraries -
#np# make sure to list the correct categories with each library.
#np# replace new_package with the name of the package being autotool'ed here
# ------------------------------------------------------------------------
# Pamgen library specifications
# ------------------------------------------------------------------------
#np# replace new_package with the name of the package being autotool'ed here
PAMGEN_LIB = libpamgen.a

#np# replace new_package with the name of the package being autotool'ed here
PAMGEN_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(ENUM_H) \
	$(UNIPROC_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) \
	$(NEWP_SWAHILI_H)

#np# replace new_package with the name of the package being autotool'ed here
libpamgen_a_SOURCES = \
	$(CORE) \
	$(PARSER) \
	$(ASRC) \
	$(RTCOMPILER) \
	$(MS_LT) \
	$(UTIL) \
	$(UNIPROC) \
	$(SERIAL) \
	$(DIST) \
	$(COMP) \
	$(MPI) \
	$(NEWP_SWAHILI)

#np# replace new_package with the name of the package being autotool'ed here
EXTRA_libpamgen_a_SOURCES = \
	$(EXTRA_MPI) \
	$(EXTRA_MPI_H) \
	$(EXTRA_NEWP_SWAHILI) \
	$(EXTRA_NEWP_SWAHILI_H)

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
libpamgen_a_AR = $(ALTERNATE_AR)
else

#np# replace new_package with the name of the package being autotool'ed here
libpamgen_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
libpamgen_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
lib_LIBRARIES = $(PAMGEN_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
include_HEADERS = $(PAMGEN_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
# SPECIAL NOTE: Pamgen_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(includedir)
#np# replace new_package with the name of the package being autotool'ed here
nodist_execinclude_HEADERS = Pamgen_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






