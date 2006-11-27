## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#                           Sacado Package
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
# Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
# (etphipp@sandia.gov).
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside Sacado. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

SACADO_FAD_H = \
	$(srcdir)/Sacado_Fad_GeneralFad.hpp \
	$(srcdir)/Sacado_Fad_GeneralFadImp.hpp \
	$(srcdir)/Sacado_Fad_GeneralFadExpr.hpp \
	$(srcdir)/Sacado_Fad_ExpressionTraits.hpp \
	$(srcdir)/Sacado_Fad_Expression.hpp \
	$(srcdir)/Sacado_Fad_Ops.hpp \
	$(srcdir)/Sacado_Fad_DynamicStorage.hpp \
	$(srcdir)/Sacado_Fad_DFad.hpp \
	$(srcdir)/Sacado_Fad_DFadTraits.hpp \
	$(srcdir)/Sacado_Fad_StaticStorage.hpp \
	$(srcdir)/Sacado_Fad_SLFad.hpp \
	$(srcdir)/Sacado_Fad_SLFadTraits.hpp \
	$(srcdir)/Sacado_Fad_SFad.hpp \
	$(srcdir)/Sacado_Fad_SFadImp.hpp \
	$(srcdir)/Sacado_Fad_SFadTraits.hpp \
	$(srcdir)/Sacado_Fad_MemPool.hpp \
	$(srcdir)/Sacado_Fad_MemPoolImp.hpp \
	$(srcdir)/Sacado_Fad_MemPoolManager.hpp \
	$(srcdir)/Sacado_Fad_MemPoolManagerImp.hpp \
	$(srcdir)/Sacado_Fad_MemPoolStorage.hpp \
	$(srcdir)/Sacado_Fad_DMFad.hpp \
	$(srcdir)/Sacado_Fad_DMFadTraits.hpp

SACADO_CACHEFAD_H = \
	$(srcdir)/Sacado_CacheFad_GeneralFad.hpp \
	$(srcdir)/Sacado_CacheFad_GeneralFadImp.hpp \
	$(srcdir)/Sacado_CacheFad_GeneralFadExpr.hpp \
	$(srcdir)/Sacado_CacheFad_Expression.hpp \
	$(srcdir)/Sacado_CacheFad_Ops.hpp \
	$(srcdir)/Sacado_CacheFad_DFadTraits.hpp \
	$(srcdir)/Sacado_CacheFad_DFad.hpp 

SACADO_TAYLOR_H = \
	$(srcdir)/Sacado_Taylor_DTaylorTraits.hpp \
	$(srcdir)/Sacado_Taylor_DTaylor.hpp \
	$(srcdir)/Sacado_Taylor_DTaylorImp.hpp \
	$(srcdir)/Sacado_Taylor_Expression.hpp \
	$(srcdir)/Sacado_Taylor_Ops.hpp 

SACADO_RAD_H = \
	$(srcdir)/Sacado_trad_Traits.hpp \
	$(srcdir)/Sacado_trad.hpp \
	$(srcdir)/Sacado_rad.hpp

SACADO_RAD = \
	$(srcdir)/Sacado_radops.cpp \
	$(srcdir)/uninit.c

SACADO_TEMPLATE_H = \
	$(srcdir)/Sacado_TemplateManager.hpp \
	$(srcdir)/Sacado_TemplateManagerImp.hpp \
	$(srcdir)/Sacado_TemplateIterator.hpp 

SACADO_MPL_H = \
	$(srcdir)/mpl/Sacado_mpl_is_same.hpp \
	$(srcdir)/mpl/Sacado_mpl_if.hpp \
	$(srcdir)/mpl/Sacado_mpl_none.hpp \
	$(srcdir)/mpl/Sacado_mpl_begin.hpp \
	$(srcdir)/mpl/Sacado_mpl_end.hpp \
	$(srcdir)/mpl/Sacado_mpl_next.hpp \
	$(srcdir)/mpl/Sacado_mpl_size.hpp \
	$(srcdir)/mpl/Sacado_mpl_at.hpp \
	$(srcdir)/mpl/Sacado_mpl_deref.hpp \
	$(srcdir)/mpl/Sacado_mpl_push_back.hpp \
	$(srcdir)/mpl/Sacado_mpl_find.hpp \
	$(srcdir)/mpl/Sacado_mpl_for_each.hpp \
	$(srcdir)/mpl/Sacado_mpl_vector.hpp \
	$(srcdir)/mpl/Sacado_mpl_vector_size_spec.hpp \
	$(srcdir)/mpl/Sacado_mpl_vector_at_spec.hpp \
	$(srcdir)/mpl/Sacado_mpl_vector_push_back_spec.hpp

SACADO_H = \
	$(srcdir)/Sacado_ConfigDefs.h \
	$(srcdir)/Sacado_Version.hpp \
	$(srcdir)/Sacado_Traits.hpp \
	$(srcdir)/Sacado.hpp \
	$(SACADO_FAD_H) \
	$(SACADO_CACHEFAD_H) \
	$(SACADO_TAYLOR_H) \
	$(SACADO_RAD_H) \
	$(SACADO_MPL_H) \
	$(SACADO_TEMPLATE_H)

SACADO = \
	$(SACADO_RAD)

# ------------------------------------------------------------------------
# Sacado library specifications
# ------------------------------------------------------------------------
SACADO_LIB = libsacado.a

libsacado_a_SOURCES = $(SACADO) $(SACADO_H)

EXTRA_libsacado_a_SOURCES = 

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libsacado_a_AR = $(ALTERNATE_AR)

else

libsacado_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libsacado_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(SACADO_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(SACADO_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Sacado_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory.  This directory can
# be set via --with-execincludedir=* configure option.

nodist_execinclude_HEADERS = Sacado_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 





