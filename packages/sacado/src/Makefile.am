## $Id$
## $Source$
# @HEADER
# ************************************************************************
#
#                           Sacado Package
#                 Copyright (2006) Sandia Corporation
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
# Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
# (etphipp@sandia.gov).
#
# ************************************************************************
# @HEADER

all-local:
	@echo ""
	@echo "Trilinos package sacado subdirectory src built successfully."
	@echo ""

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
	$(srcdir)/Sacado_Fad_DMFadTraits.hpp \
	$(srcdir)/Sacado_Fad_DVFad.hpp \
	$(srcdir)/Sacado_Fad_DVFadTraits.hpp \
	$(srcdir)/Sacado_Fad_Vector.hpp \
	$(srcdir)/Sacado_Fad_VectorImp.hpp \
	$(srcdir)/Sacado_Fad_VectorDynamicStorage.hpp \
	$(srcdir)/Sacado_Fad_SimpleFad.hpp \
	$(srcdir)/Sacado_Fad_SimpleFadImp.hpp \
	$(srcdir)/Sacado_Fad_SimpleFadOps.hpp \
	$(srcdir)/Sacado_Fad_SimpleFadTraits.hpp

SACADO_ELRFAD_H = \
	$(srcdir)/Sacado_ELRFad_GeneralFad.hpp \
	$(srcdir)/Sacado_ELRFad_GeneralFadImp.hpp \
	$(srcdir)/Sacado_ELRFad_GeneralFadExpr.hpp \
	$(srcdir)/Sacado_ELRFad_ExpressionTraits.hpp \
	$(srcdir)/Sacado_ELRFad_Expression.hpp \
	$(srcdir)/Sacado_ELRFad_Ops.hpp \
	$(srcdir)/Sacado_ELRFad_DFadTraits.hpp \
	$(srcdir)/Sacado_ELRFad_DFad.hpp \
	$(srcdir)/Sacado_ELRFad_SLFad.hpp \
	$(srcdir)/Sacado_ELRFad_SLFadTraits.hpp \
	$(srcdir)/Sacado_ELRFad_SFad.hpp \
	$(srcdir)/Sacado_ELRFad_SFadImp.hpp \
	$(srcdir)/Sacado_ELRFad_SFadTraits.hpp

SACADO_LFAD_H = \
	$(srcdir)/Sacado_LFad_LogicalSparse.hpp \
	$(srcdir)/Sacado_LFad_LogicalSparseImp.hpp \
	$(srcdir)/Sacado_LFad_LogicalSparseOps.hpp \
	$(srcdir)/Sacado_LFad_LogicalSparseTraits.hpp \
	$(srcdir)/Sacado_LFad_ExpressionTraits.hpp

SACADO_CACHEFAD_H = \
	$(srcdir)/Sacado_CacheFad_GeneralFad.hpp \
	$(srcdir)/Sacado_CacheFad_GeneralFadImp.hpp \
	$(srcdir)/Sacado_CacheFad_GeneralFadExpr.hpp \
	$(srcdir)/Sacado_CacheFad_Expression.hpp \
	$(srcdir)/Sacado_CacheFad_Ops.hpp \
	$(srcdir)/Sacado_CacheFad_DFadTraits.hpp \
	$(srcdir)/Sacado_CacheFad_DFad.hpp

SACADO_TAYLOR_H = \
	$(srcdir)/Sacado_Tay_CacheTaylorTraits.hpp \
	$(srcdir)/Sacado_Tay_CacheTaylor.hpp \
	$(srcdir)/Sacado_Tay_CacheTaylorImp.hpp \
	$(srcdir)/Sacado_Tay_CacheTaylorExpr.hpp \
	$(srcdir)/Sacado_Tay_CacheTaylorOps.hpp \
	$(srcdir)/Sacado_Tay_TaylorTraits.hpp \
	$(srcdir)/Sacado_Tay_Taylor.hpp \
	$(srcdir)/Sacado_Tay_TaylorImp.hpp

SACADO_RAD_H = \
	$(srcdir)/Sacado_trad_Traits.hpp \
	$(srcdir)/Sacado_trad.hpp \
	$(srcdir)/Sacado_trad2_Traits.hpp \
	$(srcdir)/Sacado_trad2.hpp \
	$(srcdir)/Sacado_tradvec_Traits.hpp \
	$(srcdir)/Sacado_tradvec.hpp \
	$(srcdir)/Sacado_rad.hpp \
	$(srcdir)/Sacado_rad2.hpp

SACADO_RAD = \
	$(srcdir)/Sacado_radops.cpp \
	$(srcdir)/Sacado_radops2.cpp \
	$(srcdir)/uninit.c \
	$(srcdir)/RAD_Const_Warn.c

SACADO_UTIL_H = \
	$(srcdir)/Sacado_DynamicArrayTraits.hpp \
	$(srcdir)/Sacado_StaticArrayTraits.hpp \
	$(srcdir)/Sacado_Handle.hpp

if ENABLE_SACADO_TEUCHOS

SACADO_BLAS_H = \
	$(srcdir)/Sacado_Fad_BLAS.hpp \
	$(srcdir)/Sacado_Fad_BLASImp.hpp \
	$(srcdir)/Sacado_Fad_ScalarTraitsImp.hpp

SACADO_TEMPLATE_H = \
	$(srcdir)/template/Sacado_TemplateManager.hpp \
	$(srcdir)/template/Sacado_TemplateManagerImp.hpp \
	$(srcdir)/template/Sacado_TemplateIterator.hpp

SACADO_PCE_H = \
	$(srcdir)/pce/Sacado_PCE_StandardPoly.hpp \
	$(srcdir)/pce/Sacado_PCE_StandardPolyImp.hpp \
	$(srcdir)/pce/Sacado_PCE_HermiteBasis.hpp \
	$(srcdir)/pce/Sacado_PCE_HermiteBasisImp.hpp \
	$(srcdir)/pce/Sacado_PCE_HermiteEBasis.hpp \
	$(srcdir)/pce/Sacado_PCE_HermiteEBasisImp.hpp \
	$(srcdir)/pce/Sacado_PCE_UnitHermiteBasis.hpp \
	$(srcdir)/pce/Sacado_PCE_UnitHermiteBasisImp.hpp \
	$(srcdir)/pce/Sacado_PCE_TripleProduct.hpp \
	$(srcdir)/pce/Sacado_PCE_TripleProductImp.hpp \
	$(srcdir)/pce/Sacado_PCE_Workspace.hpp \
	$(srcdir)/pce/Sacado_PCE_WorkspaceImp.hpp \
	$(srcdir)/pce/Sacado_PCE_Hermite.hpp \
	$(srcdir)/pce/Sacado_PCE_HermiteImp.hpp \
	$(srcdir)/pce/Sacado_PCE_HermiteTraits.hpp

if BUILD_PARAMLIB

SACADO_PARAMETER_H = \
	$(srcdir)/parameter/Sacado_ParameterFamilyBase.hpp \
	$(srcdir)/parameter/Sacado_ParameterFamilyBaseImp.hpp \
	$(srcdir)/parameter/Sacado_ParameterLibraryBase.hpp \
	$(srcdir)/parameter/Sacado_ParameterLibraryBaseImp.hpp \
	$(srcdir)/parameter/Sacado_ParameterVectorBase.hpp \
	$(srcdir)/parameter/Sacado_ScalarParameterEntry.hpp \
	$(srcdir)/parameter/Sacado_ScalarParameterFamily.hpp \
	$(srcdir)/parameter/Sacado_ScalarParameterLibrary.hpp \
	$(srcdir)/parameter/Sacado_ScalarParameterVector.hpp

endif

endif

if ENABLE_SACADO_STOKHOS

SACADO_SG_H = \
	$(srcdir)/pce/Sacado_PCE_OrthogPoly.hpp \
	$(srcdir)/pce/Sacado_PCE_OrthogPolyImp.hpp \
	$(srcdir)/pce/Sacado_PCE_OrthogPolyTraits.hpp

endif

SACADO_MPL_H = \
	$(srcdir)/mpl/Sacado_mpl_is_same.hpp \
	$(srcdir)/mpl/Sacado_mpl_disable_if.hpp \
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
	$(srcdir)/mpl/Sacado_mpl_vector_push_back_spec.hpp \
	$(srcdir)/mpl/Sacado_mpl_integral_c.hpp \
	$(srcdir)/mpl/Sacado_mpl_range_c.hpp \
	$(srcdir)/mpl/Sacado_mpl_apply.hpp \
	$(srcdir)/mpl/Sacado_mpl_apply_wrap.hpp \
	$(srcdir)/mpl/Sacado_mpl_bind.hpp \
	$(srcdir)/mpl/Sacado_mpl_has_type.hpp \
	$(srcdir)/mpl/Sacado_mpl_lambda.hpp \
	$(srcdir)/mpl/Sacado_mpl_placeholders.hpp \
	$(srcdir)/mpl/Sacado_mpl_quote.hpp \
	$(srcdir)/mpl/Sacado_mpl_type_wrap.hpp \
	$(srcdir)/mpl/Sacado_mpl_is_placeholder.hpp

SACADO_SFC_H = \
	$(srcdir)/Sacado_ScalarFlopCounter.hpp \
	$(srcdir)/Sacado_ScalarFlopCounterTraits.hpp

SACADO_SFC = \
	$(srcdir)/Sacado_ScalarFlopCounter.cpp

SACADO_H = \
	$(srcdir)/Sacado_ConfigDefs.h \
	$(srcdir)/Sacado_Version.hpp \
	$(srcdir)/Sacado_Traits.hpp \
	$(srcdir)/Sacado.hpp \
	$(srcdir)/Sacado_MathFunctions.hpp \
	$(srcdir)/Sacado_cmath.hpp \
	$(SACADO_FAD_H) \
	$(SACADO_ELRFAD_H) \
	$(SACADO_LFAD_H) \
	$(SACADO_CACHEFAD_H) \
	$(SACADO_TAYLOR_H) \
	$(SACADO_RAD_H) \
	$(SACADO_MPL_H) \
	$(SACADO_BLAS_H) \
	$(SACADO_TEMPLATE_H) \
	$(SACADO_PCE_H) \
	$(SACADO_SG_H) \
	$(SACADO_PARAMETER_H) \
	$(SACADO_SFC_H) \
	$(SACADO_UTIL_H)

SACADO = \
	$(SACADO_RAD) \
	$(SACADO_SFC)

include $(top_builddir)/Makefile.export.sacado

# ------------------------------------------------------------------------
# Sacado library specifications
# ------------------------------------------------------------------------
SACADO_LIB = libsacado.a

libsacado_a_SOURCES = $(SACADO) $(SACADO_H)

libsacado_a_CXXFLAGS = $(SACADO_INCLUDES)

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





