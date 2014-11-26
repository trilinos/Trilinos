#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# FAD
SET(HEADERS ${HEADERS}
  Sacado_Fad_GeneralFad.hpp
  Sacado_Fad_GeneralFadImp.hpp
  Sacado_Fad_GeneralFadExpr.hpp
  Sacado_Fad_ExpressionTraits.hpp
  Sacado_Fad_Expression.hpp
  Sacado_Fad_Ops.hpp
  Sacado_Fad_DynamicStorage.hpp
  Sacado_Fad_DFad.hpp
  Sacado_Fad_DFadTraits.hpp
  Sacado_Fad_StaticStorage.hpp
  Sacado_Fad_SLFad.hpp
  Sacado_Fad_SLFadTraits.hpp
  Sacado_Fad_SFad.hpp
  Sacado_Fad_SFadImp.hpp
  Sacado_Fad_SFadTraits.hpp
  Sacado_Fad_MemPool.hpp
  Sacado_Fad_MemPoolImp.hpp
  Sacado_Fad_MemPoolManager.hpp
  Sacado_Fad_MemPoolManagerImp.hpp
  Sacado_Fad_MemPoolStorage.hpp
  Sacado_Fad_DMFad.hpp
  Sacado_Fad_DMFadTraits.hpp
  Sacado_Fad_DVFad.hpp
  Sacado_Fad_DVFadTraits.hpp
  Sacado_Fad_Vector.hpp
  Sacado_Fad_VectorImp.hpp
  Sacado_Fad_VectorDynamicStorage.hpp
  Sacado_Fad_SimpleFad.hpp
  Sacado_Fad_SimpleFadImp.hpp
  Sacado_Fad_SimpleFadOps.hpp
  Sacado_Fad_SimpleFadTraits.hpp
  Sacado_Fad_ViewStorage.hpp
  Sacado_Fad_ViewFad.hpp
  Sacado_Fad_ViewFadTraits.hpp
  Sacado_Fad_DFad_tmpl.hpp
  Sacado_Fad_SFad_tmpl.hpp
  Sacado_Fad_SLFad_tmpl.hpp
  Sacado_Fad_ViewFad_tmpl.hpp
  )

# ELRFAD
SET(HEADERS ${HEADERS}
  Sacado_ELRFad_GeneralFad.hpp
  Sacado_ELRFad_GeneralFadImp.hpp
  Sacado_ELRFad_GeneralFadExpr.hpp
  Sacado_ELRFad_ExpressionTraits.hpp
  Sacado_ELRFad_Expression.hpp
  Sacado_ELRFad_Ops.hpp
  Sacado_ELRFad_DFadTraits.hpp
  Sacado_ELRFad_DFad.hpp
  Sacado_ELRFad_SLFad.hpp
  Sacado_ELRFad_SLFadTraits.hpp
  Sacado_ELRFad_SFad.hpp
  Sacado_ELRFad_SFadImp.hpp
  Sacado_ELRFad_SFadTraits.hpp
  Sacado_ELRFad_ViewFad.hpp
  Sacado_ELRFad_ViewFadTraits.hpp
  )

# LFAD
SET(HEADERS ${HEADERS}
  Sacado_LFad_LogicalSparse.hpp
  Sacado_LFad_LogicalSparseImp.hpp
  Sacado_LFad_LogicalSparseOps.hpp
  Sacado_LFad_LogicalSparseTraits.hpp
  Sacado_LFad_ExpressionTraits.hpp
  )

# CACHEFAD
SET(HEADERS ${HEADERS}
  Sacado_CacheFad_GeneralFad.hpp
  Sacado_CacheFad_GeneralFadImp.hpp
  Sacado_CacheFad_GeneralFadExpr.hpp
  Sacado_CacheFad_ExpressionTraits.hpp
  Sacado_CacheFad_Expression.hpp
  Sacado_CacheFad_Ops.hpp
  Sacado_CacheFad_DFadTraits.hpp
  Sacado_CacheFad_DFad.hpp
  Sacado_CacheFad_SFadTraits.hpp
  Sacado_CacheFad_SFad.hpp
  Sacado_CacheFad_SFadImp.hpp
  Sacado_CacheFad_SLFadTraits.hpp
  Sacado_CacheFad_SLFad.hpp
  Sacado_CacheFad_ViewFad.hpp
  Sacado_CacheFad_ViewFadTraits.hpp
  )

# ELRCacheFAD
SET(HEADERS ${HEADERS}
  Sacado_ELRCacheFad_GeneralFad.hpp
  Sacado_ELRCacheFad_GeneralFadImp.hpp
  Sacado_ELRCacheFad_GeneralFadExpr.hpp
  Sacado_ELRCacheFad_ExpressionTraits.hpp
  Sacado_ELRCacheFad_Expression.hpp
  Sacado_ELRCacheFad_Ops.hpp
  Sacado_ELRCacheFad_DFadTraits.hpp
  Sacado_ELRCacheFad_DFad.hpp
  Sacado_ELRCacheFad_SLFad.hpp
  Sacado_ELRCacheFad_SLFadTraits.hpp
  Sacado_ELRCacheFad_SFad.hpp
  Sacado_ELRCacheFad_SFadImp.hpp
  Sacado_ELRCacheFad_SFadTraits.hpp
  Sacado_ELRCacheFad_ViewFad.hpp
  Sacado_ELRCacheFad_ViewFadTraits.hpp
  )

#TAYLOR
SET(HEADERS ${HEADERS}
  Sacado_Tay_CacheTaylorTraits.hpp
  Sacado_Tay_CacheTaylor.hpp
  Sacado_Tay_CacheTaylorImp.hpp
  Sacado_Tay_CacheTaylorExpr.hpp
  Sacado_Tay_CacheTaylorOps.hpp
  Sacado_Tay_TaylorTraits.hpp
  Sacado_Tay_Taylor.hpp
  Sacado_Tay_TaylorImp.hpp
  )

# RAD
SET(HEADERS ${HEADERS}
  Sacado_trad_Traits.hpp
  Sacado_trad.hpp
  Sacado_trad2_Traits.hpp
  Sacado_trad2.hpp
  Sacado_tradvec_Traits.hpp
  Sacado_tradvec.hpp
  Sacado_rad.hpp
  Sacado_rad2.hpp
  )

# RAD Sources
SET(SOURCES ${SOURCES}
  Sacado_radops.cpp
  Sacado_radops2.cpp
  RAD_Const_Warn.c
  )
IF(Sacado_ENABLE_UNINIT)
  SET(SOURCES ${SOURCES}
    uninit.c
  )
ENDIF()

# UTILS
SET(HEADERS ${HEADERS}
  Sacado_DynamicArrayTraits.hpp
  Sacado_StaticArrayTraits.hpp
  Sacado_Handle.hpp
  )

IF (${PACKAGE_NAME}_ENABLE_Teuchos)

  # BLAS
  SET(HEADERS ${HEADERS}
    Sacado_Fad_BLAS.hpp
    Sacado_Fad_BLASImp.hpp
    Sacado_Fad_LAPACK.hpp
    Sacado_Fad_LAPACKImp.hpp
    Sacado_Fad_ScalarTraitsImp.hpp
    Sacado_Tay_ScalarTraitsImp.hpp
    )

  # TEMPLATE
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/template)
  SET(HEADERS ${HEADERS}
    template/Sacado_TemplateManager.hpp
    template/Sacado_TemplateManagerImp.hpp
    template/Sacado_TemplateIterator.hpp
    )

  IF (${PACKAGE_NAME}_ENABLE_PARAMLIB)

    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/parameter)
    SET(HEADERS ${HEADERS}
      parameter/Sacado_ParameterFamilyBase.hpp
      parameter/Sacado_ParameterFamilyBaseImp.hpp
      parameter/Sacado_ParameterLibraryBase.hpp
      parameter/Sacado_ParameterLibraryBaseImp.hpp
      parameter/Sacado_ParameterVectorBase.hpp
      parameter/Sacado_ScalarParameterEntry.hpp
      parameter/Sacado_ScalarParameterFamily.hpp
      parameter/Sacado_ScalarParameterLibrary.hpp
      parameter/Sacado_ScalarParameterVector.hpp
      parameter/Sacado_ParameterAccessor.hpp
      parameter/Sacado_ParameterAccessorImp.hpp
      parameter/Sacado_ParameterRegistration.hpp
      )

  ENDIF(${PACKAGE_NAME}_ENABLE_PARAMLIB)

ENDIF(${PACKAGE_NAME}_ENABLE_Teuchos)

#mpl
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/mpl)
SET(HEADERS ${HEADERS}
  mpl/Sacado_mpl_is_same.hpp
  mpl/Sacado_mpl_is_convertible.hpp
  mpl/Sacado_mpl_enable_if.hpp
  mpl/Sacado_mpl_disable_if.hpp
  mpl/Sacado_mpl_if.hpp
  mpl/Sacado_mpl_none.hpp
  mpl/Sacado_mpl_begin.hpp
  mpl/Sacado_mpl_end.hpp
  mpl/Sacado_mpl_next.hpp
  mpl/Sacado_mpl_size.hpp
  mpl/Sacado_mpl_at.hpp
  mpl/Sacado_mpl_deref.hpp
  mpl/Sacado_mpl_push_back.hpp
  mpl/Sacado_mpl_find.hpp
  mpl/Sacado_mpl_for_each.hpp
  mpl/Sacado_mpl_vector.hpp
  mpl/Sacado_mpl_vector_size_spec.hpp
  mpl/Sacado_mpl_vector_at_spec.hpp
  mpl/Sacado_mpl_vector_push_back_spec.hpp
  mpl/Sacado_mpl_integral_c.hpp
  mpl/Sacado_mpl_range_c.hpp
  mpl/Sacado_mpl_apply.hpp
  mpl/Sacado_mpl_apply_wrap.hpp
  mpl/Sacado_mpl_bind.hpp
  mpl/Sacado_mpl_has_type.hpp
  mpl/Sacado_mpl_lambda.hpp
  mpl/Sacado_mpl_placeholders.hpp
  mpl/Sacado_mpl_quote.hpp
  mpl/Sacado_mpl_type_wrap.hpp
  mpl/Sacado_mpl_is_placeholder.hpp
  )

# Scalar flop counters
SET(HEADERS ${HEADERS}
  Sacado_ScalarFlopCounter.hpp
  Sacado_ScalarFlopCounterTraits.hpp
  )
SET(SOURCES ${SOURCES}
  Sacado_ScalarFlopCounter.cpp
  )

# Main library
SET(HEADERS ${HEADERS}
  Sacado_ConfigDefs.h
  Sacado_Version.hpp
  Sacado_Traits.hpp
  Sacado.hpp
  Sacado_Kokkos.hpp
  Kokkos_View_Fad.hpp
  Kokkos_AnalyzeSacadoShape.hpp
  Sacado_DisableKokkosCuda.hpp
  Sacado_MathFunctions.hpp
  Sacado_cmath.hpp
  Sacado_dummy_arg.hpp
  Sacado_Base.hpp
  Sacado_SFINAE_Macros.hpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  sacado
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
