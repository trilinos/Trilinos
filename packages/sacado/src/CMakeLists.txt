
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

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
  Sacado_Fad_FadVector.hpp
  Sacado_Fad_FadVectorImp.hpp
  Sacado_Fad_GeneralVFad.hpp
  Sacado_Fad_GeneralVFadExpr.hpp
  Sacado_Fad_GeneralVFadImp.hpp
  Sacado_Fad_VectorDynamicStorage.hpp
  Sacado_Fad_SimpleFad.hpp
  Sacado_Fad_SimpleFadImp.hpp
  Sacado_Fad_SimpleFadOps.hpp
  Sacado_Fad_SimpleFadTraits.hpp
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
  Sacado_CacheFad_Expression.hpp
  Sacado_CacheFad_Ops.hpp
  Sacado_CacheFad_DFadTraits.hpp
  Sacado_CacheFad_DFad.hpp
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
  Sacado_trad2.hpp
  Sacado_tradvec.hpp
  Sacado_rad.hpp
  Sacado_rad2.hpp
  )

# RAD Sources
SET(SOURCES ${SOURCES}
  Sacado_radops.cpp
  Sacado_radops2.cpp
  uninit.c
  RAD_Const_Warn.c
  )

# UTILS
SET(HEADERS ${HEADERS}
  Sacado_DynamicArrayTraits.hpp
  Sacado_StaticArrayTraits.hpp
  Sacado_Handle.hpp
  )

IF (${PACKAGE_NAME}_ENABLE_Teuchos)

  # TEMPLATE
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/template)
  SET(HEADERS ${HEADERS}
    template/Sacado_TemplateManager.hpp
    template/Sacado_TemplateManagerImp.hpp
    template/Sacado_TemplateIterator.hpp
    )

  # PCE
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/pce)
  SET(HEADERS ${HEADERS}
    pce/Sacado_PCE_StandardPoly.hpp
    pce/Sacado_PCE_StandardPolyImp.hpp
    pce/Sacado_PCE_HermiteBasis.hpp
    pce/Sacado_PCE_HermiteBasisImp.hpp
    pce/Sacado_PCE_HermiteEBasis.hpp
    pce/Sacado_PCE_HermiteEBasisImp.hpp
    pce/Sacado_PCE_UnitHermiteBasis.hpp
    pce/Sacado_PCE_UnitHermiteBasisImp.hpp
    pce/Sacado_PCE_TripleProduct.hpp
    pce/Sacado_PCE_TripleProductImp.hpp
    pce/Sacado_PCE_Workspace.hpp
    pce/Sacado_PCE_WorkspaceImp.hpp
    pce/Sacado_PCE_Hermite.hpp
    pce/Sacado_PCE_HermiteImp.hpp
    pce/Sacado_PCE_HermiteTraits.hpp
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
      )

  ENDIF(${PACKAGE_NAME}_ENABLE_PARAMLIB)
  
ENDIF(${PACKAGE_NAME}_ENABLE_Teuchos)

# Stokhos
IF (${PACKAGE_NAME}_ENABLE_Stokhos)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/pce)
  SET(HEADERS ${HEADERS}
    pce/Sacado_PCE_OrthogPoly.hpp
    pce/Sacado_PCE_OrthogPolyImp.hpp
    pce/Sacado_PCE_OrthogPolyTraits.hpp
    )
ENDIF()

#mpl
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/mpl)
SET(HEADERS ${HEADERS}
  mpl/Sacado_mpl_is_same.hpp
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
  Sacado_MathFunctions.hpp
  Sacado_cmath.hpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  sacado
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
