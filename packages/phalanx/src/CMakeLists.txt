
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# FAD
APPEND_SET(HEADERS
  Phalanx.hpp
  Phalanx_Allocator_Contiguous.hpp
  Phalanx_Allocator_New.hpp
  Phalanx_ConfigDefs.hpp
  Phalanx_DataContainer_Base_Def.hpp
  Phalanx_DataContainer_Base.hpp
  Phalanx_DataContainer_Def.hpp
  Phalanx_DataContainer.hpp
  Phalanx_DataContainer_TemplateManager.hpp
  Phalanx_DataLayout_MDALayout_Def.hpp
  Phalanx_DataLayout_MDALayout.hpp
  Phalanx_DataLayout.hpp
  Phalanx_DoxygenDocumentation.hpp
  Phalanx_EvaluationContainer_Base_Def.hpp
  Phalanx_EvaluationContainer_Base.hpp
  Phalanx_EvaluationContainer_Def.hpp
  Phalanx_EvaluationContainer.hpp
  Phalanx_EvaluationContainer_TemplateBuilder.hpp
  Phalanx_EvaluationContainer_TemplateManager.hpp
  Phalanx_Evaluator_Base.hpp
  Phalanx_Evaluator_Derived.hpp
  Phalanx_Evaluator_Factory_Def.hpp
  Phalanx_Evaluator_Factory.hpp
  Phalanx_Evaluator_Factory_UFO.hpp
  Phalanx_Evaluator.hpp
  Phalanx_Evaluator_Macros.hpp
  Phalanx_Evaluator_Manager_Def.hpp
  Phalanx_Evaluator_Manager.hpp
  Phalanx_Evaluator_TemplateBuilder.hpp
  Phalanx_Evaluator_TemplateManager.hpp
  Phalanx_Evaluator_Utilities.hpp
  Phalanx_Evaluator_WithBaseImpl_Def.hpp
  Phalanx_Evaluator_WithBaseImpl.hpp
  Phalanx_ExpressionTemplates_Array.hpp
  Phalanx_ExpressionTemplates_ExprScalar.hpp
  Phalanx_ExpressionTemplates_Operands.hpp
  Phalanx_ExpressionTemplates_Operators.hpp
  Phalanx_ExpressionTemplates_Traits.hpp
  Phalanx_Field_Def.hpp
  Phalanx_Field.hpp
  Phalanx_FieldManager_Def.hpp
  Phalanx_FieldManager.hpp
  Phalanx_FieldTag.hpp
  Phalanx_FieldTag_STL_Functors.hpp
  Phalanx_FieldTag_Tag_Def.hpp
  Phalanx_FieldTag_Tag.hpp
  Phalanx_MDField_Def.hpp
  Phalanx_MDField.hpp
  Phalanx_Print_Utilities.hpp
  Phalanx_TemplateIterator.hpp
  Phalanx_TemplateManager_Def.hpp
  Phalanx_TemplateManager.hpp
  Phalanx_Traits_Base.hpp
  Phalanx_TypeStrings.hpp
  Phalanx_Version.hpp
  ) 

APPEND_SET(SOURCES
  Phalanx_DataLayout.cpp
  Phalanx_FieldTag.cpp
  )
#
# C) Define the targets for package's library(s)
#

IF (Phalanx_ENABLE_Thrust)
  # FINISH: REMOVE THIS IF IT ISN'T NEEDED
  # INCLUDE_DIRECTORIES(${Kokkos_INCLUDE_DIRS})
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/gpu)
  PACKAGE_ADD_CUDA_LIBRARY(
    phalanx_test_cuda 
    SOURCES gpu/Phalanx_GPU.cu
  )
  SET(Phalanx_CUDA_LIBRARY phalanx_test_cuda)
ELSE()
  SET(Phalanx_CUDA_LIBRARY "")
ENDIF()


PACKAGE_ADD_LIBRARY(
  phalanx
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS ${Phalanx_CUDA_LIBRARY}
  )

