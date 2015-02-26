INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddExecutableAndTest)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


APPEND_SET(HEADERS
  CellData.hpp
  Traits.hpp
  FactoryTraits.hpp
  evaluators/Evaluator_Constant.hpp
  evaluators/Evaluator_Constant_Def.hpp
  evaluators/Evaluator_Density.hpp
  evaluators/Evaluator_Density_Def.hpp
  evaluators/Evaluator_EnergyFlux_Fourier.hpp
  evaluators/Evaluator_EnergyFlux_Fourier_Def.hpp
  evaluators/Evaluator_FEInterpolation.hpp
  evaluators/Evaluator_FEInterpolation_Def.hpp
  evaluators/Evaluator_NonlinearSource.hpp
  evaluators/Evaluator_NonlinearSource_Def.hpp
  )


IF(PHX_KOKKOS_DEVICE_TYPE_CUDA2)
  IF( TPL_ENABLE_MPI )
      CUDA_INCLUDE_DIRECTORIES( ${MPI_BASE_DIR}/include )
  ENDIF()

  CUDA_INCLUDE_DIRECTORIES( "." )
  CUDA_INCLUDE_DIRECTORIES( ${KokkosCore_INCLUDE_DIRS}
                            ${Phalanx_INCLUDE_DIRS} )
  CUDA_COMPILE( OBJECTS_CUDA CellData_CU.cu )

  LIST( APPEND SOURCES ${OBJECTS_CUDA} )

#SET_PROPERTY(TARGET phalanx_test_utilities APPEND PROPERTY LINKER_LANGUAGE CXX)

 TRIBITS_ADD_LIBRARY(
  phalanx_test_utilities
  ${OBJECTS_CUDA}
  HEADERS ${HEADERS}
  TESTONLY
  )
 
 SET_PROPERTY(TARGET phalanx_test_utilities APPEND PROPERTY LINKER_LANGUAGE CXX)

ELSE()


APPEND_SET(SOURCES
  CellData.cpp
  )

#ENDIF ()

TRIBITS_ADD_LIBRARY(
  phalanx_test_utilities
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  TESTONLY
  )
ENDIF()
