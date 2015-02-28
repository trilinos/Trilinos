

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/evaluators)

IF(PHX_KOKKOS_DEVICE_TYPE_CUDA2)
  IF( TPL_ENABLE_MPI )
      CUDA_INCLUDE_DIRECTORIES( ${MPI_BASE_DIR}/include )
  ENDIF()
 

  CUDA_INCLUDE_DIRECTORIES( "." )
  CUDA_INCLUDE_DIRECTORIES( ${KokkosCore_INCLUDE_DIRS} 
                            ${Phalanx_INCLUDE_DIRS})
  CUDA_COMPILE( OBJECTS_CUDA_DFAD Example_DFAD_CU.cu )
  CUDA_COMPILE( OBJECTS_CUDA_SFAD Example_SFAD_CU.cu )
  CUDA_COMPILE( OBJECTS_CUDA_SLFAD Example_SLFAD_CU.cu )
  CUDA_COMPILE( CELL_OBJECTS_CUDA Cell_CU.cu )
#  LIST( APPEND example_SOURCES ${OBJECTS_CUDA})
SET(example_SOURCES
  FactoryTraits.hpp
  Dimension.hpp
  Dimension.cpp
  Workset.hpp
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

SET(example_DFAD_SOURCES ${OBJECTS_CUDA_DFAD} TraitsDFAD.hpp ${example_SOURCES})
SET(example_SFAD_SOURCES ${OBJECTS_CUDA_SFAD} TraitsSFAD.hpp ${example_SOURCES})
SET(example_SLFAD_SOURCES ${OBJECTS_CUDA_SFAD} TraitsSLFAD.hpp ${example_SOURCES})

ELSE()

SET(example_SOURCES
  FactoryTraits.hpp
  Cell.hpp
  Cell.cpp
  Dimension.hpp
  Dimension.cpp
  Workset.hpp
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

SET(example_DFAD_SOURCES Example_DFAD.cpp TraitsDFAD.hpp ${example_SOURCES})
SET(example_SFAD_SOURCES Example_SFAD.cpp TraitsSFAD.hpp ${example_SOURCES})
SET(example_SLFAD_SOURCES Example_SLFAD.cpp TraitsSLFAD.hpp ${example_SOURCES})

ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_dfad
  SOURCES ${example_DFAD_SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Run has completed successfully!"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_sfad
  SOURCES ${example_SFAD_SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Run has completed successfully!"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_slfad
  SOURCES ${example_SLFAD_SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Run has completed successfully!"
  )
