include_tribits_build()

bob_config_header("${CMAKE_CURRENT_BINARY_DIR}/Compadre_Config.h")

set(COMPADRE_HEADERS
  "${CMAKE_CURRENT_BINARY_DIR}/Compadre_Config.h"
  Compadre_Functors.hpp
  Compadre_GMLS.hpp
  Compadre_GMLS_ApplyTargetEvaluations.hpp
  Compadre_GMLS_Basis.hpp
  Compadre_GMLS_Targets.hpp
  Compadre_KokkosParser.hpp
  Compadre_LinearAlgebra_Declarations.hpp
  Compadre_LinearAlgebra_Definitions.hpp
  Compadre_Misc.hpp
  Compadre_Operators.hpp
  Compadre_ParallelManager.hpp
  Compadre_PointCloudSearch.hpp
  Compadre_Quadrature.hpp
  Compadre_Typedefs.hpp
  Compadre_Utilities.hpp
  Compadre_Evaluator.hpp
  USER_StandardTargetFunctionals.hpp
  USER_ManifoldTargetFunctionals.hpp
  basis/DivergenceFree3D.hpp
  basis/CreateConstraints.hpp
  )


set(COMPADRE_SOURCES
  Compadre_GMLS.cpp
  Compadre_KokkosParser.cpp
  Compadre_LinearAlgebra.cpp
  )

add_library(compadre ${COMPADRE_SOURCES})
bob_library_includes(compadre)

# openmp if needed
if (Compadre_USE_OPENMP)
  find_package(OpenMP QUIET)
  if(OpenMP_CXX_FOUND)
      target_link_libraries(compadre PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()

# link to Kokkos
if (KOKKOS_EXISTING_ELSEWHERE) # user points to Kokkos installation
  target_link_libraries(compadre PUBLIC Kokkos::kokkoscore)
elseif (KOKKOS_IN_TRILINOS) # Kokkos found in Trilinos installation
  target_include_directories(compadre PUBLIC ${Trilinos_INCLUDE_DIRS})
  target_link_libraries(compadre PUBLIC kokkoscore)
elseif (KOKKOS_BUILT_FOR_USER) # Kokkos built for user by this project
  target_link_libraries(compadre PUBLIC Kokkos::kokkoscore)
endif()

if (KOKKOSKERNELS_EXISTING_ELSEWHERE) # user points to KokkosKernels installation
  target_link_libraries(compadre PUBLIC Kokkos::kokkoskernels)
elseif (KOKKOSKERNELS_IN_TRILINOS) # KokkosKernels found in Trilinos installation
  target_link_libraries(compadre PUBLIC kokkoskernels)
elseif (KOKKOSKERNELS_BUILT_FOR_USER) # Kokkos-Kernels built for user by this project
  target_link_libraries(compadre PUBLIC Kokkos::kokkoskernels)
endif()

# link to CUDA
if (Compadre_USE_CUDA AND NOT(KOKKOS_IN_TRILINOS))
  if (TARGET Kokkos::CUDA)
    target_link_libraries(compadre PUBLIC "cusolver") # comes from Kokkos::CUBLAS
  else()
    MESSAGE(FATAL_ERROR "Compadre_USE_CUDA is ON, but CUBLAS not found by KokkosKernels.")
  endif()
endif()

# link to MPI
if (Compadre_USE_MPI)
  if (MPI_CXX_INCLUDE_PATH)
    target_include_directories(compadre PUBLIC "${MPI_CXX_INCLUDE_PATH}")
  endif()
  if (MPI_CXX_LIBRARIES)
    target_link_libraries(compadre PUBLIC "${MPI_CXX_LIBRARIES}")
  endif()
endif()

bob_export_target(compadre)

set(COMPADRE_TPL
    tpl/nanoflann.hpp
)

install(FILES ${COMPADRE_TPL} DESTINATION include/tpl)

install(FILES ${COMPADRE_HEADERS} DESTINATION include)

bob_end_subdir()
