include_tribits_build()

bob_config_header("${CMAKE_CURRENT_BINARY_DIR}/Compadre_Config.h")

MACRO(COMPADRE_APPEND_GLOB VAR)
  FILE(GLOB LOCAL_TMP_VAR ${ARGN})
  LIST(APPEND ${VAR} ${LOCAL_TMP_VAR})
ENDMACRO()

set(COMPADRE_HEADERS "")
set(COMPADRE_SOURCES "")
set(COMPADRE_HEADERS ${COMPADRE_HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/Compadre_Config.h
  )
compadre_append_glob(COMPADRE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
list(REMOVE_ITEM COMPADRE_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/Compadre_Manifold_Functions.hpp
  )
compadre_append_glob(COMPADRE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/basis/*.hpp)
compadre_append_glob(COMPADRE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/constraints/*.hpp)
compadre_append_glob(COMPADRE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/tpl/*.hpp)
compadre_append_glob(COMPADRE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_library(compadre ${COMPADRE_SOURCES})
bob_library_includes(compadre)

# allows us to use flat directory includes when building, since that will be the file structure once installed
target_include_directories(compadre PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/basis>)
target_include_directories(compadre PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/constraints>)
target_include_directories(compadre PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/tpl>)

# link to Kokkos
if (KOKKOS_EXISTING_ELSEWHERE) # user points to Kokkos installation
  target_link_libraries(compadre PUBLIC Kokkos::kokkoscore)
elseif (KOKKOS_IN_TRILINOS) # Kokkos found in Trilinos installation
  target_include_directories(compadre PUBLIC ${Trilinos_INCLUDE_DIRS})
  if(Trilinos_VERSION VERSION_LESS "13.5")
    target_link_libraries(compadre PUBLIC kokkoscore)
  elseif(Trilinos_VERSION VERSION_LESS "14.1")
    target_link_libraries(compadre PUBLIC KokkosCore::all_libs)
  else()
    target_link_libraries(compadre PUBLIC Kokkos::all_libs)
  endif()
elseif (KOKKOS_BUILT_FOR_USER) # Kokkos built for user by this project
  target_link_libraries(compadre PUBLIC Kokkos::kokkoscore)
endif()

if (KOKKOSKERNELS_EXISTING_ELSEWHERE) # user points to KokkosKernels installation
  target_link_libraries(compadre PUBLIC Kokkos::kokkoskernels)
elseif (KOKKOSKERNELS_IN_TRILINOS) # KokkosKernels found in Trilinos installation
  if(Trilinos_VERSION VERSION_LESS "13.5")
    target_link_libraries(compadre PUBLIC kokkoskernels)
  else()
    target_link_libraries(compadre PUBLIC KokkosKernels::all_libs)
  endif()
elseif (KOKKOSKERNELS_BUILT_FOR_USER) # Kokkos-Kernels built for user by this project
  target_link_libraries(compadre PUBLIC Kokkos::kokkoskernels)
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
    tpl/KokkosBatched_SolveUTV_Decl_Compadre.hpp
    tpl/KokkosBatched_SolveUTV_TeamVector_Impl_Compadre.hpp
    tpl/KokkosBatched_SolveUTV_TeamVector_Internal_Compadre.hpp
)

install(FILES ${COMPADRE_TPL} DESTINATION include/tpl)

install(FILES ${COMPADRE_HEADERS} DESTINATION include)

bob_end_subdir()
