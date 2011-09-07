INCLUDE(PackageLibraryMacros)
INCLUDE(AddOptionAndDefine)

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINRARY_DIR})

SET(HEADERS ${HEADERS}
 ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
    Kokkos_ArrayBounds.hpp
    Kokkos_MDArrayIndexMapLeft_macros.hpp
    Kokkos_MDArrayIndexMapRight_macros.hpp
    Kokkos_MDArrayView_macros.hpp
    Kokkos_MemoryInfo.hpp
    Kokkos_MemoryView_macros.hpp
    Kokkos_MultiVectorView_macros.hpp
    Kokkos_Preprocessing_macros.hpp
    Kokkos_StaticAssert.hpp
    Kokkos_ValueView_macros
    Kokkos_ViewTracker.hpp    
  )

APPEND_SET(SOURCES
    Kokkos_ArrayBounds.cpp
    Kokkos_MemoryInfo.cpp
    Kokkos_ViewTracker.cpp
    Kokkos_ViewMirror.cpp
  )

PACKAGE_ADD_LIBRARY(
    kokkosimpl
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS kokkos
    )

# No longer needed; Trilinos does this at the top level.
# IF (Kokkos_ENABLE_MDArray_Cuda)
#   FIND_PACKAGE(CUDA)
# ENDIF()
