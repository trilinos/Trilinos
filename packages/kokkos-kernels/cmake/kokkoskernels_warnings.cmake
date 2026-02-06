macro(kokkoskernels_enable_warnings)
  set(COMMON_WARNINGS
      "-Wall"
      "-Wextra"
      "-Wextra-semi"
      "-Wunused-parameter"
      "-Wshadow"
      "-pedantic"
      "-Wsign-compare"
      "-Wtype-limits"
      "-Wuninitialized"
      "-Wsuggest-override"
  )

  # NOTE KOKKOS_ prefixed variable (all uppercase) is not set yet because TPLs are processed after ARCH
  if (Kokkos_ENABLE_LIBQUADMATH)
    # warning: non-standard suffix on floating constant [-Wpedantic]
    list(REMOVE_ITEM COMMON_WARNINGS "-pedantic")
  endif()

  # NVHPC compiler does not support -Wtype-limits.
  if (KOKKOS_ENABLE_OPENACC)
    if (Kokkos_CXX_COMPILER_ID STREQUAL NVHPC)
      list(REMOVE_ITEM COMMON_WARNINGS "-Wtype-limits")
    endif()
  endif()

  # nvcc raises internal warnings about extra semicolons
  if (Kokkos_CXX_COMPILER_ID STREQUAL NVIDIA)
    list(REMOVE_ITEM COMMON_WARNINGS "-Wextra-semi")
  endif()

  if (Kokkos_CXX_COMPILER_ID STREQUAL Clang)
    list(APPEND COMMON_WARNINGS "-Wimplicit-fallthrough")
  endif()

  set(GNU_WARNINGS "-Wempty-body" "-Wignored-qualifiers" ${COMMON_WARNINGS})
  if (Kokkos_CXX_COMPILER_ID STREQUAL GNU)
    list(APPEND GNU_WARNINGS "-Wimplicit-fallthrough")
  endif()

  # Not using COMPILER_SPECIFIC_FLAGS function so the warning flags are not passed downstream
  if (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
    string(REPLACE ";" " " WARNING_FLAGS "${GNU_WARNINGS}")
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL NVHPC)
    # FIXME_NVHPC
  else()
    string(REPLACE ";" " " WARNING_FLAGS "${COMMON_WARNINGS}")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_FLAGS}")

  if (KokkosKernels_ENABLE_WARNINGS_AS_ERRORS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
  endif()
endmacro()
