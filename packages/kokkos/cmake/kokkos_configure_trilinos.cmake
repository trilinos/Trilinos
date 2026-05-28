if(CMAKE_PROJECT_NAME STREQUAL "Trilinos")

  # Trilinos sets some TPL dependency options using different names.
  # We handle the translation here.
  # This can cause issues if Trilinos is built against an external Kokkos.

  ASSERT_DEFINED(Kokkos_ENABLE_quadmath)
  set(Kokkos_ENABLE_LIBQUADMATH ${Kokkos_ENABLE_quadmath} CACHE BOOL "Whether to enable the LIBQUADMATH library" FORCE)

  ASSERT_DEFINED(Kokkos_ENABLE_DLlib)
  set(Kokkos_ENABLE_LIBDL ${Kokkos_ENABLE_DLlib} CACHE BOOL "Whether to enable the LIBDL library" FORCE)

endif()
