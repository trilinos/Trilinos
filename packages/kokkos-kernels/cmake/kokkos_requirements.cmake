
if(NOT KOKKOSKERNELS_HAS_TRILINOS AND NOT KOKKOSKERNELS_HAS_PARENT)

  if(DEFINED KokkosKernels_REQUIRE_DEVICES)
    string(REPLACE "," ";" REQUIRE_DEVICES "${KokkosKernels_REQUIRE_DEVICES}")
    kokkos_check(DEVICES ${REQUIRE_DEVICES})
  endif()

  if(DEFINED KokkosKernels_REQUIRE_OPTIONS)
    string(REPLACE "," ";" REQUIRE_OPTIONS "${KokkosKernels_REQUIRE_OPTIONS}")
    kokkos_check(OPTIONS ${REQUIRE_OPTIONS})
  endif()

  if(DEFINED KokkosKernels_REQUIRE_ARCH)
    string(REPLACE "," ";" REQUIRE_ARCH "${KokkosKernels_REQUIRE_ARCH}")
    kokkos_check(ARCH ${REQUIRE_ARCH})
  endif()

  if(DEFINED KokkosKernels_REQUIRE_TPLS)
    string(REPLACE "," ";" REQUIRE_TPLS "${KokkosKernels_REQUIRE_TPLS}")
    kokkos_check(TPLS ${REQUIRE_TPLS})
  endif()

endif()
