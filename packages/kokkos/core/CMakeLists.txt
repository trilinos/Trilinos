
INCLUDE(TribitsSubPackageMacros)

TRIBITS_SUBPACKAGE(Core)

# placeholder for future device...
TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_Winthread
  KOKKOS_HAVE_WINTHREAD
  "Enable Winthread support in Kokkos."
  "${TPL_ENABLE_Winthread}"
  )

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(unit_test)
TRIBITS_ADD_TEST_DIRECTORIES(perf_test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(usecases)

TRIBITS_SUBPACKAGE_POSTPROCESS()

