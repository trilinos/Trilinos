if(CMAKE_PROJECT_NAME STREQUAL "Seacas" OR CMAKE_PROJECT_NAME STREQUAL "SEACAS" )
TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_REQUIRED_TPLS fmt
  LIB_OPTIONAL_PACKAGES SEACASExodus Zoltan
  LIB_OPTIONAL_TPLS HDF5 Pamgen CGNS ParMETIS Faodel Cereal DLlib Pthread ADIOS2 Catalyst2 ${SEACAS_GTest_TPL_name} Kokkos DataWarp Catch2
)
else()
# Typically for Trilinos since don't have fmt as TPL, but instead have embedded versions.
TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_OPTIONAL_PACKAGES SEACASExodus Pamgen Zoltan Kokkos
  LIB_OPTIONAL_TPLS HDF5 CGNS ParMETIS Faodel Cereal DLlib Pthread DataWarp ADIOS2 Catalyst2 ${SEACAS_GTest_TPL_name}
)
endif()

TRIBITS_TPL_TENTATIVELY_ENABLE(DLlib)
