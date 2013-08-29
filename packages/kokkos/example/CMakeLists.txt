
INCLUDE(TribitsSubPackageMacros)

# Subpackage name must match what appears in kokkos/cmake/Dependencies.cmake
#
TRIBITS_SUBPACKAGE(Example)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(query_device)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(fixture)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(feint)

TRIBITS_SUBPACKAGE_POSTPROCESS()

