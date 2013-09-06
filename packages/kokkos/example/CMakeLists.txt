
INCLUDE(TribitsSubPackageMacros)

# Subpackage name must match what appears in kokkos/cmake/Dependencies.cmake
#
TRIBITS_SUBPACKAGE(Example)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(query_device)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(fixture)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(feint)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(md_skeleton)

TRIBITS_SUBPACKAGE_POSTPROCESS()

