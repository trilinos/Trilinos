
INCLUDE(TribitsSubPackageMacros)

#
# A) Define the (sub)package
#

TRIBITS_SUBPACKAGE(Compat)

#
# B) Set up package-specific options
#

TRIBITS_ADD_DEBUG_OPTION()

# TRIBITS_ADD_OPTION_AND_DEFINE( 
#   KokkosCompat_ENABLE_TBB
#   HAVE_KOKKOSCOMPAT_TBB
#   "Enable TBB support in KokkosCompat."
#   ${${PROJECT_NAME}_ENABLE_TBB}
#   )

# TRIBITS_ADD_OPTION_AND_DEFINE( 
#   KokkosCompat_ENABLE_Cusp
#   HAVE_KOKKOSCOMPAT_CUSP
#   "Enable Cusp support in KokkosCompat."
#   ${${PROJECT_NAME}_ENABLE_CUSP}
#   )

# TRIBITS_ADD_OPTION_AND_DEFINE( 
#   KokkosCompat_ENABLE_MKL
#   HAVE_KOKKOSCOMPAT_MKL
#   "Enable MKL support in KokkosCompat."
#   ${${PROJECT_NAME}_ENABLE_MKL}
#   )

#ADD_SUBDIRECTORY(src)
#TRIBITS_ADD_TEST_DIRECTORIES(perf_test)
#TRIBITS_ADD_TEST_DIRECTORIES(unit_test)
#TRIBITS_ADD_EXAMPLE_DIRECTORIES(usecases)

TRIBITS_SUBPACKAGE_POSTPROCESS()
