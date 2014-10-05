
#
# Define the subpackage
#

TRIBITS_SUBPACKAGE(Compat)

#
# Set up subpackage-specific configuration macros
#

IF(NOT ${PROJECT_NAME}_ENABLE_TBB)
  SET(${PROJECT_NAME}_ENABLE_TBB OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE( KokkosCompat_ENABLE_TBB
  HAVE_KOKKOSCOMPAT_TBB
  "Enable TBB support in KokkosCompat."
  ${${PROJECT_NAME}_ENABLE_TBB}
  )

IF(NOT ${PROJECT_NAME}_ENABLE_Cusp)
  SET(${PROJECT_NAME}_ENABLE_Cusp OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE( 
  KokkosCompat_ENABLE_Cusp
  HAVE_KOKKOSCOMPAT_CUSP
  "Enable Cusp support in KokkosCompat."
  ${${PROJECT_NAME}_ENABLE_Cusp}
  )

IF(NOT ${PROJECT_NAME}_ENABLE_MKL)
  SET(${PROJECT_NAME}_ENABLE_MKL OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE( 
  KokkosCompat_ENABLE_MKL
  HAVE_KOKKOSCOMPAT_MKL
  "Enable MKL support in KokkosCompat."
  ${${PROJECT_NAME}_ENABLE_MKL}
  )

#-----------------------------------------------------------------------------

ADD_SUBDIRECTORY(src)
TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()
