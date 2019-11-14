#
# Trilinos configuration options that, togther with ATDMDevEnv.cmake create an
# install of Trilinos that works for EMPIRE (and only for EMPIRE).
#

ATDM_SET_CACHE(Gtest_SKIP_INSTALL TRUE CACHE BOOL)
# NOTE: Above, we have to set Gtest_SKIP_INSTALL=TRUE here just for EMPIRE for
# now and not in ATDMDevEnvSetting.cmake because SPARC currently must have
# gtest installed from Trilinos to build any tests.  See Trilinos GitHub
# #5341.  (Once we can get SPARC to have its own gtest, then we can move this
# setting to ATDMDevEnvSettings.cmake. See SPAR-614.)

ATDM_SET_CACHE(Trilinos_ENABLE_Gtest FALSE CACHE BOOL)

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/EMPIRETrilinosPackagesList.cmake")
FOREACH(TRIBITS_PACKAGE ${EMPIRE_Trilinos_Packages})
  SET(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} ON)
ENDFOREACH()
