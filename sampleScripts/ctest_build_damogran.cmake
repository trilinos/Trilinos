set(CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}")
set(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/Dashboards/Trilinos")
set(CTEST_BINARY_DIRECTORY "$ENV{HOME}/Dashboards/tb2")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_SITE "damogran.kitware")
set(CTEST_BUILD_NAME "MacOSX-make")
set(CTEST_BUILD_FLAGS -j3)
set(ENV{FC} "/Users/davidcole/Dashboards/Support/gfortran/bin/gfortran")


include(${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake)
include(${CTEST_SCRIPT_DIRECTORY}/TrilinosArrakisCTestConfig.cmake)
set(CTEST_USE_LAUNCHERS 1)


find_package(CVS)
set(CTEST_UPDATE_COMMAND "${CVS_EXECUTABLE}")
message("CTEST_UPDATE_COMMAND='${CTEST_UPDATE_COMMAND}'")


set(Trilinos_PACKAGES
#Pamgen
#ThreadPool
#Teuchos
#Shards
#Kokkos
${CTEST_PROJECT_SUBPROJECTS}
)
set(Trilinos_FAILED_PACKAGES) # store set of failed packages


ctest_empty_binary_directory("${CTEST_BINARY_DIRECTORY}")


#ctest_start("Experimental")
ctest_start("Nightly")
message("updating...")
ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)
message("updated... res='${res}'")


# Tell CDash about the latest subproject dependencies:
#
ctest_submit(FILES "${CTEST_SOURCE_DIRECTORY}/cmake/python/data/CDashSubprojectDependencies.xml"
  RETURN_VALUE res)
message("submitted subproject dependencies... res='${res}'")


# loop over all Trilinos packages
message("begin loop")

set(Trilinos_LAST_WORKING_PACKAGE)

foreach(TRIBITS_PACKAGE ${Trilinos_PACKAGES})
  set_property(GLOBAL PROPERTY SubProject ${TRIBITS_PACKAGE})
  set_property(GLOBAL PROPERTY Label ${TRIBITS_PACKAGE})
  message("TRIBITS_PACKAGE='${TRIBITS_PACKAGE}'")

  # CONFIGURE STEP

  # create CONFIGURE_OPTIONS for this TRIBITS_PACKAGE
  set(CONFIGURE_OPTIONS
    "-DCTEST_USE_LAUNCHERS:BOOL=${CTEST_USE_LAUNCHERS}"
    "-DTrilinos_ENABLE_${TRIBITS_PACKAGE}:BOOL=ON"
    "-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON"
    "-DTrilinos_ENABLE_TESTS:BOOL=ON")
  if(DEFINED Trilinos_LAST_WORKING_PACKAGE)
    list(APPEND CONFIGURE_OPTIONS
      "-DTrilinos_ENABLE_${Trilinos_LAST_WORKING_PACKAGE}:BOOL=")
    set(Trilinos_LAST_WORKING_PACKAGE)
  endif()
  foreach(FAILED_PACKAGE ${Trilinos_FAILED_PACKAGES})
    list(APPEND CONFIGURE_OPTIONS
      "-DTrilinos_ENABLE_${FAILED_PACKAGE}:BOOL=OFF")
  endforeach(FAILED_PACKAGE)

  # new option for ctest_configure OPTIONS a list of command
  # line arguments to give to cmake: -D stuff for example
  message("Configure TRIBITS_PACKAGE='${TRIBITS_PACKAGE}'")
  message("CONFIGURE_OPTIONS = '${CONFIGURE_OPTIONS}'")
  ctest_configure(
    BUILD "${CTEST_BINARY_DIRECTORY}"
    OPTIONS "${CONFIGURE_OPTIONS}"
    RETURN_VALUE res
    )

  # if the configure failed add the package to the list
  # of failed packages
  if(NOT "${res}" EQUAL "0")
    message("${TRIBITS_PACKAGE} FAILED to configure")
    list(APPEND Trilinos_FAILED_PACKAGES ${TRIBITS_PACKAGE})
  else()
    # load target properties and test keywords
    ctest_read_custom_files(BUILD "${CTEST_BINARY_DIRECTORY}")
  endif()

  include(${CTEST_SCRIPT_DIRECTORY}/TrilinosArrakisCTestConfig.cmake)

  # submit configure results
  ctest_submit( PARTS configure notes ) # submit notes and configure at this point


  # If configure failed, move to next iteration:
  #
  if("${res}" EQUAL "0")

    set(CTEST_BUILD_TARGET ${TRIBITS_PACKAGE}_libs)
    message("Build: '${CTEST_BUILD_TARGET}'")
    ctest_build (
      BUILD "${CTEST_BINARY_DIRECTORY}"
      RETURN_VALUE res
      NUMBER_ERRORS numerrors
      APPEND
      )
    # set this to false
    set(build_success FALSE)
    # since make -i is used the res might be 0, but
    # if there are errors the build should fail, so
    # both res and numerrors should be 0 for a good build
    # and for the all target to be built and tests run
    if("${numerrors}" EQUAL "0" AND "${res}" EQUAL "0")
      set(build_success TRUE)
    endif()
    ctest_submit( PARTS build ) # submit the build
    message("return from ctest_build ${res}")

    # check to see if the build worked
    if(NOT build_success)
      message("FAILED: build '${TRIBITS_PACKAGE}'")
      list(APPEND Trilinos_FAILED_PACKAGES ${TRIBITS_PACKAGE})
    else()
      # inside here the build worked
      # Now build ALL target, but append the results to the last build.xml
      set(CTEST_BUILD_TARGET)
      message("build all for '${TRIBITS_PACKAGE}'")
      ctest_build (BUILD "${CTEST_BINARY_DIRECTORY}"
        NUMBER_ERRORS numerrors
        RETURN_VALUE res
        APPEND
        )
      # submit the build for all target
      ctest_submit( PARTS build )  
      # now run the tests that match the ${TRIBITS_PACKAGE} name 
      message("test for '${TRIBITS_PACKAGE}'")
      ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}"
        INCLUDE "^${TRIBITS_PACKAGE}" APPEND
        )
      #ctest_memcheck(BUILD "${CTEST_BINARY_DIRECTORY}"
      #  KEYWORDS ${TRIBITS_PACKAGE})
      #ctest_coverage(BUILD "${CTEST_BINARY_DIRECTORY}")
      # submit test, memcheck and coverage results
      ctest_submit(PARTS Test)
      set(Trilinos_LAST_WORKING_PACKAGE "${TRIBITS_PACKAGE}")
    endif()
  endif()
endforeach(TRIBITS_PACKAGE)

if(Trilinos_FAILED_PACKAGES)
  message("Failed packages! ${Trilinos_FAILED_PACKAGES}")
endif()

message("Done with build.")
