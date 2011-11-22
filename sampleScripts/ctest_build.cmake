set(CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}")
set(CTEST_ROOT_DIRECTORY "$ENV{HOME}/Dashboards/MyTests")
#set(CTEST_ROOT_DIRECTORY "$ENV{HOME}/Projects/Trilinos")
set(CTEST_SOURCE_DIRECTORY "${CTEST_ROOT_DIRECTORY}/Trilinos")
set(CTEST_BINARY_DIRECTORY "${CTEST_ROOT_DIRECTORY}/Trilinos-build")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_SITE "corrin.kitware")
set(CTEST_BUILD_NAME "Windows-make")
set(CTEST_BUILD_FLAGS -j2)

include(${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake)


set(Trilinos_PACKAGES
 ${CTEST_PROJECT_SUBPROJECTS}
 )
set(Trilinos_FAILED_PACKAGES) # store set of failed packages


ctest_empty_binary_directory("${CTEST_BINARY_DIRECTORY}")


ctest_start("Experimental")
set(CTEST_DROP_SITE "arrakis.kitwarein.com")
set(CTEST_DROP_LOCATION "/CDash/submit.php?project=Trilinos")
set(CTEST_DROP_SITE_CDASH TRUE)

# loop over all Trilinos packages
message("begin loop")
ctest_submit(FILES 
  ${CTEST_SOURCE_DIRECTORY}/cmake/python/data/CDashSubprojectDependencies.xml)

foreach(TRIBITS_PACKAGE ${Trilinos_PACKAGES})
  set_property(GLOBAL PROPERTY SubProject ${TRIBITS_PACKAGE})
  set_property(GLOBAL PROPERTY Label ${TRIBITS_PACKAGE})
  message("TRIBITS_PACKAGE='${TRIBITS_PACKAGE}'")
  
  # CONFIGURE STEP
  
  # create CONFIGURE_OPTIONS for this TRIBITS_PACKAGE
  set(CONFIGURE_OPTIONS
    "-DTrilinos_ENABLE_${TRIBITS_PACKAGE}:BOOL=ON"
    "-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON"
    "-DTrilinos_ENABLE_TESTS:BOOL=ON")
  foreach(FAILED_PACKAGE ${Trilinos_FAILED_PACKAGES})
    list(APPEND CONFIGURE_OPTIONS
      "-DTrilinos_ENABLE_${FAILED_PACKAGE}:BOOL=OFF")
  endforeach(FAILED_PACKAGE)

  # new option for ctest_configure OPTIONS a list of command
  # line arguments to give to cmake: -D stuff for example
  ctest_configure(
    BUILD "${CTEST_BINARY_DIRECTORY}"
    OPTIONS "${CONFIGURE_OPTIONS}"
    RETURN_VALUE res
    )

  # if the configure failed add the package to the list
  # of failed packages
  if(NOT "${res}" EQUAL "0")
    list(APPEND FAILED_PACKAGE ${TRIBITS_PACKAGE})
  else()
    # load target properties and test keywords
    ctest_read_custom_files(BUILD "${CTEST_BINARY_DIRECTORY}")
  endif()
  set(CTEST_DROP_SITE "arrakis.kitwarein.com")

  # submit configure results
  ctest_submit(  ) # submit notes and configure at this point

  # If configure passed, build and test
  if("${res}" EQUAL "0")
    # set the target to build, build it, and submit it
    set(CTEST_BUILD_TARGET ${TRIBITS_PACKAGE}_libs)
    ctest_build (BUILD "${CTEST_BINARY_DIRECTORY}"
      RETURN_VALUE res  APPEND)
    ctest_submit( PARTS build ) # submit the build

    # if the build failed add it to the failed packages
    if(NOT "${res}" EQUAL "0")
      list(APPEND FAILED_PACKAGE ${TRIBITS_PACKAGE})
    else()
      # inside here the build worked
      # Now build ALL target, but append the results to the last build.xml
      set(CTEST_BUILD_TARGET)
      ctest_build (BUILD "${CTEST_BINARY_DIRECTORY}"
        RETURN_VALUE res APPEND )
      # submit the build for all target
      # submit should take a list of things to submit
      ctest_submit( PARTS build )  
      # now run the tests that match the ${TRIBITS_PACKAGE} name
      ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}"
        INCLUDE ${TRIBITS_PACKAGE}
        )
      #ctest_memcheck(BUILD "${CTEST_BINARY_DIRECTORY}"
      #  KEYWORDS ${TRIBITS_PACKAGE})
      #ctest_coverage(BUILD "${CTEST_BINARY_DIRECTORY}")
      # submit test, memcheck and coverage results
      ctest_submit(PARTS Test)
    endif()
  endif()
endforeach(TRIBITS_PACKAGE)

message("end loop")
