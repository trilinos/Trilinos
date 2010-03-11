project(TrilinosDriver NONE)
cmake_minimum_required(VERSION 2.7.20090924)

macro(download_prebuilt_cmake _type)
  if(NOT _type)
    set(_type "release")
  endif(NOT _type)
  
  message("Downloading prebuilt CMake binaries...")
  execute_process(COMMAND python "${TrilinosDriver_SOURCE_DIR}/../../python/download-cmake.py"
                  "--install-dir=${TrilinosDriver_BINARY_DIR}/cmake-prebuilt/"
                  "--installer-type=${_type}"
                  WORKING_DIRECTORY "${TrilinosDriver_BINARY_DIR}"
                  ERROR_VARIABLE err
                  RESULT_VARIABLE result
                  OUTPUT_QUIET)
                  
  #file(REMOVE_RECURSE "${TrilinosDriver_SOURCE_DIR}/download_area")
                  
  if(${result} EQUAL 0)
    message("Finished downloading prebuilt CMake binaries.")
  else()
    message(FATAL_ERROR "Error occurred while downloading prebuilt CMake: ${err}")
  endif()
endmacro()

set(DOWNLOADED_DIR "${TrilinosDriver_SOURCE_DIR}/cmake-prebuilt/bin")
set(DOWNLOADED_CMAKE "${DOWNLOADED_DIR}/cmake")
set(DOWNLOADED_CTEST "${DOWNLOADED_DIR}/ctest")

# Figure out which machine we're on
site_name(site)
message("site=${site}")
set(subdir "${site}")
if("${site}" MATCHES "corrin")
  set(subdir "corrin")
endif()
  
add_subdirectory("${subdir}")


#add_subdirectory(corrin)
#add_subdirectory(dash2win64)
#add_subdirectory(exetazo)
#add_subdirectory(gabriel)
#add_subdirectory(godel)
#add_subdirectory(kallikrates)
#add_subdirectory(s903186)
#add_subdirectory(s909348)
#add_subdirectory(thumper)
#add_subdirectory(trilinos-test)
#add_subdirectory(trilinos-test2)
