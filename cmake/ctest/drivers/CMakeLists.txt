project(TrilinosDriver NONE)
cmake_minimum_required(VERSION 2.7.20090924)

include(CTest)
enable_testing()

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
  
  set(DOWNLOADED_DIR "${TrilinosDriver_BINARY_DIR}/cmake-prebuilt/bin")
  set(DOWNLOADED_CMAKE "${DOWNLOADED_DIR}/cmake")
  set(DOWNLOADED_CTEST "${DOWNLOADED_DIR}/ctest")
endmacro()

set(SITES
corrin
dash2win64
exetazo
gabriel
godel
kallikrates
s903186
s909348
thumper
trilinos-test
trilinos-test2
)

# Figure out which machine we're on
site_name(site)
set(subdir "${site}")

foreach(sitename ${SITES})
  if("${site}" MATCHES "${sitename}")
    set(subdir "${sitename}")
  endif()
endforeach()
  
add_subdirectory("${subdir}")
