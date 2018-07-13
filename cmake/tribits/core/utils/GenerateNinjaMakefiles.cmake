set(_GenerateNinjaMakefiles_DIR ${CMAKE_CURRENT_LIST_DIR})
function(generate_ninja_makefiles srcdir)
  # Recurse into subdirectories.
  get_property(dirs DIRECTORY ${srcdir} PROPERTY SUBDIRECTORIES)
  foreach(d IN LISTS dirs)
    generate_ninja_makefiles(${d})
  endforeach()

  # Get the targets from this directory.
  get_property(targets DIRECTORY ${srcdir} PROPERTY BUILDSYSTEM_TARGETS)

  # Accumulate a global list of all targets.
  set_property(GLOBAL APPEND PROPERTY _GenerateNinjaMakefiles_TARGETS ${targets})

  # Compute paths to/from this directory and the top.
  get_property(bindir DIRECTORY ${srcdir} PROPERTY BINARY_DIR)
  if("${bindir}" STREQUAL "${CMAKE_BINARY_DIR}")
    set(topdir ".")
    set(subdir ".")
    # The top level provides all targets.
    get_property(targets GLOBAL PROPERTY _GenerateNinjaMakefiles_TARGETS)
  else()
    file(RELATIVE_PATH subdir ${CMAKE_BINARY_DIR} ${bindir})
    string(REGEX REPLACE "[^/]+" ".." topdir "${subdir}")
  endif()

  # Write the Makefile for this directory.
  string(REPLACE ";" " " NINJA_MAKEFILE_TARGETS "${targets}")
  file(TO_NATIVE_PATH "${topdir}" NINJA_MAKEFILE_TOPDIR)
  file(TO_NATIVE_PATH "${subdir}" NINJA_MAKEFILE_SUBDIR)
  file(TO_NATIVE_PATH "${_GenerateNinjaMakefiles_DIR}/NinjaMakefileCommon.make" NINJA_MAKEFILE_COMMON)
  configure_file(${_GenerateNinjaMakefiles_DIR}/NinjaMakefile.in ${bindir}/Makefile)
endfunction()
