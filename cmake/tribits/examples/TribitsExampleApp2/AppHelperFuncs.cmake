include(CMakePrintHelpers)


# Find TribitsExProj2 package(s), load compilers and compiler options, and get
# CMake lib targets (which have include dirs also) to link against.
#
# On return, sets the vars in the current scope:
#
# * TribitsExProj2_SELECTED_PACKAGE_LIST: List of all of the packages pulled in
# * from the TribitsExProj2.
#
# * APP_DEPS_LIB_TARGETS: List of all of the IMPORTED CMake targets that 'app'
#   must link against
#
# * CMAKE_<LANG>_COMPILER and CMAKE_<LANG>_FLAGS pulled in from the
#   TribitsExProj2Config.config or a <Package>Config.cmake file
#
macro(getTribitsExProj2StuffForApp)

  set(${PROJECT_NAME}_FIND_INDIVIDUAL_PACKAGES OFF CACHE BOOL
    "Set to TRUE to find individual packages and OFF to find project TribitsExProj2")

  if (${PROJECT_NAME}_FIND_INDIVIDUAL_PACKAGES)
    getTribitsExProj2StuffForAppByPackage()
  else()
    getTribitsExProj2StuffForAppByProject()
  endif()

endmacro()


# Get TribitsExProj2 stuff with find_package(<Package>) for each
# package/component independently.
#
macro(getTribitsExProj2StuffForAppByPackage)

  # Find each package and gather up all the <Package>::all_libs targets
  set(APP_DEPS_LIB_TARGETS "")
  foreach (packageName IN LISTS ${PROJECT_NAME}_USE_COMPONENTS)
    find_package(${packageName} CONFIG REQUIRED)
    message("Found ${packageName}!")
    list(APPEND APP_DEPS_LIB_TARGETS ${packageName}::all_libs)
  endforeach()
  print_var(APP_DEPS_LIB_TARGETS)

  # Set TribitsExProj2_SELECTED_PACKAGE_LIST
  set(TribitsExProj2_SELECTED_PACKAGE_LIST ${${PROJECT_NAME}_USE_COMPONENTS})
  # NOTE: We are setting his here since TribitsExProj2Config.cmake is not being
  # read in in this case.

  # Get compilers from first package listed
  list(GET ${PROJECT_NAME}_USE_COMPONENTS 0 firstPkg)
  setCompilersForAppFromConfigFileCompilers(${firstPkg})

endmacro()


# Get TribitsExProj2 stuff from find_package(TribitsExProj2)
#
macro(getTribitsExProj2StuffForAppByProject)

  find_package(TribitsExProj2 CONFIG REQUIRED
    COMPONENTS ${${PROJECT_NAME}_USE_COMPONENTS})

  message("\nFound TribitsExProj2!  Here are the details: ")
  message("   TribitsExProj2_DIR = ${TribitsExProj2_DIR}")
  message("   TribitsExProj2_VERSION = ${TribitsExProj2_VERSION}")
  message("   TribitsExProj2_PACKAGE_LIST = ${TribitsExProj2_PACKAGE_LIST}")
  message("   TribitsExProj2_TPL_LIST = ${TribitsExProj2_TPL_LIST}")
  message("   TribitsExProj2_BUILD_SHARED_LIBS = ${TribitsExProj2_BUILD_SHARED_LIBS}")
  message("End of TribitsExProj2 details\n")

  # Make sure to use same compilers and flags as TribitsExProj2
  setCompilersForAppFromConfigFileCompilers(TribitsExProj2)

  # Get the libraries for building and linking
  if (${PROJECT_NAME}_USE_COMPONENTS)
    set(APP_DEPS_LIB_TARGETS TribitsExProj2::all_selected_libs)
  else()
    set(APP_DEPS_LIB_TARGETS TribitsExProj2::all_libs)
  endif()

endmacro()


# Get compilers and compiler flags from the imported
# ``TribitsExProj2Config.cmake`` or ``<Package>Config.cmake`` file.
#
# Here ``prefix`` is the prefix for the variables read in from the
# *Config.cmake file.
#
macro(setCompilersForAppFromConfigFileCompilers prefix)

  message("-- Setting compilers and flags read in from '${prefix}Config.cmake' file:")

  set(CMAKE_CXX_COMPILER ${${prefix}_CXX_COMPILER} )
  set(CMAKE_C_COMPILER ${${prefix}_C_COMPILER} )
  set(CMAKE_Fortran_COMPILER ${${prefix}_Fortran_COMPILER} )

  set(CMAKE_CXX_FLAGS "${${prefix}_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "${${prefix}_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${${prefix}_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")

  cmake_print_variables(CMAKE_CXX_COMPILER)
  cmake_print_variables(CMAKE_C_COMPILER)
  cmake_print_variables(CMAKE_Fortran_COMPILER)
  cmake_print_variables(CMAKE_CXX_FLAGS)
  cmake_print_variables(CMAKE_C_FLAGS)
  cmake_print_variables(CMAKE_Fortran_FLAGS)

endmacro()


# Add compiler defines to the ``app`` target for optionally supported packages
# from upstream TribitExProj
#
function(addAppDepCompileDefines)
  addAppDepCompileDefine("Package1")
  addAppDepCompileDefine("Package2")
  addAppDepCompileDefine("Package3")
endfunction()


function(addAppDepCompileDefine componentName)
  if (${componentName} IN_LIST TribitsExProj2_SELECTED_PACKAGE_LIST)
    string(TOUPPER "${componentName}" componentNameUpper)
    target_compile_definitions(app PRIVATE TRIBITSEXAPP2_HAVE_${componentNameUpper})
  endif()
endfunction()


# Return the extended dependency string from the app at runtime given the
# enabled packages from TribitsExProj2.
#
function(getExpectedAppDepsStr expectedDepsStrOut)

  set(tpl1 "tpl1")  # ToDo: Change to Tpl1{no deps}?
  set(tpl2a "Tpl2a{${tpl1}}")
  set(tpl2b "Tpl2b{no deps}")
  set(tpl3 "Tpl3{${tpl2a}, ${tpl2b}}")
  set(tpl4 "Tpl4{${tpl3}, ${tpl2a}, ${tpl2b}}")

  set(package1 "Package1{${tpl1}}")

  if (Package2_ENABLE_Tpl3)
    set(package2_tpl3_deps_str ", ${tpl3}")
  else()
    set(package2_tpl3_deps_str "")
  endif()
  set(package2 "Package2{${package1}${package2_tpl3_deps_str}}")

  if (Package3_ENABLE_Package2)
    set(package3_package2_deps_str "${package2}, ")
  else()
    set(package3_package2_deps_str "")
  endif()
  if (Package3_ENABLE_Tpl4)
    set(package3_tpl4_deps_str "${tpl4}, ")
  else()
    set(package3_tpl4_deps_str "")
  endif()
  set(package3
    "Package3{${package3_package2_deps_str}${package1}, ${package3_tpl4_deps_str}${tpl2a}, ${tpl2b}}")

  if (TARGET Package3::all_libs)
    set(depsStr "${package3}")
  elseif (TARGET Package2::all_libs)
    set(depsStr "${package2}")
  elseif (TARGET Package1::all_libs)
    set(depsStr "${package1}")
  else()
    message("No package targets exist!")
    set(depsStr "no deps")
  endif()

  set(${expectedDepsStrOut} "${depsStr}" PARENT_SCOPE)

endfunction()


function(print_var varName)
  message("-- ${varName} = '${${varName}}'")
endfunction()
