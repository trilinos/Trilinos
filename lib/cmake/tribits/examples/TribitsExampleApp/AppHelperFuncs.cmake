include(CMakePrintHelpers)


# Find TribitsExProj package(s), load compilers and compiler options, and get
# CMake lib targets (which have include dirs also) to link against.
#
# On return, sets the vars in the current scope:
#
# * TribitsExProj_SELECTED_PACKAGE_LIST: List of all of the packages pulled in
# * from the TribtsExProj.
#
# * APP_DEPS_LIB_TARGETS: List of all of the IMPORTED CMake targets that 'app'
#   must link against
#
# * CMAKE_<LANG>_COMPILER and CMAKE_<LANG>_FLAGS pulled in from the
#   TribitsExProjConfig.config or a <Package>Config.cmake file
#
macro(getTribitsExProjStuffForApp)

  set(${PROJECT_NAME}_FIND_INDIVIDUAL_PACKAGES OFF CACHE BOOL
    "Set to TRUE to find individual packages and OFF to find project TribitsExProj")

  if (${PROJECT_NAME}_FIND_INDIVIDUAL_PACKAGES)
    getTribitsExProjStuffForAppByPackage()
  else()
    getTribitsExProjStuffForAppByProject()
  endif()

endmacro()


# Get TribitsExProj stuff with find_package(<Package>) for each
# package/component independently.
#
macro(getTribitsExProjStuffForAppByPackage)

  # Find each package and gather up all the <Package>::all_libs targets
  set(APP_DEPS_LIB_TARGETS "")
  foreach (packageName IN LISTS ${PROJECT_NAME}_USE_COMPONENTS)
    find_package(${packageName} REQUIRED)
    message("Found ${packageName}!")
    list(APPEND APP_DEPS_LIB_TARGETS ${packageName}::all_libs)
  endforeach()
  print_var(APP_DEPS_LIB_TARGETS)

  # Set TribitsExProj_SELECTED_PACKAGE_LIST
  set(TribitsExProj_SELECTED_PACKAGE_LIST ${${PROJECT_NAME}_USE_COMPONENTS})
  # NOTE: We are setting his here since TribitsExProjConfig.cmake is not being
  # read in in this case.

  # Get compilers from first package listed
  list(GET ${PROJECT_NAME}_USE_COMPONENTS 0 firstPkg)
  setCompilersForAppFromConfigFileCompilers(${firstPkg})

endmacro()


# Get TribitsExProj stuff from find_package(TribitsExProj)
#
macro(getTribitsExProjStuffForAppByProject)

  find_package(TribitsExProj REQUIRED COMPONENTS ${${PROJECT_NAME}_USE_COMPONENTS})

  message("\nFound TribitsExProj!  Here are the details: ")
  message("   TribitsExProj_DIR = ${TribitsExProj_DIR}")
  message("   TribitsExProj_VERSION = ${TribitsExProj_VERSION}")
  message("   TribitsExProj_PACKAGE_LIST = ${TribitsExProj_PACKAGE_LIST}")
  message("   TribitsExProj_BUILD_SHARED_LIBS = ${TribitsExProj_BUILD_SHARED_LIBS}")
  message("End of TribitsExProj details\n")

  # Make sure to use same compilers and flags as TribitsExProj
  setCompilersForAppFromConfigFileCompilers(TribitsExProj)

  # Get the libraries for building and linking
  if (${PROJECT_NAME}_USE_COMPONENTS)
    set(APP_DEPS_LIB_TARGETS TribitsExProj::all_selected_libs)
  else()
    set(APP_DEPS_LIB_TARGETS TribitsExProj::all_libs)
  endif()

endmacro()


# Get compilers and compiler flags from the imported
# ``TribitsExProjConfig.cmake`` or ``<Package>Config.cmake`` file.
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
  addAppDepCompileDefine("SimpleCxx")
  addAppDepCompileDefine("MixedLang")
  addAppDepCompileDefine("WithSubpackagesA")
  addAppDepCompileDefine("WithSubpackagesB")
  addAppDepCompileDefine("WithSubpackagesC")
endfunction()


function(addAppDepCompileDefine componentName)
  if (TARGET ${componentName}::all_libs)
    string(TOUPPER "${componentName}" componentNameUpper)
    target_compile_definitions(app PRIVATE TRIBITSEXAPP_HAVE_${componentNameUpper})
  endif()
endfunction()
# NOTE: Above, we look to see if the 'all_libs' target for a package is
# defined as a way to know if that package is enabled.  That will determine if
# the package is enabled even if it is just implicitly enabled and therefore
# is not listed in the list of COMPONENTS passed to find_package().


# Return the extended dependency string from the app at runtime given the
# enabled packages from TribitsExProj.
#
function(getExpectedAppDepsStr expectedDepsStrOut)

  if (TARGET SimpleTpl::all_libs)
    set(simpleCxxDeps "simpletpl ")
  else()
    set(simpleCxxDeps "")
  endif()
  set(simpleCxxDeps "${simpleCxxDeps}headeronlytpl")

  set(withSubpackagesADeps "SimpleCxx ${simpleCxxDeps}")
  set(withSubpackagesBDeps "A ${withSubpackagesADeps} SimpleCxx ${simpleCxxDeps}")
  set(withSubpackagesCDeps "B ${withSubpackagesBDeps} A ${withSubpackagesADeps}")

  set(depsStr "")
  appendExpectedAppDepsStr("WithSubpackagesC"
    "WithSubpackagesC:${withSubpackagesCDeps}"
    depsStr)
  appendExpectedAppDepsStr("WithSubpackagesB"
    "WithSubpackagesB:${withSubpackagesBDeps}"
    depsStr)
  appendExpectedAppDepsStr("WithSubpackagesA"
    "WithSubpackagesA:${withSubpackagesADeps}"
    depsStr)
  appendExpectedAppDepsStr("MixedLang" "MixedLang:Mixed Language" depsStr)
  appendExpectedAppDepsStr("SimpleCxx" "SimpleCxx:${simpleCxxDeps}" depsStr)

  set(${expectedDepsStrOut} "${depsStr}" PARENT_SCOPE)

endfunction()


function(appendExpectedAppDepsStr componentName str depsStrOut)
  set(depsStr "${${depsStrOut}}")  # Should be value of var in parent scope!
  #message("-- depsStr (inner) = '${depsStr}'")
  if (TARGET ${componentName}::all_libs)
    if (depsStr)
      set(depsStr "${depsStr}[;] ${str}")
    else()
      set(depsStr "${str}")
    endif()
  endif()
  set(${depsStrOut} "${depsStr}" PARENT_SCOPE)
endfunction()


function(print_var varName)
  message("-- ${varName} = '${${varName}}'")
endfunction()
