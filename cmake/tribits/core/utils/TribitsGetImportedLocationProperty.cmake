# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include_guard()

include(CMakeBuildTypesList)


# Return the value of the target property IMPORTED_LOCATION_<CONFIG> for every
# known <CONFIG> and returns the fist non-null value.  Otherwise, returns
# empty.
#
function(tribits_get_imported_location_property  targetName  importedLocationOut)

  set(importedLocation "")
  foreach (cmakeBuildType IN LISTS CMAKE_CONFIG_FILE_BUILD_TYPES_LIST)
    string(TOUPPER ${cmakeBuildType} cmakeBuildTypeUpper)
    get_target_property(importedLoc ${targetName}
      IMPORTED_LOCATION_${cmakeBuildTypeUpper})
    if (importedLoc)
      set(importedLocation "${importedLoc}")
      break()
    endif()
  endforeach()
  set(${importedLocationOut} ${importedLocation} PARENT_SCOPE)
endfunction()
