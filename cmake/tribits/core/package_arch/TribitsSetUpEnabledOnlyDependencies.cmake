# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(TribitsDeprecatedHelpers)


# Function that creates enable-only dependency data-structures
#
# For each enabled package `<Package>`, this function sets up the global list
# var::
#
#   <Package>_FULL_ENABLED_DEP_PACKAGES
#
# If ``${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_PACKAGES`` is
# set, then ``<Package>_FULL_ENABLED_DEP_PACKAGES`` will only be sets for
# those packages.  Otherwise, ``<Package>_FULL_ENABLED_DEP_PACKAGES`` will be
# set for all packages listed in `${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES`_.
#
# NOTE: The modern TriBITS implementation does not need this full list of
# dependencies for each package.  Only the function
# `tribits_find_most_recent_file_timestamp()` needs this.  (Therefore, this
# could not be striped out of TriBITS because there are still some projects
# that use this function.)
#
function(tribits_set_up_enabled_only_dependencies)

  set(GENERATE_EXPORT_DEPENDENCIES ${${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES})
  if ("${${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_PACKAGES}" STREQUAL ""
      AND NOT
      "${${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_SE_PACKAGES}" STREQUAL ""
    )
    tribits_deprecated(
      "The cache var"
      " ${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_SE_PACKAGES"
      "='${${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_SE_PACKAGES}'"
      " is deprecated!  Please instead set"
      " ${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_PACKAGES"
      "='${${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_SE_PACKAGES}'")
    set(${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_PACKAGES
      ${${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_SE_PACKAGES} )
  endif()

  # Determine lastExportTribitsPackage if not to generate any of these full
  # dependency lists
  set(lastExportTribitsPackage "")
  if (GENERATE_EXPORT_DEPENDENCIES
      AND ${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_PACKAGES
    )
    # Find the last enabled package for which an export file is requested.
    set(LAST_PKG_IDX -1)
    set(LAST_PKG)
    foreach(tribitsPkg  IN LISTS
        ${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_PACKAGES
      )
      set(PKG_IDX ${${tribitsPkg}_PKG_IDX})
      if (NOT "${PKG_IDX}" STREQUAL "")
        # The listed package is enabled so we will consider it
        if (PKG_IDX GREATER ${LAST_PKG_IDX})
          set(LAST_PKG_IDX ${PKG_IDX})
          set(LAST_PKG ${tribitsPkg})
        endif()
      endif()
    endforeach()
    if (LAST_PKG)
      # At least one listed package was enabled
      set(lastExportTribitsPackage ${LAST_PKG})
    else()
      # None of the listed packages were enabled so don't bother generating
      # any export dependencies
      set(GENERATE_EXPORT_DEPENDENCIES FALSE)
    endif()

  endif()

  if (GENERATE_EXPORT_DEPENDENCIES)

    if (lastExportTribitsPackage)
      message("\nSetting up export dependencies up through ${lastExportTribitsPackage} ...\n")
    else()
      message("\nSetting up export dependencies for all enabled packages ...\n")
    endif()

    foreach(tribitsPackage  IN LISTS  ${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES)
      tribits_package_set_full_enabled_dep_packages(${tribitsPackage})
      if (${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES)
        set(printedVar FALSE)
        print_nonempty_var_with_spaces(${tribitsPackage}_FULL_ENABLED_DEP_PACKAGES
          printedVar)
        if (NOT printedVar)
          message("-- ${tribitsPackage}: No library dependencies!")
        endif()
      endif()
      if ("${lastExportTribitsPackage}" STREQUAL "${tribitsPackage}")
        break()
      endif()
    endforeach()

  endif()

endfunction()


# Function that sets up the full package dependencies for the given *internal*
# enabled package ``${packageName}``, including all of its indirect upstream
# *internal* package dependencies.
#
# After running, this function sets the internal cache var:
#
#  * ``${packageName}_FULL_ENABLED_DEP_PACKAGES``
#
# NOTE: This function must be called for all of the upstream internal packages
# before calling it for this package.
#
# NOTE: The complexity of this function is O(<numPackages>^2) due to the
# sorting algorithm.  That is why it would be good to get rid of this function
# at some point (or refactor it to have a better complexity).
#
function(tribits_package_set_full_enabled_dep_packages  packageName)

  tribits_package_build_unsorted_full_enabled_dep_packages(${packageName}
    packageFullDepsList)

  tribits_package_sort_full_enabled_dep_packages(packageFullDepsList
    orderedPackageFullDepsList)

  global_set(${packageName}_FULL_ENABLED_DEP_PACKAGES ${orderedPackageFullDepsList})

endfunction()


# Helper function that builds the full list of internal upstream dep packages
# (with no duplicates) for a given internal package.
#
function(tribits_package_build_unsorted_full_enabled_dep_packages  packageName
    packageFullDepsListOut
  )

  set(packageFullDepsList "")
  foreach(depPkg  IN LISTS  ${packageName}_LIB_DEFINED_DEPENDENCIES)
    if ((${depPkg}_PACKAGE_BUILD_STATUS STREQUAL "INTERNAL")
        AND ((${packageName}_LIB_DEP_REQUIRED AND ${PROJECT_NAME}_ENABLE_${depPkg})
          OR ((NOT ${packageName}_LIB_DEP_REQUIRED) AND ${packageName}_ENABLE_${depPkg}))
      )
      list(APPEND  packageFullDepsList  ${depPkg})
    endif()
  endforeach()

  if(packageFullDepsList)
    list(REMOVE_DUPLICATES  packageFullDepsList)

    foreach(DEP_PACKAGE  IN LISTS  packageFullDepsList)
      list(APPEND packageFullDepsList  ${${DEP_PACKAGE}_FULL_ENABLED_DEP_PACKAGES})
    endforeach()

    list(REMOVE_DUPLICATES packageFullDepsList)
  endif()

  set(${packageFullDepsListOut} ${packageFullDepsList} PARENT_SCOPE)

endfunction()


# Helper function to sort the full set of upstream dep packages for a given
# internal package.
#
function(tribits_package_sort_full_enabled_dep_packages  packageFullDepsListName
    orderedPackageFullDepsListOut
  )

  set(orderedPackageFullDepsList "")

  foreach(depPkg  IN LISTS  ${packageFullDepsListName})

    set(depPkgIdx ${${depPkg}_PKG_IDX})

    set(sortedIndex 0)
    set(insertedDepPkg FALSE)

    foreach(sortedPackage  IN LISTS  orderedPackageFullDepsList)

      set(sortedPackageIdx ${${sortedPackage}_PKG_IDX})

      if (${depPkgIdx} GREATER ${sortedPackageIdx})
        list(INSERT  orderedPackageFullDepsList  ${sortedIndex}  ${depPkg})
        set(insertedDepPkg TRUE)
        break()
      endif()

      math(EXPR sortedIndex ${sortedIndex}+1)

    endforeach()

    if(NOT insertedDepPkg)
      list(APPEND  orderedPackageFullDepsList  ${depPkg})
    endif()

  endforeach()

  set(${orderedPackageFullDepsListOut} ${orderedPackageFullDepsList} PARENT_SCOPE)

endfunction()
