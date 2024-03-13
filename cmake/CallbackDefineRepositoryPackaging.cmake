MACRO(TRIBITS_REPOSITORY_DEFINE_PACKAGING)
 
  # We need to make sure that these excludes only apply to Trilinos, not the global
  # project.
  SET(Trilinos_SOURCE_EXCLUDE_DIR ${Trilinos_SOURCE_DIR})

  SET(CPACK_SOURCE_IGNORE_FILES
    ${CPACK_SOURCE_IGNORE_FILES}
    ${Trilinos_SOURCE_EXCLUDE_DIR}/.*[.]pyc$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/.*classicMakefile$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/sparse_checkout.sh$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/cmake/tribits/common_tools/git/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/cmake/CMakeKitwareBacklog.txt$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/cmake/TODO$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/cmake/ctest/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/ITAPS/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/external/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/jpetra/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/cmmlib/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/configure.ac$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/configure$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/Makefile.am$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/Makefile.in$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/doc/[^b]
    ${Trilinos_SOURCE_EXCLUDE_DIR}/README_old$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/sampleScripts/old_autotools/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/sampleScripts/git-profiles/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/SIERRA/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/test/coverage/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/test/harness/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/test/utilities/README$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/test/utilities/dependencies/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/test/utilities/packages/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/test/utilities/r.*
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/scripts/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/release/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/common/DoxyfilePackageTemplate
    ${Trilinos_SOURCE_EXCLUDE_DIR}/stamp-h.in$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/configure.ac$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/aclocal.m4$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/configure$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/Makefile.am$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/Makefile.in$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/bootstrap$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/config/
    )

  if (NOT ${CMAKE_PROJECT_NAME}_ENABLE_TrilinosBuildStats)
    # Don't strip out all of commonTools/build_stats/, just the TriBITS
    # package-related files.  We want to keep these, even for a tarball
    # release.
    set(buildStatsPkgDir "${Trilinos_SOURCE_DIR}/commonTools/build_stats")
    list(REMOVE_ITEM CPACK_SOURCE_IGNORE_FILES "${buildStatsPkgDir}/")
    list(APPEND CPACK_SOURCE_IGNORE_FILES
      "${buildStatsPkgDir}/cmake/Dependencies.cmake"
      "${buildStatsPkgDir}/CMakeLists.txt" )
  endif()

  APPEND_SET(TRIBITS_CPACK_PACKAGES_TO_NOT_IGNORE TriBITS)
  
ENDMACRO()
