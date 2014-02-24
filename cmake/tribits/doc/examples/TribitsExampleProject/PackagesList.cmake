

INCLUDE(TribitsListHelpers)

#
# Define the ${PROJECT_NAME} package names, directories, and classification.
#
# Package classifications are:
#
#   PS: Primary Stable Package
#
#     Primary Stable Packages have at least some Primary Stable Code which is
#     expected to be fully tested before every push to the global repo.  The
#     default enable for PS packages is empty "" which allows the PS package
#     to be enabled implicitly based on other criteria.  The option
#     ${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON will cause all PS packages to be enabled
#     unless they are explicitly disabled.
#
#   SS: Secondary Stable Package
#
#     Secondary Stable Packages have no PS code or they would be classified as
#     PS packages.  A package must be classified as SS if it has a required
#     dependency on another SS package or SS TPL.  A package may also be
#     declared SS to avoid requiring it to be tested before every push to the
#     global repo.  For example, a package that does not provide any
#     significant functionally like Didasko is classified as a SS package even
#     through it could be classified as PS just based on its required package
#     and TPL dependencies.  SS packages will have their default enables set
#     to empty "".  This allows them to be enabled implicilty.  When
#     ${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON but
#     ${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE=OFF, the SS packages will not be
#     enabled.  However, when ${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON and
#     ${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE=ON, then SS packages will be
#     enabled if they are not explicitly disabled.  Packages that are SS but
#     not PS must be disabled in pre-push testing.  However, SS packages are
#     tested by the post-push CI and nightly testing processes.
#
#   EX: Experimental Package
#
#     Experimental packages are those packages that contain no PS or SS
#     code. The default enable for EX packages is always OFF which requires
#     that they be explicitly enabled in order to be turned on. EX packages
#     must be disabled in pre-push testing and are not tested as part of the
#     post-push CI or nightly testing processes.  However, package developers
#     of EX pacakges are encouraged to set up their own nightly testing for
#     thier EX packages.
#
# NOTE: These packages must be listed in strictly assending order in terms of
# package dependencies.  If you get the order wrong, then an error message
# will be printed during configuration with CMake.
#

SET( ${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
  SimpleCxx               packages/simple_cxx                 PS
  MixedLanguage           packages/mixed_language             PS
  PackageWithSubpackages  packages/package_with_subpackages   PS
  WrapExternal            packages/wrap_external              EX
  )
