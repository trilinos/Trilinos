TRIBITS_REPOSITORY_DEFINE_PACKAGES(
  SimpleCxx               packages/simple_cxx                 PT
  MixedLanguage           packages/mixed_language             PT
  PackageWithSubpackages  packages/package_with_subpackages   PT
  WrapExternal            packages/wrap_external              ST
  )

TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(WrapExternal Windows)
