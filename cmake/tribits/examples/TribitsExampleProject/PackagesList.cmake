TRIBITS_REPOSITORY_DEFINE_PACKAGES(
  SimpleCxx          packages/simple_cxx         PT
  MixedLang          packages/mixed_lang         PT
  ExternalPkg        ExternalPkg                 ST
  WithSubpackages    packages/with_subpackages   PT
  WrapExternal       packages/wrap_external      ST
  )

TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(WrapExternal Windows)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(ExternalPkg)
