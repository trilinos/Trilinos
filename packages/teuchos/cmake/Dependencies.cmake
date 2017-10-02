TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    Core          core          PS  REQUIRED
    Parser        parser        PS  REQUIRED
    ParameterList parameterlist PS  REQUIRED
    Comm          comm          PS  REQUIRED
    Numerics      numerics      PS  REQUIRED
    Remainder     remainder     PS  REQUIRED
    KokkosCompat  kokkoscompat  PS  OPTIONAL
    KokkosComm    kokkoscomm    PS  OPTIONAL
  )
