
SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  CoreLibs  src  PS  REQUIRED
  GoodStuff  src/good_stuff  SS  OPTIONAL
  CrazyStuff  src/crazy_stuff  EX  OPTIONAL
  Epetra  adapters/epetra  PS  OPTIONAL
  EpetraExt  adapters/epetraext  PS  OPTIONAL
  Tpetra  adapters/tpetra  PS  OPTIONAL
  )

# NOTE: The above subpackages automatically become required and optional LIB
# package dependencies (prefixed by 'Thyra') added to the variables below.
# There is no need to add them again!

SET(LIB_REQUIRED_DEP_PACKAGES)
SET(LIB_OPTIONAL_DEP_PACKAGES MissingPackage)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)

ALLOW_MISSING_EXTERNAL_PACKAGES(MissingPackage)

SET(REGRESSION_EMAIL_LIST thyra-boneheads@gmail.com)
