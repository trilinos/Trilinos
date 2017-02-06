# Here we list Epetra as a required dependence so that the Thyra/Epetra
# adapters will get enabled by default.  Note that Belos and ML do *not* have
# a required dependence on Epetra but the Stratimikos Belos and ML adapters
# need the Thyra/Epetra adapters.
TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_OPTIONAL_PACKAGES Amesos AztecOO Belos Ifpack ML EpetraExt ThyraEpetraAdapters
  TEST_REQUIRED_PACKAGES ThyraEpetraAdapters
  TEST_OPTIONAL_PACKAGES Triutils Ifpack2 ThyraTpetraAdapters
  )
# Note: EpetraExt is used by the AztecOO adapters
