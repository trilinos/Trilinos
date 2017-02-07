# Here we list Epetra as a required dependence so that the Thyra/Epetra
# adapters will get enabled by default.  Note that Belos and ML do *not* have
# a required dependence on Epetra but the Stratimikos Belos and ML adapters
# need the Thyra/Epetra adapters.
SET(LIB_REQUIRED_DEP_PACKAGES)
SET(LIB_OPTIONAL_DEP_PACKAGES Amesos AztecOO Belos Ifpack ML EpetraExt ThyraEpetraAdapters)
SET(TEST_REQUIRED_DEP_PACKAGES ThyraEpetraAdapters)
SET(TEST_OPTIONAL_DEP_PACKAGES Triutils Ifpack2 ThyraTpetraAdapters)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)

# Note: EpetraExt is used by the AztecOO adapters
