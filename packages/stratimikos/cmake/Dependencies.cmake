# Here we list Epetra as a required dependence so that the Thyra/Epetra
# adapters will get enabled by default.  Note that Belos and ML do *not* have
# a required dependence on Epetra but the Stratimikos Belos and ML adpaters
# need the Thyra/Epetra adapters.
SET(LIB_REQUIRED_DEP_PACKAGES ThyraEpetraAdapters)
SET(LIB_OPTIONAL_DEP_PACKAGES Amesos AztecOO Belos Ifpack ML)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES Triutils EpetraExt)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)
