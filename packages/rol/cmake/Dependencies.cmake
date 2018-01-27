#SET(SUBPACKAGES_DIRS_CLASSIFICATION_OPTREQS
  # SubPackageName       Directory       Class       Req/Opt
#  xrol                   experimental    EX          OPTIONAL
#)      

SET(LIB_REQUIRED_DEP_PACKAGES Teuchos)
SET(LIB_OPTIONAL_DEP_PACKAGES Belos Epetra Tpetra Thyra Sacado Intrepid MiniTensor Shards Amesos Amesos2 Ifpack2 MueLu TriKota Tempus)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES Gtest)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS Boost ArrayFireCPU Eigen)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)

