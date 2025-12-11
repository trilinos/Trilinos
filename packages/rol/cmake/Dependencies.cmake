#SET(SUBPACKAGES_DIRS_CLASSIFICATION_OPTREQS
  # SubPackageName       Directory       Class       Req/Opt
#  xrol                   experimental    EX          OPTIONAL
#)

SET(LIB_REQUIRED_DEP_PACKAGES Teuchos)
SET(LIB_OPTIONAL_DEP_PACKAGES
    Belos
    Tpetra
    Thyra
    Sacado
    Intrepid2
    MiniTensor
    Shards
    Amesos2
    Ifpack2
    MueLu
    Tempus
)
IF(Trilinos_MAJOR_MINOR_VERSION LESS 160200)
  LIST(APPEND LIB_OPTIONAL_DEP_PACKAGES
       Epetra
       Intrepid
       Amesos
       TriKota
  )
ENDIF()
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS Boost ArrayFireCPU Eigen pebbl)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS gtest)

