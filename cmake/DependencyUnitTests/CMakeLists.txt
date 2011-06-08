
INCLUDE(PackageAddAdvancedTest)

INCLUDE(ParseVariableArguments)

#
# Unit tests for just creating the depencency XML files
#

FUNCTION(CREATE_DEPENCENCY_XML_DUMP_TEST_CASE TEST_NAME)

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "ARGS;PASS_REGULAR_EXPRESSION_ALL;COMPARE_XML_FILE"
     #options
     ""
     ${ARGN}
     )

  SET(TEST_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/DepXmlDumpTests_${TEST_NAME})

  PACKAGE_ADD_ADVANCED_TEST( DepXmlDumpTests_${TEST_NAME}
    OVERALL_WORKING_DIRECTORY TEST_NAME
    TEST_0 CMND ${CMAKE_COMMAND}
      ARGS
        -DPROJECT_NAME=Trilinos
        -DTrilinos_CMAKE_DIR=${TRILINOS_HOME_DIR}/cmake
        -DTrilinos_DEPS_HOME_DIR:PATH=${TRILINOS_HOME_DIR}/cmake/DependencyUnitTests/MockTrilinos
        -DTrilinos_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR=${TEST_OUTPUT_DIR}
        ${PARSE_ARGS}
        -P ${TRILINOS_HOME_DIR}/cmake/package_arch/PackageArchDumpDepsXmlScript.cmake
      PASS_REGULAR_EXPRESSION_ALL ${PARSE_PASS_REGULAR_EXPRESSION_ALL}
    TEST_1 CMND ${CMAKE_COMMAND}
      ARGS
        -E compare_files
        ${TEST_OUTPUT_DIR}/TrilinosPackageDependencies.xml
        ${TRILINOS_HOME_DIR}/cmake/python/UnitTests/${PARSE_COMPARE_XML_FILE}
    )

ENDFUNCTION()


CREATE_DEPENCENCY_XML_DUMP_TEST_CASE(
  BasePackagesTpls
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_PACKAGES='22'"
    "Trilinos_NUM_TPLS='19'"
  COMPARE_XML_FILE
    TrilinosPackageDependencies.gold.xml
  )


CREATE_DEPENCENCY_XML_DUMP_TEST_CASE(
  WithPreCopyrightTrilinos
  ARGS
    -DTrilinos_EXTRA_REPOSITORIES=preCopyrightTrilinos
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_PACKAGES='24'"
    "Trilinos_NUM_TPLS='20'"
  COMPARE_XML_FILE
    TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml
  )


CREATE_DEPENCENCY_XML_DUMP_TEST_CASE(
  WithPreCopyrightTrilinosExtraTrilinosRepo
  ARGS
    -DTrilinos_EXTRA_REPOSITORIES=preCopyrightTrilinos,extraTrilinosRepo
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_PACKAGES='25'"
    "Trilinos_NUM_TPLS='22'"
  COMPARE_XML_FILE
    TrilinosPackageDependencies.preCopyrightTrilinos.extraTrilinosRepo.gold.xml
  )

#
# Unit tests for full dependency handling
#

FUNCTION(CREATE_DEPENCENCY_HANDLING_TEST_CASE TEST_NAME)

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "ARGS;PASS_REGULAR_EXPRESSION_ALL"
     #options
     ""
     ${ARGN}
     )

  PACKAGE_ADD_ADVANCED_TEST( DepTests_${TEST_NAME}
    OVERALL_WORKING_DIRECTORY TEST_NAME
    TEST_0 CMND ${CMAKE_COMMAND}
      ARGS
        -DTrilinos_DEPS_HOME_DIR:PATH=${TRILINOS_HOME_DIR}/cmake/DependencyUnitTests/MockTrilinos
        -DTrilinos_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING:BOOL=ON
        -DTrilinos_DEPS_XML_OUTPUT_FILE:STRING=
        -DTrilinos_DEPS_HTML_OUTPUT_FILE:STRING=
        -DTrilinos_ENABLE_C:BOOL=OFF
        -DTrilinos_ENABLE_CXX:BOOL=OFF
        -DTrilinos_ENABLE_Fortran:BOOL=OFF
        ${PARSE_ARGS}
        ${TRILINOS_HOME_DIR}
      PASS_REGULAR_EXPRESSION_ALL ${PARSE_PASS_REGULAR_EXPRESSION_ALL}
    )

ENDFUNCTION()


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableAllPackages
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 18"
    "WARNING: MissingPackage is being ignored since its directory is missing and MissingPackage_ALLOW_MISSING_EXTERNAL_PACKAGE = TRUE"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableTests
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON -DTrilinos_ENABLE_TESTS:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 18"
    "Setting Teuchos_ENABLE_TESTS=ON"
    "Setting RTOp_ENABLE_TESTS=ON"
    "Setting Epetra_ENABLE_TESTS=ON"
    "Setting ML_ENABLE_TESTS=ON"
    "Setting Stratimikos_ENABLE_TESTS=ON"
    "Setting RBGen_ENABLE_TESTS=ON"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  DisableEpetra_EnableAllPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=OFF -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Zoltan Shards Tpetra Thyra Intrepid ML Stratimikos RBGen 10"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableEpetra_DisableTeuchos_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=ON -DTrilinos_ENABLE_Teuchos:BOOL=OFF -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Setting Trilinos_ENABLE_RTOp=OFF because RTOp has a required library dependence on disabled package Teuchos"
    "Final set of enabled packages:  Epetra Zoltan Triutils AztecOO ML 5"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  DisableEpetra_DisableTeuchos_EnableAllPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=OFF -DTrilinos_ENABLE_Teuchos:BOOL=OFF -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Zoltan Shards ML 3"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableTeuchos_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Teuchos:BOOL=ON -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 18"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableEpetraExt_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 17"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableEpetraExt
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra EpetraExt 3"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableEpetraExt_EnableAllOptionalPackages
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan Triutils EpetraExt 5"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableEpetraExt_DisableAztecOO_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON -DTrilinos_ENABLE_AztecOO:BOOL=OFF
    -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: Not Setting Trilinos_ENABLE_Belos=ON even through Stratimikos_ENABLE_Belos=ON because Trilinos_ENABLE_SECONDARY_STABLE_CODE=OFF"
    "WARNING: Not Setting Trilinos_ENABLE_Sacado=ON even through Intrepid_ENABLE_Sacado=ON because Trilinos_ENABLE_SECONDARY_STABLE_CODE=OFF"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils EpetraExt Thyra Isorropia Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 16"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableAmesos_DisableZoltan_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Amesos:BOOL=ON -DTrilinos_ENABLE_Zoltan:BOOL=OFF
    -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Shards Triutils EpetraExt Thyra AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos 14"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableEpetra_DisableThyra_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=ON -DTrilinos_ENABLE_Thyra:BOOL=OFF
    -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan Shards Triutils Tpetra EpetraExt Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML RBGen 15"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableAllPackages_DefaultTPLs
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK DUMMY 3"
    "Final set of non-enabled TPLs:  MPI Boost Scotch METIS ParMETIS CppUnit ADOLC ADIC TVMET y12m SuperLUDist SuperLU UMFPACK AMD PETSC MUMPS 16"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableBoost
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON -DTPL_ENABLE_Boost:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK Boost DUMMY 4"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableBoostParMETIS
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON -DTPL_ENABLE_Boost:BOOL=ON
    -DTPL_ENABLE_ParMETIS:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK Boost ParMETIS DUMMY 5"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableZoltan_ZoltanEnableParMETIS
  ARGS  -DTrilinos_ENABLE_Zoltan=ON -DZoltan_ENABLE_ParMETIS=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Setting TPL_ENABLE_ParMETIS=ON because Zoltan_ENABLE_ParMETIS=ON"
    "Final set of enabled TPLs:  ParMETIS 1"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableZoltan_ZoltanEnableParMETIS_TplDisableParMETIS
  ARGS -DTrilinos_ENABLE_Zoltan=ON -DZoltan_ENABLE_ParMETIS=ON -DTPL_ENABLE_ParMETIS=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: Setting Zoltan_ENABLE_ParMETIS=OFF which was ON since TPL_ENABLE_ParMETIS=OFF"
    "Final set of enabled TPLs:  0"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableThyra_ThyraEnableEpetraExt
  ARGS  -DTrilinos_ENABLE_Thyra=ON -DThyra_ENABLE_EpetraExt=ON
    -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra EpetraExt Thyra 5"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableSecondaryStableCode_EnableAllOptionalPackges
  ARGS  -DTrilinos_ENABLE_ALL_PACKAGES=ON -DTrilinos_ENABLE_SECONDARY_STABLE_CODE=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx 21"
    "Final set of enabled TPLs:  BLAS LAPACK Boost DUMMY 4"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableStratimikos_EnableAllOptionalPackages_BadSystem1
  ARGS  -DTrilinos_ENABLE_Stratimikos=ON
    -DTrilinos_HOSTTYPE:STRING=BadSystem1
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly disabled packages on input .+:  Stokhos Ifpack ML 3"
    "Setting Stratimikos_ENABLE_Ifpack=OFF because Stratimikos has an optional library dependence on disabled package Ifpack"
    "Setting Stratimikos_ENABLE_ML=OFF because Stratimikos has an optional library dependence on disabled package ML"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra AztecOO Galeri Amesos Stratimikos 11"
    "Final set of non-enabled packages:  Shards Tpetra Stokhos Sacado Isorropia Intrepid Ifpack ML Belos RBGen Phalanx 11"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableStratimikos_EnableAllOptionalPackages_BadSystem2
  ARGS  -DTrilinos_ENABLE_Stratimikos=ON
    -DTrilinos_HOSTTYPE:STRING=BadSystem2
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly disabled packages on input .+:  Stokhos Ifpack 2"
    "Setting Stratimikos_ENABLE_Ifpack=OFF because Stratimikos has an optional library dependence on disabled package Ifpack"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra Isorropia AztecOO Galeri Amesos ML Stratimikos 13"
    "Final set of non-enabled packages:  Shards Tpetra Stokhos Sacado Intrepid Ifpack Belos RBGen Phalanx 9"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableStratimikos_EnableIfpack_EnableAllOptionalPackages_BadSystem2
  ARGS  -DTrilinos_ENABLE_Stratimikos=ON -DTrilinos_ENABLE_Ifpack=ON
    -DTrilinos_HOSTTYPE:STRING=BadSystem2
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: User has set Trilinos_ENABLE_Ifpack=ON but the" 
    "package Ifpack is not supported on this platform type 'BadSystem2'"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos 14"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  preCopyrightTrilinos_EnableTeko
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
    -DTrilinos_ENABLE_Teko=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_TPLS='20'"
    "Explicitly enabled packages on input .* Teko 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting Trilinos_ENABLE_Teuchos=ON because Teko_ENABLE_Teuchos=ON"
    "Setting Trilinos_ENABLE_Epetra=ON because Teko_ENABLE_Epetra=ON"
    "Trilinos_ENABLE_Thyra=ON because Teko_ENABLE_Thyra=ON"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos Teko 15"
    "Final set of non-enabled packages:  Shards Tpetra Stokhos Sacado Intrepid Belos RBGen Phalanx Stalix 9"
    "Final set of enabled TPLs:  BLAS LAPACK TekoDepTPL 3"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  preCopyrightTrilinos_EnableStalix
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
    -DTrilinos_ENABLE_Stalix=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .* Stalix 1"
    "Explicitly disabled packages on input .* Stokhos 1"
    "WARNING: Not Setting Trilinos_ENABLE_Sacado=ON even through Stalix_ENABLE_Sacado=ON because Trilinos_ENABLE_SECONDARY_STABLE_CODE=OFF"
    "Setting Trilinos_ENABLE_Tpetra=ON because Stalix_ENABLE_Tpetra=ON"
    "Setting Trilinos_ENABLE_Thyra=ON because Stalix_ENABLE_Thyra=ON"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Sacado Thyra Stalix 10"
    "Final set of non-enabled packages:  Shards Stokhos Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Teko 14"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  preCopyrightTrilinos_ShowInvalidPackageNameError
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
    -DSHOW_INVALID_PACKAGE_NAME_ERROR=ON
    -DTrilinos_ASSERT_MISSING_PACKAGES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Error, the package 'Anasazi' is listed as a dependency of the package"
    "'Teko' in the list 'Trilinos_PACKAGES' but the package 'Anasazi' is either"
    "not defined or is listed later in the package order.  Check the spelling of"
    "'Anasazi' or see how it is listed in"
    "Trilinos_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS in relationship to 'Teko'"
    "Configuring incomplete, errors occurred"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  extraTrilinosRepo_EnableExtraPack
  ARGS "-DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos,extraTrilinosRepo"
    -DTrilinos_ENABLE_ExtraPack=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_TPLS='22'"
    "Explicitly enabled packages on input .* ExtraPack 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting Trilinos_ENABLE_Teuchos=ON because ExtraPack_ENABLE_Teuchos=ON"
    "Setting Trilinos_ENABLE_Teko=ON because Trilinos_ENABLE_ExtraPack=ON"
    "Final set of enabled packages:  Teuchos RTOp Epetra Triutils EpetraExt Thyra AztecOO Amesos Ifpack ML Stratimikos Teko ExtraPack 13"
    "Final set of non-enabled packages:  Zoltan Shards Tpetra Stokhos Sacado Isorropia Galeri Intrepid Belos RBGen Phalanx Stalix 12"
    "Final set of enabled TPLs:  BLAS LAPACK TekoDepTPL ExtraTPL2 4"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  extraTrilinosRepo_EnableExtraPack_EnableExtraTPL1
  ARGS "-DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos,extraTrilinosRepo"
    -DTrilinos_ENABLE_ExtraPack=ON -DTPL_ENABLE_ExtraTPL1=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_TPLS='22'"
    "Explicitly enabled packages on input .* ExtraPack 1"
    "Final set of enabled TPLs:  BLAS LAPACK TekoDepTPL ExtraTPL1 ExtraTPL2 5"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_CI_EnableTeko
  ARGS -DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=Continuous
    -DTrilinos_ENABLE_Teko=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .* Teko 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting Trilinos_ENABLE_Teuchos=ON because Teko_ENABLE_Teuchos=ON"
    "Setting Trilinos_ENABLE_Epetra=ON because Teko_ENABLE_Epetra=ON"
    "Trilinos_ENABLE_Thyra=ON because Teko_ENABLE_Thyra=ON"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos Teko 15"
    "Final set of non-enabled packages:  Shards Tpetra Stokhos Sacado Intrepid Belos RBGen Phalanx Stalix 9"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_Nightly_EnableExtraPack
  ARGS  -DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=Nightly
    -DTrilinos_ENABLE_ExtraPack=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .* ExtraPack 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting Trilinos_ENABLE_Teuchos=ON because ExtraPack_ENABLE_Teuchos=ON"
    "Setting Trilinos_ENABLE_Teko=ON because Trilinos_ENABLE_ExtraPack=ON"
    "Final set of enabled packages:  Teuchos RTOp Epetra Triutils EpetraExt Thyra AztecOO Amesos Ifpack ML Stratimikos Teko ExtraPack 13"
    "Final set of non-enabled packages:  Zoltan Shards Tpetra Stokhos Sacado Isorropia Galeri Intrepid Belos RBGen Phalanx Stalix 12"
  )
