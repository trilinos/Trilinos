
INCLUDE(TribitsAddAdvancedTest)

INCLUDE(ParseVariableArguments)

SET(MOCK_TRILINOS_DIR ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_MOCK_TRILINOS_DIR})


#
# Unit tests for just creating the depencency XML files
#

FUNCTION(CREATE_DEPENDENCY_XML_DUMP_TEST_CASE TEST_NAME)

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

  TRIBITS_ADD_ADVANCED_TEST( DepXmlDumpTests_${TEST_NAME}
    OVERALL_WORKING_DIRECTORY TEST_NAME
    TEST_0 CMND ${CMAKE_COMMAND}
      ARGS
        -DPROJECT_NAME=${PROJECT_NAME}
        -D${PROJECT_NAME}_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
        -D${PROJECT_NAME}_DEPS_HOME_DIR:PATH=${MOCK_TRILINOS_DIR}
        -D${PROJECT_NAME}_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR=${TEST_OUTPUT_DIR}
        ${PARSE_ARGS}
        -P ${${PROJECT_NAME}_TRIBITS_DIR}/package_arch/TribitsDumpDepsXmlScript.cmake
      PASS_REGULAR_EXPRESSION_ALL ${PARSE_PASS_REGULAR_EXPRESSION_ALL}
    TEST_1 CMND ${CMAKE_COMMAND}
      ARGS
        -E compare_files
        ${TEST_OUTPUT_DIR}/${PROJECT_NAME}PackageDependencies.xml
        ${${PROJECT_NAME}_TRIBITS_DIR}/python/UnitTests/${PARSE_COMPARE_XML_FILE}
    )

ENDFUNCTION()


CREATE_DEPENDENCY_XML_DUMP_TEST_CASE(
  BasePackagesTpls
  PASS_REGULAR_EXPRESSION_ALL
    "${PROJECT_NAME}_NUM_PACKAGES='23'"
    "${PROJECT_NAME}_NUM_TPLS='19'"
  COMPARE_XML_FILE
    ${PROJECT_NAME}PackageDependencies.gold.xml
  )


CREATE_DEPENDENCY_XML_DUMP_TEST_CASE(
  WithPreCopyright${PROJECT_NAME}
  ARGS
    -D${PROJECT_NAME}_EXTRA_REPOSITORIES=preCopyright${PROJECT_NAME}
  PASS_REGULAR_EXPRESSION_ALL
    "${PROJECT_NAME}_NUM_PACKAGES='25'"
    "${PROJECT_NAME}_NUM_TPLS='20'"
  COMPARE_XML_FILE
    ${PROJECT_NAME}PackageDependencies.preCopyright${PROJECT_NAME}.gold.xml
  )


CREATE_DEPENDENCY_XML_DUMP_TEST_CASE(
  WithPreCopyright${PROJECT_NAME}Extra${PROJECT_NAME}Repo
  ARGS
    -D${PROJECT_NAME}_EXTRA_REPOSITORIES=preCopyright${PROJECT_NAME},extra${PROJECT_NAME}Repo
  PASS_REGULAR_EXPRESSION_ALL
    "${PROJECT_NAME}_NUM_PACKAGES='26'"
    "${PROJECT_NAME}_NUM_TPLS='22'"
  COMPARE_XML_FILE
    ${PROJECT_NAME}PackageDependencies.preCopyright${PROJECT_NAME}.extra${PROJECT_NAME}Repo.gold.xml
  )

#
# Unit tests for full dependency handling
#

FUNCTION(CREATE_DEPENDENCY_HANDLING_TEST_CASE TEST_NAME)

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "ARGS;PASS_REGULAR_EXPRESSION_ALL"
     #options
     ""
     ${ARGN}
     )

  TRIBITS_ADD_ADVANCED_TEST( DepTests_${TEST_NAME}
    OVERALL_WORKING_DIRECTORY TEST_NAME
    TEST_0 CMND ${CMAKE_COMMAND}
      ARGS
        -D${PROJECT_NAME}_DEPS_HOME_DIR:PATH=${MOCK_TRILINOS_DIR}
        -D${PROJECT_NAME}_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING:BOOL=ON
        -D${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE:STRING=
        -D${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE:STRING=
        -D${PROJECT_NAME}_ENABLE_C:BOOL=OFF
        -D${PROJECT_NAME}_ENABLE_CXX:BOOL=OFF
        -D${PROJECT_NAME}_ENABLE_Fortran:BOOL=OFF
        ${PARSE_ARGS}
        ${TRILINOS_HOME_DIR}
      PASS_REGULAR_EXPRESSION_ALL ${PARSE_PASS_REGULAR_EXPRESSION_ALL}
    )

ENDFUNCTION()


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages
  ARGS -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  ${PROJECT_NAME}Framework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 19"
    "WARNING: MissingPackage is being ignored since its directory is missing and MissingPackage_ALLOW_MISSING_EXTERNAL_PACKAGE = TRUE"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableTests
  ARGS -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=ON -D${PROJECT_NAME}_ENABLE_TESTS:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  ${PROJECT_NAME}Framework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 19"
    "Setting Teuchos_ENABLE_TESTS=ON"
    "Setting RTOp_ENABLE_TESTS=ON"
    "Setting Epetra_ENABLE_TESTS=ON"
    "Setting ML_ENABLE_TESTS=ON"
    "Setting Stratimikos_ENABLE_TESTS=ON"
    "Setting RBGen_ENABLE_TESTS=ON"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  DisableEpetra_EnableAllPackages
  ARGS -D${PROJECT_NAME}_ENABLE_Epetra:BOOL=OFF -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  ${PROJECT_NAME}Framework Teuchos RTOp Zoltan Shards Tpetra Thyra Intrepid ML RBGen 10"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetra_DisableTeuchos_EnableAllForwardPackages
  ARGS -D${PROJECT_NAME}_ENABLE_Epetra:BOOL=ON -D${PROJECT_NAME}_ENABLE_Teuchos:BOOL=OFF
    -D${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Setting ${PROJECT_NAME}_ENABLE_RTOp=OFF because RTOp has a required library dependence on disabled package Teuchos"
    "Final set of enabled packages:  Epetra Zoltan Triutils AztecOO ML 5"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  DisableEpetra_DisableTeuchos_EnableAllPackages
  ARGS -D${PROJECT_NAME}_ENABLE_Epetra:BOOL=OFF -D${PROJECT_NAME}_ENABLE_Teuchos:BOOL=OFF -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  ${PROJECT_NAME}Framework Zoltan Shards ML 4"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableTeuchos_EnableAllForwardPackages
  ARGS -D${PROJECT_NAME}_ENABLE_Teuchos:BOOL=ON -D${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 18"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetraExt_EnableAllForwardPackages
  ARGS -D${PROJECT_NAME}_ENABLE_EpetraExt:BOOL=ON -D${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 18"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetraExt
  ARGS -D${PROJECT_NAME}_ENABLE_EpetraExt:BOOL=ON -D${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra EpetraExt 3"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetraExt_EnableAllOptionalPackages
  ARGS -D${PROJECT_NAME}_ENABLE_EpetraExt:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan Triutils EpetraExt 5"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetraExt_DisableAztecOO_EnableAllForwardPackages
  ARGS -D${PROJECT_NAME}_ENABLE_EpetraExt:BOOL=ON -D${PROJECT_NAME}_ENABLE_AztecOO:BOOL=OFF
    -D${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: Not Setting ${PROJECT_NAME}_ENABLE_Belos=ON even though Stratimikos_ENABLE_Belos=ON because ${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE=OFF"
    "WARNING: Not Setting ${PROJECT_NAME}_ENABLE_Sacado=ON even though Intrepid_ENABLE_Sacado=ON because ${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE=OFF"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 17"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAmesos_DisableZoltan_EnableAllForwardPackages
  ARGS -D${PROJECT_NAME}_ENABLE_Amesos:BOOL=ON -D${PROJECT_NAME}_ENABLE_Zoltan:BOOL=OFF
    -D${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Shards Triutils EpetraExt Thyra AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos 14"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyraCoreLibs_EnableTests
  ARGS -D${PROJECT_NAME}_ENABLE_ThyraCoreLibs:BOOL=ON -D${PROJECT_NAME}_ENABLE_TESTS:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Setting ${PROJECT_NAME}_ENABLE_Thyra=ON because ${PROJECT_NAME}_ENABLE_ThyraCoreLibs=ON"
    "Setting Thyra_ENABLE_TESTS=ON because ThyraCoreLibs_ENABLE_TESTS=ON"
    "Setting Thyra_ENABLE_EXAMPLES=ON because ThyraCoreLibs_ENABLE_EXAMPLES=ON"
    "Final set of enabled packages:  Teuchos RTOp Thyra 3"
    "Final set of enabled SE packages:  Teuchos RTOp ThyraCoreLibs Thyra 4"    
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetra_DisableThyraCoreLibs_EnableAllForwardPackages
  ARGS -D${PROJECT_NAME}_ENABLE_Epetra:BOOL=ON -D${PROJECT_NAME}_ENABLE_ThyraCoreLibs:BOOL=OFF
    -D${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan Shards Triutils Tpetra EpetraExt Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML RBGen 15"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_DefaultTPLs
  ARGS -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK DUMMY 3"
    "Final set of non-enabled TPLs:  MPI Boost Scotch METIS ParMETIS CppUnit ADOLC ADIC TVMET y12m SuperLUDist SuperLU UMFPACK AMD PETSC MUMPS 16"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableBoost
  ARGS -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=ON -DTPL_ENABLE_Boost:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK Boost DUMMY 4"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableBoostParMETIS
  ARGS -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=ON -DTPL_ENABLE_Boost:BOOL=ON
    -DTPL_ENABLE_ParMETIS:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK Boost ParMETIS DUMMY 5"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableZoltan_ZoltanEnableParMETIS
  ARGS  -D${PROJECT_NAME}_ENABLE_Zoltan=ON -DZoltan_ENABLE_ParMETIS=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Setting TPL_ENABLE_ParMETIS=ON because Zoltan_ENABLE_ParMETIS=ON"
    "Final set of enabled TPLs:  ParMETIS 1"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableZoltan_ZoltanEnableParMETIS_TplDisableParMETIS
  ARGS -D${PROJECT_NAME}_ENABLE_Zoltan=ON -DZoltan_ENABLE_ParMETIS=ON -DTPL_ENABLE_ParMETIS=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: Setting Zoltan_ENABLE_ParMETIS=OFF which was ON since TPL_ENABLE_ParMETIS=OFF"
    "Final set of enabled TPLs:  0"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyra_ThyraEnableEpetraExt
  ARGS  -D${PROJECT_NAME}_ENABLE_Thyra=ON -D${PROJECT_NAME}_ENABLE_ThyraEpetraExt=ON
    -D${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra EpetraExt Thyra 5"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableTpetra_EnableAllForwardPackages
  ARGS  -D${PROJECT_NAME}_ENABLE_Tpetra=ON
    -D${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Thyra 8"
    "Final set of enabled SE packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt ThyraCoreLibs ThyraEpetra ThyraEpetraExt ThyraTpetra Thyra 12"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableSecondaryStableCode_EnableAllOptionalPackges
  ARGS  -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON -D${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  ${PROJECT_NAME}Framework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx 22"
    "Final set of enabled TPLs:  BLAS LAPACK Boost DUMMY 4"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableStratimikos_EnableAllOptionalPackages_BadSystem1
  ARGS  -D${PROJECT_NAME}_ENABLE_Stratimikos=ON
    -D${PROJECT_NAME}_HOSTTYPE:STRING=BadSystem1
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly disabled packages on input .+:  Stokhos Ifpack ML 3"
    "Setting Stratimikos_ENABLE_Ifpack=OFF because Stratimikos has an optional library dependence on disabled package Ifpack"
    "Setting Stratimikos_ENABLE_ML=OFF because Stratimikos has an optional library dependence on disabled package ML"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra AztecOO Galeri Amesos Stratimikos 11"
    "Final set of non-enabled packages:  ${PROJECT_NAME}Framework Shards Tpetra Stokhos Sacado Isorropia Intrepid Ifpack ML Belos RBGen Phalanx 12"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableStratimikos_EnableAllOptionalPackages_BadSystem2
  ARGS  -D${PROJECT_NAME}_ENABLE_Stratimikos=ON
    -D${PROJECT_NAME}_HOSTTYPE:STRING=BadSystem2
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly disabled packages on input .+:  Stokhos Ifpack 2"
    "Setting Stratimikos_ENABLE_Ifpack=OFF because Stratimikos has an optional library dependence on disabled package Ifpack"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra Isorropia AztecOO Galeri Amesos ML Stratimikos 13"
    "Final set of enabled SE packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt ThyraCoreLibs ThyraEpetra ThyraEpetraExt Thyra Isorropia AztecOO Galeri Amesos ML Stratimikos 16"
    "Final set of non-enabled packages:  ${PROJECT_NAME}Framework Shards Tpetra Stokhos Sacado Intrepid Ifpack Belos RBGen Phalanx 10"
    "Final set of non-enabled SE packages:  ${PROJECT_NAME}Framework Shards Tpetra Stokhos Sacado ThyraGoodStuff ThyraCrazyStuff ThyraTpetra Intrepid Ifpack Belos RBGen Phalanx 13"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableStratimikos_EnableIfpack_EnableAllOptionalPackages_BadSystem2
  ARGS  -D${PROJECT_NAME}_ENABLE_Stratimikos=ON -D${PROJECT_NAME}_ENABLE_Ifpack=ON
    -D${PROJECT_NAME}_HOSTTYPE:STRING=BadSystem2
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: User has set ${PROJECT_NAME}_ENABLE_Ifpack=ON but the" 
    "package Ifpack is not supported on this platform type 'BadSystem2'"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos 14"
    "Final set of enabled SE packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt ThyraCoreLibs ThyraEpetra ThyraEpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos 17"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  preCopyright${PROJECT_NAME}_EnableTeko
  ARGS -D${PROJECT_NAME}_EXTRA_REPOSITORIES:STRING=preCopyright${PROJECT_NAME}
    -D${PROJECT_NAME}_ENABLE_Teko=ON
  PASS_REGULAR_EXPRESSION_ALL
    "${PROJECT_NAME}_NUM_TPLS='20'"
    "Explicitly enabled packages on input .* Teko 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting ${PROJECT_NAME}_ENABLE_Teuchos=ON because Teko_ENABLE_Teuchos=ON"
    "Setting ${PROJECT_NAME}_ENABLE_Epetra=ON because Teko_ENABLE_Epetra=ON"
    "${PROJECT_NAME}_ENABLE_Thyra=ON because Teko_ENABLE_Thyra=ON"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos Teko 16"
    "Final set of non-enabled packages:  ${PROJECT_NAME}Framework Shards Stokhos Sacado Intrepid Belos RBGen Phalanx Stalix 9"
    "Final set of enabled TPLs:  BLAS LAPACK TekoDepTPL 3"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  preCopyright${PROJECT_NAME}_EnableStalix
  ARGS -D${PROJECT_NAME}_EXTRA_REPOSITORIES:STRING=preCopyright${PROJECT_NAME}
    -D${PROJECT_NAME}_ENABLE_Stalix=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .* Stalix 1"
    "Explicitly disabled packages on input .* Stokhos 1"
    "WARNING: Not Setting ${PROJECT_NAME}_ENABLE_Sacado=ON even though Stalix_ENABLE_Sacado=ON because ${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE=OFF"
    "Setting ${PROJECT_NAME}_ENABLE_Tpetra=ON because Stalix_ENABLE_Tpetra=ON"
    "Setting ${PROJECT_NAME}_ENABLE_Thyra=ON because Stalix_ENABLE_Thyra=ON"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Sacado Thyra Stalix 10"
    "Final set of non-enabled packages:  ${PROJECT_NAME}Framework Shards Stokhos Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Teko 15"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  preCopyright${PROJECT_NAME}_ShowInvalidPackageNameError
  ARGS -D${PROJECT_NAME}_EXTRA_REPOSITORIES:STRING=preCopyright${PROJECT_NAME}
    -DSHOW_INVALID_PACKAGE_NAME_ERROR=ON
    -D${PROJECT_NAME}_ASSERT_MISSING_PACKAGES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Error, the package 'Anasazi' is listed as a dependency of the package"
    "'Teko' is in the list '${PROJECT_NAME}_SE_PACKAGES' but the package 'Anasazi' is"
    "either not defined or is listed later in the package order.  Check the"
    "pelling of 'Anasazi' or see how it is listed in"
    "${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS in relationship to 'Teko'"
    "Configuring incomplete, errors occurred"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  extra${PROJECT_NAME}Repo_EnableExtraPack
  ARGS "-D${PROJECT_NAME}_EXTRA_REPOSITORIES:STRING=preCopyright${PROJECT_NAME},extra${PROJECT_NAME}Repo"
    -D${PROJECT_NAME}_ENABLE_ExtraPack=ON
  PASS_REGULAR_EXPRESSION_ALL
    "${PROJECT_NAME}_NUM_TPLS='22'"
    "Explicitly enabled packages on input .* ExtraPack 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting ${PROJECT_NAME}_ENABLE_Teuchos=ON because ExtraPack_ENABLE_Teuchos=ON"
    "Setting ${PROJECT_NAME}_ENABLE_Teko=ON because ${PROJECT_NAME}_ENABLE_ExtraPack=ON"
    "Final set of enabled packages:  Teuchos RTOp Epetra Triutils EpetraExt Thyra AztecOO Amesos Ifpack ML Stratimikos Teko ExtraPack 13"
    "Final set of non-enabled packages:  ${PROJECT_NAME}Framework Zoltan Shards Tpetra Stokhos Sacado Isorropia Galeri Intrepid Belos RBGen Phalanx Stalix 13"
    "Final set of enabled TPLs:  BLAS LAPACK TekoDepTPL ExtraTPL2 4"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  extra${PROJECT_NAME}Repo_EnableExtraPack_EnableExtraTPL1
  ARGS "-D${PROJECT_NAME}_EXTRA_REPOSITORIES:STRING=preCopyright${PROJECT_NAME},extra${PROJECT_NAME}Repo"
    -D${PROJECT_NAME}_ENABLE_ExtraPack=ON -DTPL_ENABLE_ExtraTPL1=ON
  PASS_REGULAR_EXPRESSION_ALL
    "${PROJECT_NAME}_NUM_TPLS='22'"
    "Explicitly enabled packages on input .* ExtraPack 1"
    "Final set of enabled TPLs:  BLAS LAPACK TekoDepTPL ExtraTPL1 ExtraTPL2 5"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_CI_EnableTeko
  ARGS -D${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=Continuous
    -D${PROJECT_NAME}_ENABLE_Teko=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .* Teko 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting ${PROJECT_NAME}_ENABLE_Teuchos=ON because Teko_ENABLE_Teuchos=ON"
    "Setting ${PROJECT_NAME}_ENABLE_Epetra=ON because Teko_ENABLE_Epetra=ON"
    "${PROJECT_NAME}_ENABLE_Thyra=ON because Teko_ENABLE_Thyra=ON"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos Teko 16"
    "Final set of non-enabled packages:  ${PROJECT_NAME}Framework Shards Stokhos Sacado Intrepid Belos RBGen Phalanx Stalix 9"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_Nightly_EnableExtraPack
  ARGS  -D${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=Nightly
    -D${PROJECT_NAME}_ENABLE_ExtraPack=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .* ExtraPack 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting ${PROJECT_NAME}_ENABLE_Teuchos=ON because ExtraPack_ENABLE_Teuchos=ON"
    "Setting ${PROJECT_NAME}_ENABLE_Teko=ON because ${PROJECT_NAME}_ENABLE_ExtraPack=ON"
    "Final set of enabled packages:  Teuchos RTOp Epetra Triutils EpetraExt Thyra AztecOO Amesos Ifpack ML Stratimikos Teko ExtraPack 13"
    "Final set of non-enabled packages:  ${PROJECT_NAME}Framework Zoltan Shards Tpetra Stokhos Sacado Isorropia Galeri Intrepid Belos RBGen Phalanx Stalix 13"
  )
