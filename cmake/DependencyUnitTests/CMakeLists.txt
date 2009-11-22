
INCLUDE(PackageAddAdvancedTest)

INCLUDE(ParseVariableArguments)


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
       -DTrilinos_PACKAGES_FILE:STRING=DependencyUnitTests/TrilinosPackagesForUnitTests
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
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra EpetraExt 3"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableEpetraExt_EnableAllOptionalPackages
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan Triutils EpetraExt 5"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableEpetraExt_DisableAztecOO_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON -DTrilinos_ENABLE_AztecOO:BOOL=OFF -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: Not Setting Trilinos_ENABLE_Belos=ON even through Stratimikos_ENABLE_Belos=ON because Trilinos_ENABLE_SECONDARY_STABLE_CODE=OFF"
    "WARNING: Not Setting Trilinos_ENABLE_Sacado=ON even through Intrepid_ENABLE_Sacado=ON because Trilinos_ENABLE_SECONDARY_STABLE_CODE=OFF"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils EpetraExt Thyra Isorropia Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 16"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableAmesos_DisableZoltan_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Amesos:BOOL=ON -DTrilinos_ENABLE_Zoltan:BOOL=OFF -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Shards Triutils EpetraExt Thyra AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos 14"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableEpetra_DisableThyra_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=ON -DTrilinos_ENABLE_Thyra:BOOL=OFF -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan Shards Triutils Tpetra EpetraExt Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML RBGen 15"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableAllPackages_DefaultTPLs
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK 2"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableBoost
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON -DTPL_ENABLE_Boost:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK Boost 3"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableBoostParMETIS
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON -DTPL_ENABLE_Boost:BOOL=ON -DTPL_ENABLE_ParMETIS:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK Boost ParMETIS 4"
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
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra EpetraExt Thyra 5"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableSecondaryStableCode_EnableAllOptionalPackges
  ARGS  -DTrilinos_ENABLE_ALL_PACKAGES=ON -DTrilinos_ENABLE_SECONDARY_STABLE_CODE=ON
    -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx 21"
    "Final set of enabled TPLs:  BLAS LAPACK Boost 3"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableStratimikos_EnableAllOptionalPackages_BadSystem1
  ARGS  -DTrilinos_ENABLE_Stratimikos=ON -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON
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
  ARGS  -DTrilinos_ENABLE_Stratimikos=ON -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON
    -DTrilinos_HOSTTYPE:STRING=BadSystem2
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly disabled packages on input .+:  Stokhos Ifpack 2"
    "Setting Stratimikos_ENABLE_Ifpack=OFF because Stratimikos has an optional library dependence on disabled package Ifpack"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra Isorropia AztecOO Galeri Amesos ML Stratimikos 13"
    "Final set of non-enabled packages:  Shards Tpetra Stokhos Sacado Intrepid Ifpack Belos RBGen Phalanx 9"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableStratimikos_EnableIfpack_EnableAllOptionalPackages_BadSystem2
  ARGS  -DTrilinos_ENABLE_Stratimikos=ON -DTrilinos_ENABLE_Ifpack=ON -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON
    -DTrilinos_HOSTTYPE:STRING=BadSystem2
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: User has set Trilinos_ENABLE_Ifpack=ON but the" 
    "package Ifpack is not supported on this platform type 'BadSystem2'"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos 14"
  )
