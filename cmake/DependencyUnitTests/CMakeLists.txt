
INCLUDE(PackageAddAdvancedTest)

INCLUDE(ParseVariableArguments)


SET(PACKAGE_NAME TrilinosDep)


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

  PACKAGE_ADD_ADVANCED_TEST( ${TEST_NAME}
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
    "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils Tpetra EpetraExt Isorropia Thyra AztecOO Galeri Amesos Ifpack ML Belos RBGen Stratimikos Sacado 18"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  DisableEpetra_EnableAllPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=OFF -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Zoltan RTOp Tpetra Thyra ML RBGen Stratimikos Sacado 9"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableEpetra_DisableTeuchos_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=ON -DTrilinos_ENABLE_Teuchos:BOOL=OFF -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Setting Trilinos_ENABLE_RTOp=OFF because RTOp has a required library dependence on disabled package Teuchos"
    "Final set of enabled packages:  Epetra Zoltan Triutils AztecOO ML Sacado 6"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  DisableEpetra_DisableTeuchos_EnableAllPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=OFF -DTrilinos_ENABLE_Teuchos:BOOL=OFF -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Zoltan ML Sacado 3"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableTeuchos_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Teuchos:BOOL=ON -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils Tpetra EpetraExt Isorropia Thyra AztecOO Galeri Amesos Ifpack ML Belos RBGen Stratimikos Sacado 18"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableEpetraExt_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils Tpetra EpetraExt Isorropia Thyra AztecOO Galeri Amesos Ifpack ML Belos RBGen Stratimikos Sacado 18"
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
    "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils Tpetra EpetraExt Isorropia Thyra Galeri Amesos Ifpack ML Belos RBGen Stratimikos Sacado 17"

  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableAmesos_DisableZoltan_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Amesos:BOOL=ON -DTrilinos_ENABLE_Zoltan:BOOL=OFF -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra RTOp Triutils Tpetra EpetraExt Thyra AztecOO Galeri Amesos Ifpack ML Belos Stratimikos 14"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableEpetra_DisableThyra_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=ON -DTrilinos_ENABLE_Thyra:BOOL=OFF -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan Triutils Tpetra EpetraExt Isorropia AztecOO Galeri Amesos Ifpack ML Belos RBGen Sacado 15"
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
    "Final set of enabled packages:  Teuchos Epetra RTOp EpetraExt Thyra 5"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableSecondaryStableCode_EnableAllOptionalPackges
  ARGS  -DTrilinos_ENABLE_ALL_PACKAGES=ON -DTrilinos_ENABLE_SECONDARY_STABLE_CODE=ON -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils Tpetra EpetraExt Isorropia Thyra AztecOO Galeri Amesos Ifpack ML Belos RBGen Stratimikos Sacado Phalanx 19"
    "Final set of enabled TPLs:  BLAS LAPACK Boost 3"
  )
