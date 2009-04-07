
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
  TrilinosDep_EnableAllPackages
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils Tpetra EpetraExt Isorropia Thyra AztecOO Galeri Amesos Ifpack ML Belos RBGen Stratimikos 17"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_DisableEpetra_EnableAllPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=OFF -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Zoltan RTOp Tpetra Thyra ML RBGen Stratimikos 8"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetra_DisableTeuchos_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=ON -DTrilinos_ENABLE_Teuchos:BOOL=OFF -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Setting Trilinos_ENABLE_RTOp=OFF because RTOp has a required library dependence on disabled package Teuchos"
    "Final set of enabled packages:  Epetra Zoltan Triutils AztecOO ML 5"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_DisableEpetra_DisableTeuchos_EnableAllPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=OFF -DTrilinos_ENABLE_Teuchos:BOOL=OFF -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Zoltan ML 2"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableTeuchos_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Teuchos:BOOL=ON -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils Tpetra EpetraExt Isorropia Thyra AztecOO Galeri Amesos Ifpack ML Belos RBGen Stratimikos 17"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetraExt_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils Tpetra EpetraExt Isorropia Thyra AztecOO Galeri Amesos Ifpack ML Belos RBGen Stratimikos 17"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetraExt
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra EpetraExt 3"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetraExt_EnableAllOptionalPackages
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan Triutils EpetraExt 5"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetraExt_DisableAztecOO_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON -DTrilinos_ENABLE_AztecOO:BOOL=OFF -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils Tpetra EpetraExt Isorropia Thyra Galeri Amesos Ifpack ML Belos RBGen Stratimikos 16"

  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableAmesos_DisableZoltan_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Amesos:BOOL=ON -DTrilinos_ENABLE_Zoltan:BOOL=OFF -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra RTOp Triutils Tpetra EpetraExt Thyra AztecOO Galeri Amesos Ifpack ML Belos Stratimikos 14"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetra_DisableThyra_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=ON -DTrilinos_ENABLE_Thyra:BOOL=OFF -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan Triutils Tpetra EpetraExt Isorropia AztecOO Galeri Amesos Ifpack ML Belos RBGen 14"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableAllPackages_DefaultTPLs
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK 2"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableAllPackages_EnableBoost
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON -DTPL_ENABLE_Boost:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK Boost 3"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableAllPackages_EnableBoostParMETIS
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON -DTPL_ENABLE_Boost:BOOL=ON -DTPL_ENABLE_ParMETIS:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK Boost ParMETIS 4"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableZoltan_ZoltanEnableParMETIS
  ARGS  -DTrilinos_ENABLE_Zoltan=ON -DZoltan_ENABLE_ParMETIS=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  ParMETIS 1"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableThyra_ThyraEnableEpetraExt
  ARGS  -DTrilinos_ENABLE_Thyra=ON -DThyra_ENABLE_EpetraExt=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra RTOp EpetraExt Thyra 5"
  )
