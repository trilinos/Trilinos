
INCLUDE(Multiline_Set)
INCLUDE(Generate_Script_File_Target)

FUNCTION(CREATE_DEPENCENCY_HANDLING_TEST_CASE TEST_NAME TEST_ARGUMENTS TEST_PASS_REGEX)

  MULTILINE_SET(TEST_CMND
    "#!/bin/bash\n"
    "if [ -e CMakeCache.txt ]\n"
    "  then rm CMakeCache.txt\n"
    "fi\n"
    "${CMAKE_COMMAND}"
    "  -C ${TRILINOS_BUILD_DIR}/CMakeCache.txt"
    "  -D Trilinos_PACKAGES_OVERRIDE:BOOL=ON"
    "  -D Trilinos_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING:BOOL=ON"
    "  -D Trilinos_ENABLE_CXX:BOOL=OFF"
    "  -D Trilinos_ENABLE_Fortran:BOOL=OFF"
    "   ${TEST_ARGUMENTS}"
    "  ${TRILINOS_HOME_DIR}\n"
    )
  
  SET(TEST_SCRIPT_FILE "${TEST_NAME}.sh")

  GENERATE_SCRIPT_FILE_TARGET(${TEST_SCRIPT_FILE} ${TEST_CMND}
    ${TRILINOS_HOME_DIR}/CMakeLists.txt
    ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
    )
  
  ADD_TEST( ${TEST_NAME} ${TEST_SCRIPT_FILE} )
  
  SET_TESTS_PROPERTIES( ${TEST_NAME}
    PROPERTIES PASS_REGULAR_EXPRESSION ${TEST_PASS_REGEX} )

ENDFUNCTION()


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableAllPackages
  "-D Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils EpetraExt Isorropia Thyra AztecOO Galeri Amesos Ifpack ML Belos RBGen Stratimikos 16"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_DisableEpetra_EnableAllPackages
  "-D Trilinos_ENABLE_Epetra:BOOL=OFF -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos Zoltan RTOp Thyra ML RBGen Stratimikos 7"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetra_DisableTeuchos_EnableAllForwardPackages
  "-D Trilinos_ENABLE_Epetra:BOOL=ON -D Trilinos_ENABLE_Teuchos:BOOL=OFF -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Epetra Triutils AztecOO ML 4"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_DisableEpetra_DisableTeuchos_EnableAllPackages
  "-D Trilinos_ENABLE_Epetra:BOOL=OFF -D Trilinos_ENABLE_Teuchos:BOOL=OFF -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Zoltan ML 2"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableTeuchos_EnableAllForwardPackages
  "-D Trilinos_ENABLE_Teuchos:BOOL=ON -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils EpetraExt Isorropia Thyra AztecOO Galeri Amesos Ifpack ML Belos RBGen Stratimikos 16"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetraExt_EnableAllForwardPackages
  "-D Trilinos_ENABLE_EpetraExt:BOOL=ON -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils EpetraExt Isorropia Thyra AztecOO Galeri Amesos Ifpack ML Belos Stratimikos 15"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetraExt
  "-D Trilinos_ENABLE_EpetraExt:BOOL=ON"
  "Final set of enabled packages:  Teuchos Epetra EpetraExt 3"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetraExt_EnableAllOptionalPackages
  "-D Trilinos_ENABLE_EpetraExt:BOOL=ON -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos Epetra Triutils EpetraExt 4"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetraExt_DisableAztecOO_EnableAllForwardPackages
  "-D Trilinos_ENABLE_EpetraExt:BOOL=ON -D Trilinos_ENABLE_AztecOO:BOOL=OFF -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils EpetraExt Isorropia Thyra Galeri Amesos Ifpack ML Belos Stratimikos 14"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableAmesos_DisableZoltan_EnableAllForwardPackages
  "-D Trilinos_ENABLE_Amesos:BOOL=ON -D Trilinos_ENABLE_Zoltan:BOOL=OFF -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos Epetra RTOp Triutils EpetraExt Thyra AztecOO Galeri Amesos Ifpack ML Belos Stratimikos 13"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetra_DisableThyra_EnableAllForwardPackages
  "-D Trilinos_ENABLE_Epetra:BOOL=ON -D Trilinos_ENABLE_Thyra:BOOL=OFF -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos Epetra Zoltan Triutils EpetraExt Isorropia AztecOO Galeri Amesos Ifpack ML Belos RBGen 13"
  )
