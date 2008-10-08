
INCLUDE(Multiline_Set)
INCLUDE(Generate_Script_File_Target)


FUNCTION(CREATE_DEPENCENCY_HANDLING_TEST_CASE TEST_NAME TEST_ARGUMENTS TEST_PASS_REGEX)

  MULTILINE_SET(TEST_CMND
    "#!/bin/bash\n"
    "if [ -e CMakeCache.txt ]\n"
    "  then rm CMakeCache.txt\n"
    "fi\n"
    "${CMAKE_COMMAND}"
    "  -D Trilinos_PACKAGES_OVERRIDE:BOOL=ON"
    "  -D Trilinos_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING:BOOL=ON"
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
  "Final set of enabled packages:  Teuchos RTOp Epetra Triutils EpetraExt Thyra RBGen 7"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_DisableEpetra_EnableAllPackages
  "-D Trilinos_ENABLE_Epetra:BOOL=OFF -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos RTOp Thyra RBGen 4"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetra_DisableTeuchos_EnableAllForwardPackages
  "-D Trilinos_ENABLE_Epetra:BOOL=ON -D Trilinos_ENABLE_Teuchos:BOOL=OFF -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Epetra Triutils 2"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableTeuchos_EnableAllForwardPackages
  "-D Trilinos_ENABLE_Teuchos:BOOL=ON -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos RTOp Epetra Triutils EpetraExt Thyra RBGen 7"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetraExt_EnableAllForwardPackages
  "-D Trilinos_ENABLE_EpetraExt:BOOL=ON -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos RTOp Epetra Triutils EpetraExt Thyra 6"
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
