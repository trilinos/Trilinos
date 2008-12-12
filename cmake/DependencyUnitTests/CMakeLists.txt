
INCLUDE(Multiline_Set)
INCLUDE(Generate_Script_File_Target)

FUNCTION(CREATE_DEPENCENCY_HANDLING_TEST_CASE TEST_NAME TEST_ARGUMENTS TEST_PASS_REGEX)

  MULTILINE_SET(TEST_CMND
    "#!/bin/bash\n"
    "if [ ! -e ${TEST_NAME} ]\n"
    "  then mkdir ${TEST_NAME}\n"
    "fi\n"
    "cd ${TEST_NAME}\n"
    "if [ -e CMakeCache.txt ]\n"
    "  then rm CMakeCache.txt\n"
    "fi\n"
    "${CMAKE_COMMAND}"
    "  -D Trilinos_PACKAGES_FILE:STRING=DependencyUnitTests/TrilinosPackagesForUnitTests"
    "  -D Trilinos_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING:BOOL=ON"
    "  -D Trilinos_DEPS_XML_OUTPUT_FILE:STRING="
    "  -D Trilinos_DEPS_HTML_OUTPUT_FILE:STRING="
    "  -D Trilinos_ENABLE_C:BOOL=OFF"
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
  "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils Tpetra EpetraExt Isorropia Thyra AztecOO Galeri Amesos Ifpack ML Belos RBGen Stratimikos 17"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_DisableEpetra_EnableAllPackages
  "-D Trilinos_ENABLE_Epetra:BOOL=OFF -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos Zoltan RTOp Tpetra Thyra ML RBGen Stratimikos 8"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetra_DisableTeuchos_EnableAllForwardPackages
  "-D Trilinos_ENABLE_Epetra:BOOL=ON -D Trilinos_ENABLE_Teuchos:BOOL=OFF -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Epetra Zoltan Triutils AztecOO ML 5"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_DisableEpetra_DisableTeuchos_EnableAllPackages
  "-D Trilinos_ENABLE_Epetra:BOOL=OFF -D Trilinos_ENABLE_Teuchos:BOOL=OFF -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Zoltan ML 2"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableTeuchos_EnableAllForwardPackages
  "-D Trilinos_ENABLE_Teuchos:BOOL=ON -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils Tpetra EpetraExt Isorropia Thyra AztecOO Galeri Amesos Ifpack ML Belos RBGen Stratimikos 17"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetraExt_EnableAllForwardPackages
  "-D Trilinos_ENABLE_EpetraExt:BOOL=ON -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils Tpetra EpetraExt Isorropia Thyra AztecOO Galeri Amesos Ifpack ML Belos Stratimikos 16"
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
  "Final set of enabled packages:  Teuchos Epetra Zoltan RTOp Triutils Tpetra EpetraExt Isorropia Thyra Galeri Amesos Ifpack ML Belos Stratimikos 15"

  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableAmesos_DisableZoltan_EnableAllForwardPackages
  "-D Trilinos_ENABLE_Amesos:BOOL=ON -D Trilinos_ENABLE_Zoltan:BOOL=OFF -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos Epetra RTOp Triutils Tpetra EpetraExt Thyra AztecOO Galeri Amesos Ifpack ML Belos Stratimikos 14"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableEpetra_DisableThyra_EnableAllForwardPackages
  "-D Trilinos_ENABLE_Epetra:BOOL=ON -D Trilinos_ENABLE_Thyra:BOOL=OFF -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON"
  "Final set of enabled packages:  Teuchos Epetra Zoltan Triutils Tpetra EpetraExt Isorropia AztecOO Galeri Amesos Ifpack ML Belos RBGen 14"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableAllPackages_DefaultTPLs
  "-D Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON"
  "Final set of enabled TPLs:  BLAS LAPACK 2"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableAllPackages_EnableBoost
  "-D Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON -D TPL_ENABLE_Boost:BOOL=ON"
  "Final set of enabled TPLs:  BLAS LAPACK Boost 3"
  )


CREATE_DEPENCENCY_HANDLING_TEST_CASE(
  TrilinosDep_EnableAllPackages_EnableBoostMPI
  "-D Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON -D TPL_ENABLE_Boost:BOOL=ON  -D TPL_ENABLE_MPI:BOOL=ON"
  "Final set of enabled TPLs:  BLAS LAPACK MPI Boost 4"
  )
