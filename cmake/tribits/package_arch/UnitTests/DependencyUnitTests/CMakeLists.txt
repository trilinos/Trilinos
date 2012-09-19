
INCLUDE(TribitsAddAdvancedTest)

INCLUDE(ParseVariableArguments)

# The mock test project 
SET(MOCK_PROJECT_NAME Trilinos)

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
        -DPROJECT_NAME=${MOCK_PROJECT_NAME}
        -DPROJECT_SOURCE_DIR:PATH=${MOCK_TRILINOS_DIR}
        -D${MOCK_PROJECT_NAME}_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
        -D${MOCK_PROJECT_NAME}_NATIVE_REPOSITORIES=.
        -D${MOCK_PROJECT_NAME}_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR=${TEST_OUTPUT_DIR}
        ${PARSE_ARGS}
        -P ${${PROJECT_NAME}_TRIBITS_DIR}/package_arch/TribitsDumpDepsXmlScript.cmake
      PASS_REGULAR_EXPRESSION_ALL ${PARSE_PASS_REGULAR_EXPRESSION_ALL}
    TEST_1 CMND ${CMAKE_COMMAND}
      ARGS
        -E compare_files
        ${TEST_OUTPUT_DIR}/TrilinosPackageDependencies.xml
        ${${PROJECT_NAME}_TRIBITS_DIR}/python/UnitTests/${PARSE_COMPARE_XML_FILE}
    )

ENDFUNCTION()


CREATE_DEPENDENCY_XML_DUMP_TEST_CASE(
  BasePackagesTpls
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_PACKAGES='24'"
    "Trilinos_NUM_TPLS='19'"
  COMPARE_XML_FILE
    TrilinosPackageDependencies.gold.xml
  )


CREATE_DEPENDENCY_XML_DUMP_TEST_CASE(
  WithPreCopyrightTrilinos
  ARGS
    -DTrilinos_EXTRA_REPOSITORIES=preCopyrightTrilinos
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_PACKAGES='26'"
    "Trilinos_NUM_TPLS='20'"
  COMPARE_XML_FILE
    TrilinosPackageDependencies.preCopyrightTrilinos.gold.xml
  )


CREATE_DEPENDENCY_XML_DUMP_TEST_CASE(
  WithPreCopyrightTrilinosExtraTrilinosRepo
  ARGS
    -DTrilinos_EXTRA_REPOSITORIES=preCopyrightTrilinos,extraTrilinosRepo
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_PACKAGES='26'"
    "Trilinos_NUM_TPLS='22'"
  COMPARE_XML_FILE
    TrilinosPackageDependencies.preCopyrightTrilinos.extraTrilinosRepo.gold.xml
  )


CREATE_DEPENDENCY_XML_DUMP_TEST_CASE(
  WithExtraRepoOnePackage
  ARGS
    -DTrilinos_EXTRA_REPOSITORIES=extraRepoOnePackage
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_PACKAGES='24'"
    "Trilinos_NUM_TPLS='19'"
  COMPARE_XML_FILE
    TrilinosPackageDependencies.extraRepoOnePackage.gold.xml
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
        -D${MOCK_PROJECT_NAME}_TRIBITS_DIR:FILEPATH=${${PROJECT_NAME}_TRIBITS_DIR}
        -D${MOCK_PROJECT_NAME}_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING:BOOL=ON
        -D${MOCK_PROJECT_NAME}_DEPS_XML_OUTPUT_FILE:STRING=
        -D${MOCK_PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE:STRING=
        -D${MOCK_PROJECT_NAME}_ENABLE_C:BOOL=OFF
        -D${MOCK_PROJECT_NAME}_ENABLE_CXX:BOOL=OFF
        -D${MOCK_PROJECT_NAME}_ENABLE_Fortran:BOOL=OFF
        ${PARSE_ARGS}
        ${MOCK_TRILINOS_DIR}
      PASS_REGULAR_EXPRESSION_ALL ${PARSE_PASS_REGULAR_EXPRESSION_ALL}
    )

ENDFUNCTION()


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 19"
    "WARNING: MissingPackage is being ignored since its directory is missing and MissingPackage_ALLOW_MISSING_EXTERNAL_PACKAGE = TRUE"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableTests
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON -DTrilinos_ENABLE_TESTS:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 19"
    "Setting Teuchos_ENABLE_TESTS=ON"
    "Setting RTOp_ENABLE_TESTS=ON"
    "Setting Epetra_ENABLE_TESTS=ON"
    "Setting ML_ENABLE_TESTS=ON"
    "Setting Stratimikos_ENABLE_TESTS=ON"
    "Setting RBGen_ENABLE_TESTS=ON"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  DisableEpetra_EnableAllPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=OFF -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Zoltan Shards Tpetra Thyra Intrepid ML RBGen 10"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetra_DisableTeuchos_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=ON -DTrilinos_ENABLE_Teuchos:BOOL=OFF
    -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Setting Trilinos_ENABLE_RTOp=OFF because RTOp has a required library dependence on disabled package Teuchos"
    "Final set of enabled packages:  Epetra Zoltan Triutils AztecOO ML 5"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  DisableEpetra_DisableTeuchos_EnableAllPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=OFF -DTrilinos_ENABLE_Teuchos:BOOL=OFF -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  TrilinosFramework Zoltan Shards ML 4"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableTeuchos_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Teuchos:BOOL=ON -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 18"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetraExt_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 18"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetraExt
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra EpetraExt 3"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetraExt_EnableAllOptionalPackages
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan Triutils EpetraExt 5"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetraExt_DisableAztecOO_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_EpetraExt:BOOL=ON -DTrilinos_ENABLE_AztecOO:BOOL=OFF
    -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: Not Setting Trilinos_ENABLE_Belos=ON even though Stratimikos has an optional dependence on Belos because Trilinos_ENABLE_SECONDARY_STABLE_CODE=OFF"
    "WARNING: Not Setting Trilinos_ENABLE_Sacado=ON even though Intrepid has an optional dependence on Sacado because Trilinos_ENABLE_SECONDARY_STABLE_CODE=OFF"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Thyra Isorropia Galeri Amesos Intrepid Ifpack ML Stratimikos RBGen 17"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAmesos_DisableZoltan_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Amesos:BOOL=ON -DTrilinos_ENABLE_Zoltan:BOOL=OFF
    -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Shards Triutils EpetraExt Thyra AztecOO Galeri Amesos Intrepid Ifpack ML Stratimikos 14"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyraCoreLibs_EnableTests
  ARGS -DTrilinos_ENABLE_ThyraCoreLibs:BOOL=ON -DTrilinos_ENABLE_TESTS:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Setting Trilinos_ENABLE_Thyra=ON because Trilinos_ENABLE_ThyraCoreLibs=ON"
    "Setting Thyra_ENABLE_TESTS=ON because ThyraCoreLibs_ENABLE_TESTS=ON"
    "Setting Thyra_ENABLE_EXAMPLES=ON because ThyraCoreLibs_ENABLE_EXAMPLES=ON"
    "Final set of enabled packages:  Teuchos RTOp Thyra 3"
    "Final set of enabled SE packages:  Teuchos RTOp ThyraCoreLibs Thyra 4"    
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableStratimikos_DisableThyra
  ARGS -DTrilinos_ENABLE_Stratimikos=ON -DTrilinos_ENABLE_Thyra=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "Setting subpackage enable Trilinos_ENABLE_ThyraCoreLibs=OFF because parent package Trilinos_ENABLE_Thyra=OFF"
    "Setting subpackage enable Trilinos_ENABLE_ThyraGoodStuff=OFF because parent package Trilinos_ENABLE_Thyra=OFF" 
    "Setting subpackage enable Trilinos_ENABLE_ThyraEpetra=OFF because parent package Trilinos_ENABLE_Thyra=OFF"
    "Setting subpackage enable Trilinos_ENABLE_ThyraEpetraExt=OFF because parent package Trilinos_ENABLE_Thyra=OFF"
    "Setting subpackage enable Trilinos_ENABLE_ThyraTpetra=OFF because parent package Trilinos_ENABLE_Thyra=OFF"
    "WARNING: Setting Trilinos_ENABLE_Stratimikos=OFF which was 'ON' because Stratimikos has a required library dependence on disabled package ThyraCoreLibs"
    "Final set of enabled packages:  0"
    "Final set of non-enabled packages:  .+ Teuchos .+ Thyra .+ Stratimikos .+"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableEpetra_DisableThyraCoreLibs_EnableAllForwardPackages
  ARGS -DTrilinos_ENABLE_Epetra:BOOL=ON -DTrilinos_ENABLE_ThyraCoreLibs:BOOL=OFF
    -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos Epetra Zoltan Shards Triutils Tpetra EpetraExt Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML RBGen 15"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_DefaultTPLs
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK DUMMY 3"
    "Final set of non-enabled TPLs:  MPI Boost Scotch METIS ParMETIS CppUnit ADOLC ADIC TVMET y12m SuperLUDist SuperLU UMFPACK AMD PETSC MUMPS 16"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableBoost
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON -DTPL_ENABLE_Boost:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK Boost DUMMY 4"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableBoostParMETIS
  ARGS -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=ON -DTPL_ENABLE_Boost:BOOL=ON
    -DTPL_ENABLE_ParMETIS:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled TPLs:  BLAS LAPACK Boost ParMETIS DUMMY 5"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableZoltan_ZoltanEnableParMETIS
  ARGS  -DTrilinos_ENABLE_Zoltan=ON -DZoltan_ENABLE_ParMETIS=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Setting TPL_ENABLE_ParMETIS=ON because Zoltan_ENABLE_ParMETIS=ON"
    "Final set of enabled TPLs:  ParMETIS 1"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableZoltan_ZoltanEnableParMETIS_TplDisableParMETIS
  ARGS -DTrilinos_ENABLE_Zoltan=ON -DZoltan_ENABLE_ParMETIS=ON -DTPL_ENABLE_ParMETIS=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: Setting Zoltan_ENABLE_ParMETIS=OFF which was ON since TPL_ENABLE_ParMETIS=OFF"
    "Final set of enabled TPLs:  0"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableThyra_ThyraEnableEpetraExt
  ARGS  -DTrilinos_ENABLE_Thyra=ON -DTrilinos_ENABLE_ThyraEpetraExt=ON
    -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra EpetraExt Thyra 5"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnablePanzer_DefaultMPI
  ARGS  -DTrilinos_ENABLE_Panzer=ON
    -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "Setting Trilinos_ENABLE_Panzer=OFF which was 'ON' because Panzer has a required library dependence on disabled TPL MPI!"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnablePanzer_EnableMPI
  ARGS  -DTrilinos_ENABLE_Panzer=ON
    -DTPL_ENABLE_MPI:BOOL=ON
    -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  .*Panzer"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableTpetra_EnableAllForwardPackages
  ARGS  -DTrilinos_ENABLE_Tpetra=ON
    -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Thyra 8"
    "Final set of enabled SE packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt ThyraCoreLibs ThyraEpetra ThyraEpetraExt ThyraTpetra Thyra 12"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableAllPackages_EnableSecondaryStableCode_EnableAllOptionalPackges
  ARGS  -DTrilinos_ENABLE_ALL_PACKAGES=ON -DTrilinos_ENABLE_SECONDARY_STABLE_CODE=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Final set of enabled packages:  TrilinosFramework Teuchos RTOp Epetra Zoltan Shards Triutils Tpetra EpetraExt Sacado Thyra Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx 22"
    "Final set of enabled TPLs:  BLAS LAPACK Boost DUMMY 4"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableStratimikos_EnableAllOptionalPackages_BadSystem1
  ARGS  -DTrilinos_ENABLE_Stratimikos=ON
    -DTrilinos_HOSTTYPE:STRING=BadSystem1
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly disabled packages on input .+:  Stokhos Ifpack ML 3"
    "Setting Stratimikos_ENABLE_Ifpack=OFF because Stratimikos has an optional library dependence on disabled package Ifpack"
    "Setting Stratimikos_ENABLE_ML=OFF because Stratimikos has an optional library dependence on disabled package ML"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra AztecOO Galeri Amesos Stratimikos 11"
    "Final set of non-enabled packages:  TrilinosFramework Shards Tpetra Stokhos Sacado Isorropia Intrepid Ifpack ML Belos RBGen Phalanx Panzer 13"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableStratimikos_EnableAllOptionalPackages_BadSystem2
  ARGS  -DTrilinos_ENABLE_Stratimikos=ON
    -DTrilinos_HOSTTYPE:STRING=BadSystem2
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly disabled packages on input .+:  Stokhos Ifpack 2"
    "Setting Stratimikos_ENABLE_Ifpack=OFF because Stratimikos has an optional library dependence on disabled package Ifpack"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra Isorropia AztecOO Galeri Amesos ML Stratimikos 13"
    "Final set of enabled SE packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt ThyraCoreLibs ThyraEpetra ThyraEpetraExt Thyra Isorropia AztecOO Galeri Amesos ML Stratimikos 16"
    "Final set of non-enabled packages:  TrilinosFramework Shards Tpetra Stokhos Sacado Intrepid Ifpack Belos RBGen Phalanx Panzer 11"
    "Final set of non-enabled SE packages:  TrilinosFramework Shards Tpetra Stokhos Sacado ThyraGoodStuff ThyraCrazyStuff ThyraTpetra Intrepid Ifpack Belos RBGen Phalanx Panzer 14"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  EnableStratimikos_EnableIfpack_EnableAllOptionalPackages_BadSystem2
  ARGS  -DTrilinos_ENABLE_Stratimikos=ON -DTrilinos_ENABLE_Ifpack=ON
    -DTrilinos_HOSTTYPE:STRING=BadSystem2
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING: User has set Trilinos_ENABLE_Ifpack=ON but the" 
    "package Ifpack is not supported on this platform type 'BadSystem2'"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos 14"
    "Final set of enabled SE packages:  Teuchos RTOp Epetra Zoltan Triutils EpetraExt ThyraCoreLibs ThyraEpetra ThyraEpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos 17"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  preCopyrightTrilinos_EnableTeko
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
    -DTrilinos_ENABLE_Teko=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_PACKAGES='26'"
    "Trilinos_NUM_TPLS='20'"
    "Explicitly enabled packages on input .* Teko 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting Trilinos_ENABLE_Teuchos=ON because Teko has an optional dependence on Teuchos"
    "Setting Trilinos_ENABLE_Epetra=ON because Teko has an optional dependence on Epetra"
    "Setting Trilinos_ENABLE_Thyra=ON because Teko has an optional dependence on Thyra"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos Teko 16"
    "Final set of non-enabled packages:  TrilinosFramework Shards Stokhos Sacado Intrepid Belos RBGen Phalanx Panzer Stalix 10"
    "Final set of enabled TPLs:  BLAS LAPACK TekoDepTPL 3"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  preCopyrightTrilinos_EnableTeko_InconsistentExtraReposFile
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
    -DTrilinos_EXTRAREPOS_FILE=${CMAKE_CURRENT_SOURCE_DIR}/NoPackagesExtraReposFile.cmake
    -DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE=Nightly
    #-DTrilinos_VERBOSE_CONFIGURE=ON
    #-DTRIBITS_PROCESS_EXTRAREPOS_LISTS_DEBUG=ON
    -DTrilinos_ENABLE_Teko=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_PACKAGES='26'"
    "Trilinos_NUM_TPLS='20'"
    "Explicitly enabled packages on input .* Teko 1"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  preCopyrightTrilinos_EnableStalix
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
    -DTrilinos_ENABLE_Stalix=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .* Stalix 1"
    "Explicitly disabled packages on input .* Stokhos 1"
    "WARNING: Not Setting Trilinos_ENABLE_Sacado=ON even though Stalix has an optional dependence on Sacado because Trilinos_ENABLE_SECONDARY_STABLE_CODE=OFF"
    "Setting Trilinos_ENABLE_Tpetra=ON because Stalix has an optional dependence on Tpetra"
    "Setting Trilinos_ENABLE_Thyra=ON because Stalix has an optional dependence on Thyra"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Sacado Thyra Stalix 10"
    "Final set of non-enabled packages:  TrilinosFramework Shards Stokhos Isorropia AztecOO Galeri Amesos Intrepid Ifpack ML Belos Stratimikos RBGen Phalanx Panzer Teko 16"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  preCopyrightTrilinos_ShowInvalidPackageNameError
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
    -DSHOW_INVALID_PACKAGE_NAME_ERROR=ON
    -DTrilinos_ASSERT_MISSING_PACKAGES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Error, the package 'Anasazi' is listed as a dependency of the package"
    "'Teko' is in the list 'Trilinos_SE_PACKAGES' but the package 'Anasazi' is"
    "either not defined or is listed later in the package order.  Check the"
    "pelling of 'Anasazi' or see how it is listed in"
    "Trilinos_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS in relationship to 'Teko'"
    "Configuring incomplete, errors occurred"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  preCopyrightTrilinos_ShowInvalidPackageSelfDependency
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
    -DSHOW_SELF_PACKAGE_DEPENDENCY_ERROR=ON
    -DTrilinos_ASSERT_MISSING_PACKAGES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Error, the package 'Teko' is listed as a dependency itself in the list"
    " 'LIB_REQUIRED_DEP_PACKAGES'!"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  preCopyrightTrilinos_ShowInvalidSubpackageSelfDependency
  ARGS
    -DSHOW_SELF_SUBPACKAGE_DEPENDENCY_ERROR=ON
    -DTrilinos_ASSERT_MISSING_PACKAGES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Error, the package 'ThyraTpetra' is listed as a dependency itself in the"
     " list 'LIB_REQUIRED_DEP_PACKAGES'!"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  MissingExtraRepositoryFail
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=MissingExtraRepo
    -DTrilinos_IGNORE_MISSING_EXTRA_REPOSITORIES=OFF
  PASS_REGULAR_EXPRESSION_ALL
    "ERROR: Skipping missing extra repo 'MissingExtraRepo' packages list file"
    "ERROR: Skipping missing extra repo 'MissingExtraRepo' TPLs list file"
    "Configuring incomplete, errors occurred"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  MissingExtraRepositoryIgnore
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=MissingExtraRepo
    -DTrilinos_IGNORE_MISSING_EXTRA_REPOSITORIES=ON
  PASS_REGULAR_EXPRESSION_ALL
    "WARNING!  Ignoring missing extra repo 'MissingExtraRepo' packages list file ..+/tribits/package_arch/UnitTests/MockTrilinos/MissingExtraRepo/PackagesList.cmake. on request!"
    "WARNING!  Ignoring missing extra repo 'MissingExtraRepo' TPLs list file ..+/tribits/package_arch/UnitTests/MockTrilinos/MissingExtraRepo/TPLsList.cmake. on request!"
    "Configuring done"
    "Generating done"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  extraTrilinosRepo_EnableExtraPack
  ARGS "-DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos,extraTrilinosRepo"
    -DTrilinos_ENABLE_ExtraPack=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_TPLS='22'"
    "Explicitly enabled packages on input .* ExtraPack 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting Trilinos_ENABLE_Teuchos=ON because ExtraPack has an optional dependence on Teuchos"
    "Setting Trilinos_ENABLE_Teko=ON because Trilinos_ENABLE_ExtraPack=ON"
    "Final set of enabled packages:  Teuchos RTOp Epetra Triutils EpetraExt Thyra AztecOO Amesos Ifpack ML Stratimikos Teko ExtraPack 13"
    "Final set of non-enabled packages:  TrilinosFramework Zoltan Shards Tpetra Stokhos Sacado Isorropia Galeri Intrepid Belos RBGen Phalanx Panzer Stalix 14"
    "Final set of enabled TPLs:  BLAS LAPACK TekoDepTPL ExtraTPL2 4"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  extraTrilinosRepo_EnableExtraPack_EnableExtraTPL1
  ARGS "-DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos,extraTrilinosRepo"
    -DTrilinos_ENABLE_ExtraPack=ON -DTPL_ENABLE_ExtraTPL1=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Trilinos_NUM_TPLS='22'"
    "Explicitly enabled packages on input .* ExtraPack 1"
    "Final set of enabled TPLs:  BLAS LAPACK TekoDepTPL ExtraTPL1 ExtraTPL2 5"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  extraTrilinosRepo_ReplacePackagesList
  ARGS -DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
    -DpreCopyrightTrilinos_PACKAGES_LIST_FILE=cmake/Reduced_preCopyrightTrilinos_PackagesList.cmake
  PASS_REGULAR_EXPRESSION_ALL
    "Reading a list of extra packages from .*/MockTrilinos/cmake/Reduced_preCopyrightTrilinos_PackagesList.cmake"
    "Trilinos_NUM_PACKAGES='24'"
    "Final set of non-enabled SE packages:  .* Phalanx Panzer Stalix 31"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_CI_EnableTeko
  ARGS -DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=Continuous
    -DTrilinos_ENABLE_Teko=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .* Teko 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting Trilinos_ENABLE_Teuchos=ON because Teko has an optional dependence on Teuchos"
    "Setting Trilinos_ENABLE_Epetra=ON because Teko has an optional dependence on Epetra"
    "Setting Trilinos_ENABLE_Thyra=ON because Teko has an optional dependence on Thyra"
    "Final set of enabled packages:  Teuchos RTOp Epetra Zoltan Triutils Tpetra EpetraExt Thyra Isorropia AztecOO Galeri Amesos Ifpack ML Stratimikos Teko 16"
    "Final set of non-enabled packages:  TrilinosFramework Shards Stokhos Sacado Intrepid Belos RBGen Phalanx Panzer Stalix 10"
  )


CREATE_DEPENDENCY_HANDLING_TEST_CASE(
  ENABLE_KNOWN_EXTERNAL_REPOS_TYPE_Nightly_EnableExtraPack
  ARGS  -DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=Nightly
    -DTrilinos_ENABLE_ExtraPack=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Explicitly enabled packages on input .* ExtraPack 1"
    "Explicitly disabled packages on input .* Stokhos Stalix 2"
    "Setting Trilinos_ENABLE_Teuchos=ON because ExtraPack has an optional dependence on Teuchos"
    "Setting Trilinos_ENABLE_Teko=ON because Trilinos_ENABLE_ExtraPack=ON"
    "Final set of enabled packages:  Teuchos RTOp Epetra Triutils EpetraExt Thyra AztecOO Amesos Ifpack ML Stratimikos Teko ExtraPack 13"
    "Final set of non-enabled packages:  TrilinosFramework Zoltan Shards Tpetra Stokhos Sacado Isorropia Galeri Intrepid Belos RBGen Phalanx Panzer Stalix 14"
  )
