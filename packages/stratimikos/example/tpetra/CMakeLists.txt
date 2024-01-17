tribits_add_executable( simple_tpetra_stratimikos_example
  SOURCES simple_tpetra_stratimikos_example.cpp )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos2)

IF (${PACKAGE_NAME}_ENABLE_Amesos2)
ASSERT_DEFINED(Amesos2_ENABLE_KLU2)
TRIBITS_ADD_TEST( simple_tpetra_stratimikos_example
  EXCLUDE_IF_NOT_TRUE Amesos2_ENABLE_KLU2
  POSTFIX_AND_ARGS_0 amesos2_klu2
    "--linear-solver-params-file=${CMAKE_CURRENT_SOURCE_DIR}/amesos2.klu2.xml"
    "--matrix-file=${CMAKE_CURRENT_SOURCE_DIR}/FourByFour.mtx" --tol=1e-13
    --linear-solver-params-used-file=lowsf.amesos2.klu2.used.xml
  NUM_MPI_PROCS 1
  )
ENDIF()

TRIBITS_ADD_TEST( simple_tpetra_stratimikos_example
  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Belos
  POSTFIX_AND_ARGS_0 belos_noprec
    "--linear-solver-params-file=${CMAKE_CURRENT_SOURCE_DIR}/belos.noprec.xml"
    "--matrix-file=${CMAKE_CURRENT_SOURCE_DIR}/FourByFour.mtx" --tol=1e-13
    --linear-solver-params-used-file=lowsf.belos.noprec.used.xml
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST( simple_tpetra_stratimikos_example
  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Belos ${PACKAGE_NAME}_ENABLE_Ifpack2
  POSTFIX_AND_ARGS_0 belos_ifpack2
    "--linear-solver-params-file=${CMAKE_CURRENT_SOURCE_DIR}/belos.ifpack2.xml"
    "--matrix-file=${CMAKE_CURRENT_SOURCE_DIR}/FourByFour.mtx" --tol=1e-13
    --linear-solver-params-used-file=lowsf.belos.ifpack2.used.xml
  NUM_MPI_PROCS 1
  )
