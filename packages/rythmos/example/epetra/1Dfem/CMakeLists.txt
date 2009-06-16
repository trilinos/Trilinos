
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)
INCLUDE(MultilineSet)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

FILE(GLOB SOURCES *.cpp *.C)

PACKAGE_ADD_EXECUTABLE(
  1DfemTransient
  SOURCES ${SOURCES}
  DEPLIBS rythmos_test_models
  COMM serial mpi  
  )

ASSERT_DEFINED(Stratimikos_ENABLE_Amesos)
IF (Stratimikos_ENABLE_Amesos)

  SET(STANDARD_TEST_CONFIG
    COMM serial mpi  
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --echo-command-line "
    )

  PACKAGE_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_amesos_BE
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --method=BDF "
    " --echo-command-line "
    )

  PACKAGE_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_amesos_IBDF
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --method=BDF "
    " --useNOX "
    " --echo-command-line "
    )

  PACKAGE_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_amesos_IBDF_NOX
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

ENDIF()


ASSERT_DEFINED(Stratimikos_ENABLE_AztecOO)
IF (Stratimikos_ENABLE_AztecOO)

  SET(STANDARD_TEST_CONFIG
    COMM serial mpi  
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --extra-linear-solver-params=\""
    "   <ParameterList> "
    "     <Parameter name='Linear Solver Type' type='string' value='AztecOO'/> "
    "     <ParameterList name ='Linear Solver Types'> "
    "       <ParameterList name='AztecOO'> "
    "         <ParameterList name='Forward Solve'> "
    "           <Parameter name='Tolerance' type='double' value='1e-12'/> "
    "         </ParameterList> "
    "       </ParameterList> "
    "     </ParameterList> "
    "   </ParameterList> \"" 
    " --linear-solver-params-used-file=lowsf.aztecoo.used.xml "
    " --echo-command-line "
    )

  PACKAGE_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_aztecoo_BE
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --method=BDF "
    " --extra-linear-solver-params=\""
    "   <ParameterList> "
    "     <Parameter name='Linear Solver Type' type='string' value='AztecOO'/> "
    "     <ParameterList name ='Linear Solver Types'> "
    "       <ParameterList name='AztecOO'> "
    "         <ParameterList name='Forward Solve'> "
    "           <Parameter name='Tolerance' type='double' value='1e-12'/> "
    "         </ParameterList> "
    "       </ParameterList> "
    "     </ParameterList> "
    "   </ParameterList> \"" 
    " --linear-solver-params-used-file=lowsf.aztecoo.used.xml "
    " --echo-command-line "
    )

  PACKAGE_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_aztecoo_IBDF
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --method=IRK "
    " --numsteps=5 "
    " --extra-linear-solver-params=\""
    "   <ParameterList> "
    "     <Parameter name='Linear Solver Type' type='string' value='AztecOO'/> "
    "     <ParameterList name ='Linear Solver Types'> "
    "       <ParameterList name='AztecOO'> "
    "         <ParameterList name='Forward Solve'> "
    "           <Parameter name='Tolerance' type='double' value='1e-12'/> "
    "         </ParameterList> "
    "       </ParameterList> "
    "     </ParameterList> "
    "   </ParameterList> \"" 
    " --linear-solver-params-used-file=lowsf.aztecoo.used.xml "
    " --echo-command-line "
    )

  PACKAGE_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_aztecoo_IRK
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --method=BDF "
    " --useNOX "
    " --extra-linear-solver-params=\""
    "   <ParameterList> "
    "     <Parameter name='Linear Solver Type' type='string' value='AztecOO'/> "
    "     <ParameterList name ='Linear Solver Types'> "
    "       <ParameterList name='AztecOO'> "
    "         <ParameterList name='Forward Solve'> "
    "           <Parameter name='Tolerance' type='double' value='1e-12'/> "
    "         </ParameterList> "
    "       </ParameterList> "
    "     </ParameterList> "
    "   </ParameterList> \"" 
    " --linear-solver-params-used-file=lowsf.aztecoo.used.xml "
    " --echo-command-line "
    )

  PACKAGE_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_aztecoo_IBDF_NOX
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

ENDIF()


ASSERT_DEFINED(Stratimikos_ENABLE_Belos)
IF (Stratimikos_ENABLE_Belos)

  SET(STANDARD_TEST_CONFIG
    COMM serial mpi  
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --extra-linear-solver-params=\""
    "   <ParameterList> "
    "     <Parameter name='Linear Solver Type' type='string' value='Belos'/> "
    "     <ParameterList name='Linear Solver Types'> "
    "       <ParameterList name='Belos'> "
    "         <Parameter name='Solver Type' type='string' value='Block GMRES'/> "
    "         <ParameterList name='Solver Types'> "
    "           <ParameterList name='Block GMRES'> "
    "             <Parameter name='Convergence Tolerance' type='double' value='1e-12'/> "
    "           </ParameterList> "
    "         </ParameterList> "
    "       </ParameterList> "
    "     </ParameterList> "
    "   </ParameterList> \""
    "  --linear-solver-params-used-file=lowsf.belos.used.xml "
    )
  PACKAGE_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_belos_BE
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --method=BDF "
    " --extra-linear-solver-params=\""
    "   <ParameterList> "
    "     <Parameter name='Linear Solver Type' type='string' value='Belos'/> "
    "     <ParameterList name='Linear Solver Types'> "
    "       <ParameterList name='Belos'> "
    "         <Parameter name='Solver Type' type='string' value='Block GMRES'/> "
    "         <ParameterList name='Solver Types'> "
    "           <ParameterList name='Block GMRES'> "
    "             <Parameter name='Convergence Tolerance' type='double' value='1e-12'/> "
    "           </ParameterList> "
    "         </ParameterList> "
    "       </ParameterList> "
    "     </ParameterList> "
    "   </ParameterList> \""
    "  --linear-solver-params-used-file=lowsf.belos.used.xml "
    )
  PACKAGE_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_belos_IBDF
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --method=BDF "
    " --useNOX "
    " --extra-linear-solver-params=\""
    "   <ParameterList> "
    "     <Parameter name='Linear Solver Type' type='string' value='Belos'/> "
    "     <ParameterList name='Linear Solver Types'> "
    "       <ParameterList name='Belos'> "
    "         <Parameter name='Solver Type' type='string' value='Block GMRES'/> "
    "         <ParameterList name='Solver Types'> "
    "           <ParameterList name='Block GMRES'> "
    "             <Parameter name='Convergence Tolerance' type='double' value='1e-12'/> "
    "           </ParameterList> "
    "         </ParameterList> "
    "       </ParameterList> "
    "     </ParameterList> "
    "   </ParameterList> \""
    "  --linear-solver-params-used-file=lowsf.belos.used.xml "
    )
  PACKAGE_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_belos_IBDF_NOX
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

ENDIF()
