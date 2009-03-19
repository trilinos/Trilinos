INCLUDE(AddSubdirectories)

#
# Add specific test executables
#

ADD_SUBDIRECTORIES(UnitTest ConvergenceTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Stratimikos)
IF (${PACKAGE_NAME}_ENABLE_Stratimikos)
  ADD_SUBDIRECTORIES(complicatedExample)
  ADD_SUBDIRECTORIES(simpleAdjoint)
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Sacado)
IF (${PACKAGE_NAME}_ENABLE_Sacado)
  ADD_SUBDIRECTORIES(ExplicitTaylorPolynomialStepper)
ENDIF()

