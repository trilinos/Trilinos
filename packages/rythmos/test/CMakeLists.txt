INCLUDE(AddSubdirectories)

#
# Add specific test executables
#

ADD_SUBDIRECTORIES(UnitTest ConvergenceTest Charon)

ASSERT_DEFINED( ${PACKAGE_NAME}_ENABLE_Experimental )
IF ( ${PACKAGE_NAME}_ENABLE_Experimental )
  ADD_SUBDIRECTORIES(ThetaStepper)
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Stratimikos)
IF (${PACKAGE_NAME}_ENABLE_Stratimikos  AND ${PACKAGE_NAME}_ENABLE_ThyraEpetraExtAdapters)
  ADD_SUBDIRECTORIES(complicatedExample)
  ADD_SUBDIRECTORIES(simpleAdjoint)
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Sacado Thyra_ENABLE_ME_POLYNOMIAL)
IF (${PACKAGE_NAME}_ENABLE_Sacado AND
    ${PACKAGE_NAME}_ENABLE_ThryaEpetraAdapters AND
    Thyra_ENABLE_ME_POLYNOMIAL)
  ADD_SUBDIRECTORIES(ExplicitTaylorPolynomialStepper)
ENDIF()
