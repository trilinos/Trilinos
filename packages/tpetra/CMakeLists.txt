
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)
INCLUDE(CombinedOption)

#
# A) Define the package
#

ASSERT_DEFINED(Kokkos_ENABLE_TBB)
IF (Kokkos_ENABLE_TBB)
  # ENABLE_SHADOWING_WARNINGS not possible, because of TBB
  PACKAGE( Tpetra CLEANED )
ELSE()
  # ToDo: Enable this once Kokkos headers and Tpetra have shadowing
  # warnings cleaned.
  #PACKAGE( Tpetra ENABLE_SHADOWING_WARNINGS CLEANED )
  PACKAGE( Tpetra CLEANED )
ENDIF()

COMBINED_OPTION( Tpetra_ENABLE_Thrust
  DEP_OPTIONS_NAMES Kokkos_ENABLE_Thrust
  DOCSTR "Enables CUDA GPU support in Tpetra.  Requires Kokkos/Thrust"
  )

PACKAGE_ADD_EXPLICIT_INSTANTIATION_OPTION()

PACKAGE_ADD_DEBUG_OPTION()

PACKAGE_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#
# B) Set up package-specific options
#

PACKAGE_ARCH_SET_SS_FOR_DEV_MODE(SS_FOR_DEV_PS_FOR_RELEASE)

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_SS_TESTING
  HAVE_TPETRA_ENABLE_SS_TESTING
  "Enabling testing of Secondary Stable Code."
  ${${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE}
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_BUGTESTS
  HAVE_TPETRA_BUGTESTS
  "Enabling testing of previous Tpetra bugs."
  ${${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE}
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_RTI
  HAVE_TPETRA_RTI
  "Enabling Tpetra Reduction/Transformation Interface."
  ${${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE}
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_FLOAT
  HAVE_TPETRA_INST_FLOAT
  "If explicit instantiation is enabled, instantiate over scalar field float."
  ${SS_FOR_DEV_PS_FOR_RELEASE}
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_DOUBLE
  HAVE_TPETRA_INST_DOUBLE
  "If explicit instantiation is enabled, instantiate over scalar field double."
  ON
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Threaded_MKL
  HAVE_TPETRA_THREADED_MKL
  "Indicates that the linked BLAS is a threaded version of the MKL."
  OFF
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_COMPLEX_FLOAT
  HAVE_TPETRA_INST_COMPLEX_FLOAT
  "If explicit instantiation is enabled, instantiate over scalar field complex<float>."
  ${SS_FOR_DEV_PS_FOR_RELEASE}
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_INST_COMPLEX_DOUBLE
  HAVE_TPETRA_INST_COMPLEX_DOUBLE
  "If explicit instantiation is enabled, instantiate over scalar field complex<double>."
  ${SS_FOR_DEV_PS_FOR_RELEASE}
  )

ADD_OPTION_AND_DEFINE(
  Tpetra_THROW_Warnings
  HAVE_TPETRA_THROW_WARNINGS
  "Enable exception throwing for a number of warnings in Tpetra."
  OFF
  )

ADD_OPTION_AND_DEFINE(
  Tpetra_PRINT_Warnings
  HAVE_TPETRA_PRINT_WARNINGS
  "Enable printing of a number of warnings in Tpetra."
  OFF
  )

ADD_OPTION_AND_DEFINE(
  Tpetra_THROW_Efficiency_Warnings
  HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS
  "Enable exception throwing for Tpetra efficiency warnings."
  ${Tpetra_THROW_Warnings}
  )

ADD_OPTION_AND_DEFINE(
  Tpetra_PRINT_Efficiency_Warnings
  HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS
  "Enable printing of Tpetra efficiency warnings."
  ${Tpetra_PRINT_Warnings}
  )

ADD_OPTION_AND_DEFINE(
  Tpetra_THROW_Abuse_Warnings
  HAVE_TPETRA_THROW_ABUSE_WARNINGS
  "Enable exception throwing for potential Tpetra abuse warnings."
  ${Tpetra_THROW_Warnings}
  )

ADD_OPTION_AND_DEFINE(
  Tpetra_PRINT_Abuse_Warnings
  HAVE_TPETRA_PRINT_ABUSE_WARNINGS
  "Enable printing of potential Tpetra abuse warnings."
  ${Tpetra_PRINT_Warnings}
  )

ADD_OPTION_AND_DEFINE(
  Tpetra_ENABLE_TSQR
  HAVE_TPETRA_TSQR
  "Enable TSQR support in Tpetra."
  ${Kokkos_ENABLE_TSQR}
  )


#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(inout)
ADD_SUBDIRECTORY(ext)
IF (Tpetra_ENABLE_RTI)
  ADD_SUBDIRECTORY(rti)
ENDIF()

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_TEST_DIRECTORIES(example)


#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

PACKAGE_ARCH_EXCLUDE_FILES(
  doc/CodingGuidelines
  doc/TpetraDesign
  )

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
