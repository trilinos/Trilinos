INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(CombinedOption)
INCLUDE(TribitsETISupport)

#
# A) Define the package
#

ASSERT_DEFINED(KokkosClassic_ENABLE_TBB)
IF (KokkosClassic_ENABLE_TBB)
  # ENABLE_SHADOWING_WARNINGS not possible, because of TBB
  TRIBITS_PACKAGE( Tpetra CLEANED )
ELSE()
  # ToDo: Enable this once Kokkos headers and Tpetra have shadowing
  # warnings cleaned.
  #TRIBITS_PACKAGE( Tpetra ENABLE_SHADOWING_WARNINGS CLEANED )
  TRIBITS_PACKAGE( Tpetra CLEANED )
ENDIF()

COMBINED_OPTION( Tpetra_ENABLE_Thrust
  DEP_OPTIONS_NAMES KokkosClassic_ENABLE_Thrust
  DOCSTR "Enables CUDA GPU support in Tpetra.  Requires Kokkos/Thrust"
  )

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

# we use this for testing as well as library eti, so do it regardless
TRIBITS_ADD_ETI_SUPPORT()

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

SET(Tpetra_MACHINE_XML_FILE_DIR 
    ${CMAKE_CURRENT_SOURCE_DIR}/machine_files
    CACHE INTERNAL "")

#
# B) Set up package-specific options
#

TRIBITS_SET_SS_FOR_DEV_MODE(SS_FOR_DEV_PS_FOR_RELEASE)

ASSERT_DEFINED(KokkosClassic_DefaultNode)
ASSERT_DEFINED(KokkosClassic_ENABLED_NODES)
SET(Tpetra_ETI_SCALARS "double;int")
SET(Tpetra_ETI_GORDS   "int")
SET(Tpetra_ETI_LORDS   "int")
SET(Tpetra_ETI_NODES   "${KokkosClassic_ENABLED_NODES}")

# double is always built; these variables are needed for backwards # compatibility
GLOBAL_SET(HAVE_TPETRA_INST_DOUBLE ON)
GLOBAL_SET(Tpetra_INST_DOUBLE ON)

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_SS_TESTING
  HAVE_TPETRA_ENABLE_SS_TESTING
  "Enabling testing of Secondary Stable Code."
  ${${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_ENABLE_OpenMP
  HAVE_TPETRA_OPENMP
  "Enabling OpenMP in Tpetra"
  ${KokkosClassic_ENABLE_OpenMP}
  )
IF(Tpetra_ENABLE_OpenMP AND NOT KokkosClassic_ENABLE_OpenMP)
  MESSAGE(WARNING "Disabling OpenMP in Tpetra because KokkosClassic_ENABLE_OpenMP is ${KokkosClassic_ENABLE_OpenMP}")
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_BUGTESTS
  HAVE_TPETRA_BUGTESTS
  "Enabling testing of previous Tpetra bugs."
  ${${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_RTI
  HAVE_TPETRA_RTI
  "Enabling Tpetra Reduction/Transformation Interface."
  ${${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE}
  )
IF (${PACKAGE_NAME}_ENABLE_RTI)
  MESSAGE(STATUS "Enabling Tpetra/RTI")
ELSE()
  MESSAGE(STATUS "Not enabling Tpetra/RTI")
ENDIF()

IF (SS_FOR_DEV_PS_FOR_RELEASE AND HAVE_TEUCHOS_BLASFLOAT) 
  GLOBAL_SET(HAVE_TPETRA_INST_FLOAT_DEFAULT ON)
ELSE()
  GLOBAL_SET(HAVE_TPETRA_INST_FLOAT_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_FLOAT
  HAVE_TPETRA_INST_FLOAT
  "If explicit instantiation is enabled, instantiate over scalar field float."
  ${HAVE_TPETRA_INST_FLOAT_DEFAULT}
  )
IF(Tpetra_INST_FLOAT) 
  LIST(APPEND Tpetra_ETI_SCALARS "float")
ENDIF()

IF (SS_FOR_DEV_PS_FOR_RELEASE AND HAVE_COMPLEX_BLAS AND HAVE_TEUCHOS_BLASFLOAT)
  GLOBAL_SET(HAVE_TPETRA_INST_COMPLEX_FLOAT_DEFAULT ON)
ELSE()
  GLOBAL_SET(HAVE_TPETRA_INST_COMPLEX_FLOAT_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_COMPLEX_FLOAT
  HAVE_TPETRA_INST_COMPLEX_FLOAT
  "If explicit instantiation is enabled, instantiate over scalar field std::complex<float>."
  ${HAVE_TPETRA_INST_COMPLEX_FLOAT_DEFAULT}  
  )
IF(Tpetra_INST_COMPLEX_FLOAT) 
  LIST(APPEND Tpetra_ETI_SCALARS "std::complex<float>")
ENDIF()

IF (SS_FOR_DEV_PS_FOR_RELEASE AND HAVE_COMPLEX_BLAS)
  GLOBAL_SET(HAVE_TPETRA_INST_COMPLEX_DOUBLE_DEFAULT ON)
ELSE()
  GLOBAL_SET(HAVE_TPETRA_INST_COMPLEX_DOUBLE_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_COMPLEX_DOUBLE
  HAVE_TPETRA_INST_COMPLEX_DOUBLE
  "If explicit instantiation is enabled, instantiate over scalar field std::complex<double>."
  ${HAVE_TPETRA_INST_COMPLEX_DOUBLE_DEFAULT}  
  )
IF(Tpetra_INST_COMPLEX_DOUBLE) 
  LIST(APPEND Tpetra_ETI_SCALARS "std::complex<double>")
ENDIF()

IF (SS_FOR_DEV_PS_FOR_RELEASE AND Tpetra_ENABLE_QD)
  GLOBAL_SET(HAVE_TPETRA_INST_QD_REAL_DEFAULT ON)
ELSE()
  GLOBAL_SET(HAVE_TPETRA_INST_QD_REAL_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_QD_REAL
  HAVE_TPETRA_INST_QD_REAL
  "If explicit instantiation is enabled, instantiate over scalar field qd_real (requries TPL QD)."
  ${HAVE_TPETRA_INST_QD_REAL_DEFAULT}  
  )
IF(Tpetra_INST_QD_REAL) 
  LIST(APPEND Tpetra_ETI_SCALARS "qd_real")
ENDIF()

IF (SS_FOR_DEV_PS_FOR_RELEASE AND Tpetra_ENABLE_QD)
  GLOBAL_SET(HAVE_TPETRA_INST_DD_REAL_DEFAULT ON)
ELSE()
  GLOBAL_SET(HAVE_TPETRA_INST_DD_REAL_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_DD_REAL
  HAVE_TPETRA_INST_DD_REAL
  "If explicit instantiation is enabled, instantiate over scalar field dd_real (requries TPL QD)."
  ${HAVE_TPETRA_INST_DD_REAL_DEFAULT}  
  )
IF(Tpetra_INST_DD_REAL) 
  LIST(APPEND Tpetra_ETI_SCALARS "dd_real")
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_INT_LONG
  HAVE_TPETRA_INST_INT_LONG
  "If explicit instantiation is enabled, instantiate over global ordinal long."
  ON
  )
IF(Tpetra_INST_INT_LONG) 
  LIST(APPEND Tpetra_ETI_GORDS   "long")
  LIST(APPEND Tpetra_ETI_SCALARS "long")
ENDIF()

SET(Tpetra_INST_UNSIGNED_LONG
  ${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE}
  CACHE BOOL
  "If explicit instantiation is enabled, instatiate over global ordinal unsigned int."
  )
IF (Tpetra_INST_UNSIGNED_LONG)
  LIST(APPEND Tpetra_ETI_GORDS   "unsigned int")
  LIST(APPEND Tpetra_ETI_SCALARS "unsigned int")
  # NOTE: The above types are used because they help catch errors with caused
  # when GlobalOrdinal is unsigned and GlobalOrdinal is smaller than
  # global_size_t (using wrong OrdinalTraits<>::invald(). 
ENDIF()



TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Threaded_MKL
  HAVE_TPETRA_THREADED_MKL
  "Indicates that the linked BLAS is a threaded version of the MKL."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_THROW_Warnings
  HAVE_TPETRA_THROW_WARNINGS
  "Enable exception throwing for a number of warnings in Tpetra."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_PRINT_Warnings
  HAVE_TPETRA_PRINT_WARNINGS
  "Enable printing of a number of warnings in Tpetra."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_THROW_Efficiency_Warnings
  HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS
  "Enable exception throwing for Tpetra efficiency warnings."
  ${Tpetra_THROW_Warnings}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_PRINT_Efficiency_Warnings
  HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS
  "Enable printing of Tpetra efficiency warnings."
  ${Tpetra_PRINT_Warnings}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_THROW_Abuse_Warnings
  HAVE_TPETRA_THROW_ABUSE_WARNINGS
  "Enable exception throwing for potential Tpetra abuse warnings."
  ${Tpetra_THROW_Warnings}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_PRINT_Abuse_Warnings
  HAVE_TPETRA_PRINT_ABUSE_WARNINGS
  "Enable printing of potential Tpetra abuse warnings."
  ${Tpetra_PRINT_Warnings}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_ENABLE_TSQR
  HAVE_TPETRA_TSQR
  "Enable TSQR support in Tpetra."
  ${KokkosClassic_ENABLE_TSQR}
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

TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_ADD_TEST_DIRECTORIES(example)

GLOBAL_SET(Tpetra_ETI_SCALARS ${Tpetra_ETI_SCALARS})
GLOBAL_SET(Tpetra_ETI_GORDS   ${Tpetra_ETI_GORDS})
GLOBAL_SET(Tpetra_ETI_LORDS   ${Tpetra_ETI_LORDS})
GLOBAL_SET(Tpetra_ETI_NODES   ${Tpetra_ETI_NODES})

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()
TRIBITS_EXCLUDE_FILES(
  doc/CodingGuidelines
  doc/TpetraDesign
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
