INCLUDE(CombinedOption)
INCLUDE(TribitsETISupport)

#
# Declare the subpackage
#
TRIBITS_SUBPACKAGE(Core)

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

# Enable Explicit Template Instantiation (ETI) support for Tpetra.
TRIBITS_ADD_ETI_SUPPORT()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
# FIXME (mfh 17 Dec 2014) We really just want a macro, not a CMake option.
# This is mainly just for backwards compatibility.
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_ENABLE_EXPLICIT_INSTANTIATION
  HAVE_TPETRA_EXPLICIT_INSTANTIATION
  "Enable explicit template instantiation (ETI) in Tpetra"
  ${${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION}
  )

SET(Tpetra_MACHINE_XML_FILE_DIR 
    ${CMAKE_CURRENT_SOURCE_DIR}/machine_files
    CACHE INTERNAL "")

ASSERT_DEFINED(KokkosClassic_DefaultNode)
ASSERT_DEFINED(KokkosClassic_ENABLED_NODES)
SET(TpetraCore_ETI_SCALARS "double;int")
SET(TpetraCore_ETI_GORDS   "int")
SET(TpetraCore_ETI_LORDS   "int")
SET(TpetraCore_ETI_NODES   "${KokkosClassic_ENABLED_NODES}")

# double is always built; these variables are needed for backwards # compatibility
GLOBAL_SET(HAVE_TPETRA_INST_DOUBLE ON)
GLOBAL_SET(Tpetra_INST_DOUBLE ON)

# int as global ordinal is also always built
GLOBAL_SET(Tpetra_INST_INT_INT ON)
GLOBAL_SET(HAVE_TPETRA_INST_INT_INT ON)

#Set Variable which were implicitly set when the Kokkos packages were optional dependencies
IF(Tpetra_ENABLE_Kokkos_Refactor)
GLOBAL_SET(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT ON)
GLOBAL_SET(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMM ON)
GLOBAL_SET(HAVE_TPETRACORE_KOKKOSCORE ON)
GLOBAL_SET(HAVE_TPETRACORE_KOKKOSALGORITHMS ON)
GLOBAL_SET(HAVE_TPETRACORE_KOKKOSCONTAINERS ON)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_SS_TESTING
  HAVE_TPETRA_ENABLE_SS_TESTING
  "Enabling testing of Secondary Stable Code."
  ${${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE}
  )

# If this variable isn't defined, then the TpetraTSQR subpackage of
# Tpetra isn't where it should be.
ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_TpetraTSQR)
IF (${PROJECT_NAME}_ENABLE_TpetraTSQR)
  SET(${PACKAGE_NAME}_ENABLE_TSQR_DEFAULT ON)
ELSE ()
  SET(${PACKAGE_NAME}_ENABLE_TSQR_DEFAULT OFF)
ENDIF ()
# TpetraCore_ENABLE_TSQR says whether TpetraCore should build its TSQR
# interface.  You may turn this interface off if you like, even if the
# TpetraTSQR subpackage of Tpetra is enabled.
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_TSQR
  HAVE_TPETRA_TSQR
  "Whether to build Tpetra with TSQR support.  This is enabled by default if the KokkosTSQR subpackage is enabled.  You may turn off TSQR support if you like, but leaving it enabled is fine.  Please note that Epetra's TSQR adapter lives in the Tpetra package, for various historical reasons.  Therefore, in order to use TSQR with Epetra, you must enable Tpetra as well."
  ${${PACKAGE_NAME}_ENABLE_TSQR_DEFAULT}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_CLASSIC_VBR
  HAVE_TPETRA_CLASSIC_VBR
  "Whether to build Tpetra with the VbrMatrix, BlockMap, BlockCrsGraph, and BlockMultiVector classes enabled.  These classes implement variable-block-size block sparse matrices, block dense multivectors, and other related technology.  (These differ from the constant-block-size classes currently in the Tpetra::Experimental namespace.)  This option is ON (enabled) by default.  Turning OFF this option may be useful for testing compilation with new Kokkos Node types.  The \"CLASSIC\" in the name suggests that we may revisit these classes after finishing moving Tpetra to use new Kokkos."
  ON
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_BUGTESTS
  HAVE_TPETRA_BUGTESTS
  "Enabling testing of previous Tpetra bugs."
  ${${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE}
  )

# Enable explicit template instantiation (ETI) and tests for Scalar =
# float.
#
# ST_FOR_DEV_PT_FOR_RELEASE means something like, "in a debug build."
# In a debug build, we do ETI for Scalar=float, in order that this
# capability be tested.  In a release build, we turn off ETI for
# Scalar=float, to reduce build time.
IF (ST_FOR_DEV_PT_FOR_RELEASE AND HAVE_TEUCHOS_BLASFLOAT) 
  GLOBAL_SET(HAVE_TPETRA_INST_FLOAT_DEFAULT ON)
ELSE()
  GLOBAL_SET(HAVE_TPETRA_INST_FLOAT_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_FLOAT
  HAVE_TPETRA_INST_FLOAT
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes over Scalar = float.  This option is ${HAVE_TPETRA_INST_FLOAT_DEFAULT} by default."
  ${HAVE_TPETRA_INST_FLOAT_DEFAULT}
  )
IF(Tpetra_INST_FLOAT AND NOT HAVE_TEUCHOS_BLASFLOAT)
  MESSAGE(ERROR "You may not instantiate Tpetra classes with Scalar = float, if the BLAS library you are using does not support float (S) arithmetic.")
ENDIF()
IF(Tpetra_INST_FLOAT) 
  LIST(APPEND TpetraCore_ETI_SCALARS "float")
ENDIF()

# Enable explicit template instantiation (ETI) and tests for Scalar =
# std::complex<float>.
#
# ST_FOR_DEV_PT_FOR_RELEASE means something like, "in a debug build."
# In a debug build, we do ETI for Scalar=std::complex<float>, in order
# that this capability be tested.  In a release build, we turn off ETI
# for Scalar=std::complex<float>, to reduce build time.
IF (ST_FOR_DEV_PT_FOR_RELEASE AND HAVE_COMPLEX_BLAS AND HAVE_TEUCHOS_BLASFLOAT AND NOT Tpetra_ENABLE_Kokkos_Refactor)
  # mfh 09 Jan 2015: std::complex<float> does NOT currently work with
  # the Kokkos refactor version of Tpetra, due to
  # Kokkos::atomic_fetch_add build issues.  Feel free to remove "AND
  # NOT Tpetra_ENABLE_Kokkos_Refactor" from the above test once that
  # gets fixed.
  GLOBAL_SET(HAVE_TPETRA_INST_COMPLEX_FLOAT_DEFAULT ON)
ELSE()
  GLOBAL_SET(HAVE_TPETRA_INST_COMPLEX_FLOAT_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_COMPLEX_FLOAT
  HAVE_TPETRA_INST_COMPLEX_FLOAT
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes with Scalar = std::complex<float>.  This option is ${HAVE_TPETRA_INST_COMPLEX_FLOAT_DEFAULT} by default."
  ${HAVE_TPETRA_INST_COMPLEX_FLOAT_DEFAULT}  
  )
IF(Tpetra_INST_COMPLEX_FLOAT AND (NOT HAVE_COMPLEX_BLAS OR NOT HAVE_TEUCHOS_BLASFLOAT))
  MESSAGE(ERROR "You may not instantiate Tpetra classes with Scalar = std::complex<float>, if the BLAS library you are using does not support complex (C) or float (S) arithmetic.")
ENDIF()
IF(Tpetra_INST_COMPLEX_FLOAT AND Tpetra_ENABLE_Kokkos_Refactor)
  MESSAGE(ERROR "Scalar = std::complex<float> does not currently work with the Kokkos refactor version of Tpetra.  We will fix this at some point.  If you need this capability urgently, please ask the Kokkos developers to fix Kokkos::atomic_fetch_add for Kokkos::complex<float>.")
ENDIF()
IF(Tpetra_INST_COMPLEX_FLOAT) 
  LIST(APPEND TpetraCore_ETI_SCALARS "std::complex<float>")
ENDIF()

# Enable explicit template instantiation (ETI) and tests for Scalar =
# std::complex<double>.
#
# ST_FOR_DEV_PT_FOR_RELEASE means something like, "in a debug build."
# In a debug build, we do ETI for Scalar=std::complex<float>, in order
# that this capability be tested.  In a release build, we turn off ETI
# for Scalar=std::complex<float>, to reduce build time.
IF (ST_FOR_DEV_PT_FOR_RELEASE AND HAVE_COMPLEX_BLAS)
  GLOBAL_SET(HAVE_TPETRA_INST_COMPLEX_DOUBLE_DEFAULT ON)
ELSE()
  GLOBAL_SET(HAVE_TPETRA_INST_COMPLEX_DOUBLE_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_COMPLEX_DOUBLE
  HAVE_TPETRA_INST_COMPLEX_DOUBLE
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes with Scalar = std::complex<double>.  This option is ${HAVE_TPETRA_INST_COMPLEX_DOUBLE_DEFAULT} by default."
  ${HAVE_TPETRA_INST_COMPLEX_DOUBLE_DEFAULT}  
  )
IF(Tpetra_INST_COMPLEX_DOUBLE AND NOT HAVE_COMPLEX_BLAS)
  MESSAGE(ERROR "You may not instantiate Tpetra classes with Scalar = std::complex<double>, if the BLAS library you are using does not support double-complex (Z) arithmetic.")
ENDIF()
IF(Tpetra_INST_COMPLEX_DOUBLE) 
  LIST(APPEND TpetraCore_ETI_SCALARS "std::complex<double>")
ENDIF()

# Enable explicit template instantiation (ETI) and tests for Scalar =
# qd_real ("quad-double," provided by the QD TPL).
#
# ST_FOR_DEV_PT_FOR_RELEASE means something like, "in a debug build."
# In a debug build, we do ETI for Scalar=std::complex<float>, in order
# that this capability be tested.  In a release build, we turn off ETI
# for Scalar=std::complex<float>, to reduce build time.
IF (ST_FOR_DEV_PT_FOR_RELEASE AND ${PACKAGE_NAME}_ENABLE_QD AND NOT Tpetra_ENABLE_Kokkos_Refactor)
  # mfh 09 Jan 2015: qd_real does NOT currently work with the Kokkos
  # refactor version of Tpetra, due to a missing Kokkos::atomic_*
  # implementation for this class.  Feel free to remove "AND NOT
  # Tpetra_ENABLE_Kokkos_Refactor" from the above test once that gets
  # fixed.
  GLOBAL_SET(HAVE_TPETRA_INST_QD_REAL_DEFAULT ON)
ELSE()
  GLOBAL_SET(HAVE_TPETRA_INST_QD_REAL_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_QD_REAL
  HAVE_TPETRA_INST_QD_REAL
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes over Scalar = qd_real (requries TPL QD).  This option is ${HAVE_TPETRA_INST_QD_REAL_DEFAULT} by default."
  ${HAVE_TPETRA_INST_QD_REAL_DEFAULT}  
  )
IF(Tpetra_INST_QD_REAL AND Tpetra_ENABLE_Kokkos_Refactor)
  MESSAGE(ERROR "Scalar = qd_real does not currently work with the Kokkos refactor version of Tpetra.  We will fix this at some point.  Please let us know if you need this capability urgently.")
ENDIF()
IF(Tpetra_INST_QD_REAL) 
  LIST(APPEND TpetraCore_ETI_SCALARS "qd_real")
ENDIF()

# Enable explicit template instantiation (ETI) and tests for Scalar =
# qd_real ("double-double," provided by the QD TPL).
#
# ST_FOR_DEV_PT_FOR_RELEASE means something like, "in a debug build."
# In a debug build, we do ETI for Scalar=std::complex<float>, in order
# that this capability be tested.  In a release build, we turn off ETI
# for Scalar=std::complex<float>, to reduce build time.
IF (ST_FOR_DEV_PT_FOR_RELEASE AND ${PACKAGE_NAME}_ENABLE_QD AND NOT Tpetra_ENABLE_Kokkos_Refactor)
  # mfh 09 Jan 2015: dd_real does NOT currently work with the Kokkos
  # refactor version of Tpetra, due to a missing Kokkos::atomic_*
  # implementation for this class.  Feel free to remove "AND NOT
  # Tpetra_ENABLE_Kokkos_Refactor" from the above test once that gets
  # fixed.
  GLOBAL_SET(HAVE_TPETRA_INST_DD_REAL_DEFAULT ON)
ELSE()
  GLOBAL_SET(HAVE_TPETRA_INST_DD_REAL_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_DD_REAL
  HAVE_TPETRA_INST_DD_REAL
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes over Scalar = dd_real (requries TPL QD).  This option is ${HAVE_TPETRA_INST_DD_REAL_DEFAULT} by default."
  ${HAVE_TPETRA_INST_DD_REAL_DEFAULT}  
  )
IF(Tpetra_INST_DD_REAL AND Tpetra_ENABLE_Kokkos_Refactor)
  MESSAGE(ERROR "Scalar = dd_real does not currently work with the Kokkos refactor version of Tpetra.  We will fix this at some point.  Please let us know if you need this capability urgently.")
ENDIF()
IF(Tpetra_INST_DD_REAL) 
  LIST(APPEND TpetraCore_ETI_SCALARS "dd_real")
ENDIF()

#
# Set GlobalOrdinal ETI types.  
#
# Unfortunately, we have to add these types to the list of Scalar ETI
# types as well, since some downstream packages use Scalar =
# GlobalOrdinal to implement certain communication patterns.
#

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_INT_LONG
  HAVE_TPETRA_INST_INT_LONG
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes with GlobalOrdinal = long."
  ON
  )
IF(Tpetra_INST_INT_LONG) 
  LIST(APPEND TpetraCore_ETI_GORDS   "long")
  LIST(APPEND TpetraCore_ETI_SCALARS "long")
ENDIF()

# Note that long is 64 bits on some platforms, and 32 bits on other
# platforms.  On the latter, you may use long long instead, which
# C++11 / C99 guarantees to be at least 64 bits long.
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_INT_LONG_LONG
  HAVE_TPETRA_INST_INT_LONG_LONG
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes with GlobalOrdinal = long long."  
  ${Teuchos_ENABLE_LONG_LONG_INT}
  )
IF(Tpetra_INST_INT_LONG_LONG) 
  LIST(APPEND TpetraCore_ETI_GORDS   "long long")
  LIST(APPEND TpetraCore_ETI_SCALARS "long long")
ENDIF()

SET(Tpetra_INST_INT_UNSIGNED
  ${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE}
  CACHE BOOL
  "If explicit instantiation is enabled, instatiate over global ordinal unsigned int."
  )
GLOBAL_SET(HAVE_TPETRA_INST_INT_UNSIGNED ${Tpetra_INST_INT_UNSIGNED})
IF (Tpetra_INST_INT_UNSIGNED)
  LIST(APPEND TpetraCore_ETI_GORDS   "unsigned int")
  LIST(APPEND TpetraCore_ETI_SCALARS "unsigned int")
  # NOTE: The above types are used because they help catch errors with caused
  # when GlobalOrdinal is unsigned and GlobalOrdinal is smaller than
  # global_size_t (using wrong Teuchos::OrdinalTraits<>::invald(). 
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

#
# FIXME (mfh 09 Jan 2015) This shouldn't even be an option.  We should
# examine the Map at run time to determine the best hash function to
# use.
#
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_USE_MURMUR_HASH
  TPETRA_USE_MURMUR_HASH
  "Use the Murmur hash function in Tpetra::Map for global-to-local index lookups, rather than the default hash function.  Murmur hash is more expensive to evaluate, but does a better job of avoiding hash table collisions for Maps that do not look like local permutations of contiguous Maps."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_ENABLE_Kokkos_DistObject
  TPETRA_ENABLE_KOKKOS_DISTOBJECT
  "Use experimental Kokkos::View for DistObject methods, enabling efficient communication for GPU node types."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_ENABLE_MPI_CUDA_RDMA
  TPETRA_ENABLE_MPI_CUDA_RDMA
  "Enable RDMA support for MPI communication between CUDA GPUs.  Only turn this on if you know for sure your MPI library supports it.  Requires also enabling Tpetra_ENABLE_Kokkos_DistObject."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_ENABLE_MMM_Timings
  HAVE_TPETRA_MMM_TIMINGS
  "Enable the timings for the MMM kernels"
  OFF
)

# The "Kokkos refactor" version of Tpetra currently requires that
# several Kokkos subpackages be enabled.  These subpackages should be
# enabled by default (as of 08 Jan 2015).  The following test checks
# that the subpackages have been enabled.

ASSERT_DEFINED(Tpetra_ENABLE_Kokkos_Refactor)
#IF (Tpetra_ENABLE_Kokkos_Refactor AND
#    (NOT ${PACKAGE_NAME}_ENABLE_KokkosCore OR
#     NOT ${PACKAGE_NAME}_ENABLE_TeuchosKokkosCompat OR
#     NOT ${PACKAGE_NAME}_ENABLE_TpetraKernels OR
#     NOT ${PACKAGE_NAME}_ENABLE_KokkosContainers OR
#     NOT ${PACKAGE_NAME}_ENABLE_KokkosAlgorithms OR
#     NOT ${PACKAGE_NAME}_ENABLE_TeuchosKokkosComm))
#  MESSAGE(FATAL_ERROR "Dear user: You enabled the Tpetra_ENABLE_Kokkos_Refactor option, which turns on the Kokkos Refactor version of Tpetra.  However, not all of the required subpackages have been enabled.  These subpackages should be enabled by default, so it could be that you disabled them explicitly.  For now, you must allow them to be enabled, in order to use the Kokkos Refactor version of Tpetra.  Please enable subpackages KokkosCore, KokkosContainers, KokkosAlgorithms, TeuchosKokkosCompat, TeuchosKokkosComm, and TpetraKernels.  Enable these subpackages, if they are not already enabled, by adding the following to your CMake configuration:\n-D Trilinos_ENABLE_KokkosCore:BOOL=ON\n-D Trilinos_ENABLE_TeuchosKokkosCompat:BOOL=ON\n-D Trilinos_ENABLE_TpetraKernels:BOOL=ON\n-D Trilinos_ENABLE_KokkosContainers:BOOL=ON\n-D Trilinos_ENABLE_KokkosAlgorithms:BOOL=ON\n-D Trilinos_ENABLE_TeuchosKokkosComm:BOOL=ON\n\nHere is the current state of these packages:\n\n${PACKAGE_NAME}_ENABLE_KokkosCore = ${${PACKAGE_NAME}_ENABLE_KokkosCore}\n${PACKAGE_NAME}_ENABLE_TeuchosKokkosCompat = ${${PACKAGE_NAME}_ENABLE_TeuchosKokkosCompat}\n${PACKAGE_NAME}_ENABLE_TpetraKernels = ${${PACKAGE_NAME}_ENABLE_TpetraKernels}\n${PACKAGE_NAME}_ENABLE_KokkosContainers = ${${PACKAGE_NAME}_ENABLE_KokkosContainers}\n${PACKAGE_NAME}_ENABLE_KokkosAlgorithms = ${${PACKAGE_NAME}_ENABLE_KokkosAlgorithms}\n${PACKAGE_NAME}_ENABLE_TeuchosKokkosComm = ${${PACKAGE_NAME}_ENABLE_TeuchosKokkosComm}")
#ENDIF()

IF (Tpetra_ENABLE_Kokkos_Refactor)
  SET(Tpetra_Have_Kokkos_Refactor ON CACHE INTERNAL "")
ELSE()
  SET(Tpetra_Have_Kokkos_Refactor OFF CACHE INTERNAL "")
ENDIF()
# Set corresponding value of CPP macro
GLOBAL_SET(TPETRA_HAVE_KOKKOS_REFACTOR ${Tpetra_Have_Kokkos_Refactor})


TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_ENABLE_Kokkos_Refactor_Map
  TPETRA_USE_KOKKOS_REFACTOR_MAP
  "Enable Tpetra::Map specialization for refactor to new Kokkos interface.  If this option is OFF, you may still use the other Kokkos refactor classes with the old Map implementation."
  OFF
  )

IF (Tpetra_ENABLE_Kokkos_Refactor_Map AND NOT Tpetra_Have_Kokkos_Refactor)
  MESSAGE( FATAL_ERROR "If Tpetra_ENABLE_Kokkos_Refactor_Map is ON, then Tpetra's Kokkos refactor code must be enabled as well.  This requires Tpetra_ENABLE_Kokkos_Refactor to be ON, and the following subpackages to be enabled: KokkosCore, KokkosContainers, KokkosAlgorithms, TeuchosKokkosCompat, TeuchosKokkosComm, and TpetraKernels." )
ENDIF ()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_RTI
  HAVE_TPETRACORE_RTI
  "Enable Tpetra's Reduction/Transformation Interface (RTI).  It is OFF by default.  Note that RTI is currently incompatible with the new Kokkos wrapper Nodes."
  OFF
  )
IF (${PACKAGE_NAME}_ENABLE_RTI)
  MESSAGE(STATUS "Enabling Tpetra/RTI")
ELSE()
  MESSAGE(STATUS "Not enabling Tpetra/RTI")
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_TeuchosKokkosCompat AND ${PACKAGE_NAME}_ENABLE_RTI)
  MESSAGE(FATAL_ERROR "You enabled the Reduction/Transformation Interface (RTI) in Tpetra by setting Tpetra_ENABLE_RTI, and you also enabled the TeuchosKokkosCompat subpackage of Kokkos.  The new Kokkos wrapper Nodes in the TeuchosKokkosCompat subpackage are currently incompatible with RTI.  Note that this has nothing to do with the Kokkos refactor version of Tpetra; it's just that the new Kokkos wrapper Nodes are a bit less tolerant.")
ENDIF ()

#
# Add libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(inout)
ADD_SUBDIRECTORY(ext)
IF (Tpetra_ENABLE_RTI)
  ADD_SUBDIRECTORY(rti)
ENDIF()
 
TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

GLOBAL_SET(TpetraCore_ETI_SCALARS ${TpetraCore_ETI_SCALARS})
GLOBAL_SET(TpetraCore_ETI_GORDS   ${TpetraCore_ETI_GORDS})
GLOBAL_SET(TpetraCore_ETI_LORDS   ${TpetraCore_ETI_LORDS})
GLOBAL_SET(TpetraCore_ETI_NODES   ${TpetraCore_ETI_NODES})

#
# Do standard subpackage postprocessing
#
TRIBITS_SUBPACKAGE_POSTPROCESS()
