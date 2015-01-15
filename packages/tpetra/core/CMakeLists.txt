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

#-----------------------------------------------------------------------
#---------------Explicit Template Instantiation Logic-------------------
#-----------------------------------------------------------------------
SET(TpetraCore_ETI_SCALARS "")
SET(TpetraCore_ETI_GORDS   "")
SET(TpetraCore_ETI_LORDS   "int")
SET(TpetraCore_ETI_NODES   "")

#-----------------------------------------------------------------------
#---------------------Global Ordinals-----------------------------------
#-----------------------------------------------------------------------

# The logic is different depending on whether Tpetra_ENABLE_Kokkos_Refactor is on
IF( Tpetra_ENABLE_Kokkos_Refactor)

# The user gets what he asks for, and if nothing is asked for enable "KokkosSerialWrapperNode"

IF (${${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION})
  IF (Kokkos_ENABLE_Cuda) 
    GLOBAL_SET(HAVE_TPETRA_INST_SERIAL_DEFAULT OFF)
    GLOBAL_SET(HAVE_TPETRA_INST_PTHREAD_DEFAULT OFF)
    GLOBAL_SET(HAVE_TPETRA_INST_OPENMP_DEFAULT OFF)  
    GLOBAL_SET(HAVE_TPETRA_INST_CUDA_DEFAULT ON)  
  ELSE()
    IF (Kokkos_ENABLE_OpenMP) 
      GLOBAL_SET(HAVE_TPETRA_INST_SERIAL_DEFAULT OFF)
      GLOBAL_SET(HAVE_TPETRA_INST_PTHREAD_DEFAULT OFF)
      GLOBAL_SET(HAVE_TPETRA_INST_OPENMP_DEFAULT ON)  
      GLOBAL_SET(HAVE_TPETRA_INST_CUDA_DEFAULT OFF)  
    ELSE()
      #Treat Pthreads differently: by default its off as a node since pthreads is on by default
      IF (Kokkos_ENABLE_Serial)
        GLOBAL_SET(HAVE_TPETRA_INST_SERIAL_DEFAULT ON)
        GLOBAL_SET(HAVE_TPETRA_INST_PTHREAD_DEFAULT OFF)
        GLOBAL_SET(HAVE_TPETRA_INST_OPENMP_DEFAULT OFF)  
        GLOBAL_SET(HAVE_TPETRA_INST_CUDA_DEFAULT OFF)  
      ELSE()
        MESSAGE(STATUS "While pthreads are enabled they will not be explicitly instantiated in Tpetra unless you set Tpetra_INST_PTHREAD=ON")        
        GLOBAL_SET(HAVE_TPETRA_INST_SERIAL_DEFAULT OFF)
        GLOBAL_SET(HAVE_TPETRA_INST_PTHREAD_DEFAULT ON)
        GLOBAL_SET(HAVE_TPETRA_INST_OPENMP_DEFAULT OFF)  
        GLOBAL_SET(HAVE_TPETRA_INST_CUDA_DEFAULT OFF)  
      ENDIF()
    ENDIF()
  ENDIF()  
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_SERIAL
  HAVE_TPETRA_INST_SERIAL
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes over Node = Kokkos::Compat::KokkosSerialWrapperNode.  This option is ${HAVE_TPETRA_INST_SERIAL_DEFAULT} by default."
  ${HAVE_TPETRA_INST_SERIAL_DEFAULT}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_PTHREAD
  HAVE_TPETRA_INST_PTHREAD
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes over Node = Kokkos::Compat::KokkosThreadsWrapperNode.  This option is ${HAVE_TPETRA_INST_PTHREAD_DEFAULT} by default."
  ${HAVE_TPETRA_INST_PTHREAD_DEFAULT}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_OPENMP
  HAVE_TPETRA_INST_OPENMP
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes over Node = Kokkos::Compat::KokkosOpenMPWrapperNode.  This option is ${HAVE_TPETRA_INST_OPENMP_DEFAULT} by default."
  ${HAVE_TPETRA_INST_OPENMP_DEFAULT}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_CUDA
  HAVE_TPETRA_INST_CUDA
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes over Node = Kokkos::Compat::KokkosCudaWrapperNode.  This option is ${HAVE_TPETRA_INST_CUDA_DEFAULT} by default."
  ${HAVE_TPETRA_INST_CUDA_DEFAULT}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_SERIALCLASSIC
  HAVE_TPETRA_INST_SERIALCLASSIC
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes over Node = KokkosClassic::SerialNode.  This option is OFF by default."
  OFF
  )

# Append the instantiation list

IF(Tpetra_INST_SERIALCLASSIC) 
  LIST(APPEND TpetraCore_ETI_NODES "KokkosClassic::SerialNode")
  MESSAGE(STATUS "Enabling Node: KokkosClassic::SerialNode")
ENDIF()
IF(Tpetra_INST_SERIAL) 
  LIST(APPEND TpetraCore_ETI_NODES "Kokkos::Compat::KokkosSerialWrapperNode")
  MESSAGE(STATUS "Enabling Node: Kokkos::Compat::KokkosSerialWrapperNode")
ENDIF()
IF(Tpetra_INST_PTHREAD) 
  LIST(APPEND TpetraCore_ETI_NODES "Kokkos::Compat::KokkosThreadsWrapperNode")
  MESSAGE(STATUS "Enabling Node: Kokkos::Compat::KokkosThreadsWrapperNode")
ENDIF()
IF(Tpetra_INST_OPENMP) 
  LIST(APPEND TpetraCore_ETI_NODES "Kokkos::Compat::KokkosOpenMPWrapperNode")
  MESSAGE(STATUS "Enabling Node: Kokkos::Compat::KokkosOpenMPWrapperNode")
ENDIF()
IF(Tpetra_INST_CUDA) 
  LIST(APPEND TpetraCore_ETI_NODES "Kokkos::Compat::KokkosCudaWrapperNode")
  MESSAGE(STATUS "Enabling Node: Kokkos::Compat::KokkosCudaWrapperNode")
ENDIF()

ELSE()

SET(TpetraCore_ETI_NODES   "${KokkosClassic_ENABLED_NODES}")

ENDIF()


#-----------------------------------------------------------------------
#---------------------Global Ordinals-----------------------------------
#-----------------------------------------------------------------------
#The user gets what he asks for, and if nothing is asked for enable "int"
 
# Default logic
# Note that long is 64 bits on some platforms, and 32 bits on other
# platforms.  On the latter, you may use long long instead, which
# C++11 / C99 guarantees to be at least 64 bits long.

# NOTE: The unsigned int is used because they help catch errors with caused
# when GlobalOrdinal is unsigned and GlobalOrdinal is smaller than
# global_size_t (using wrong Teuchos::OrdinalTraits<>::invald(). 

IF (${${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION})
  IF (Tpetra_INST_INT_LONG_LONG OR HAVE_TEUCHOS_LONG_LONG_INT) 
    GLOBAL_SET(HAVE_TPETRA_INST_INT_LONG_LONG_DEFAULT ON)
    GLOBAL_SET(HAVE_TPETRA_INST_INT_LONG_DEFAULT OFF)
    GLOBAL_SET(HAVE_TPETRA_INST_INT_UNSIGNED_DEFAULT OFF)  
    GLOBAL_SET(HAVE_TPETRA_INST_INT_INT_DEFAULT ON)  
  ELSE()
    IF (Tpetra_INST_INT_LONG) 
      GLOBAL_SET(HAVE_TPETRA_INST_INT_LONG_LONG_DEFAULT OFF)
      GLOBAL_SET(HAVE_TPETRA_INST_INT_LONG_DEFAULT ON)
      GLOBAL_SET(HAVE_TPETRA_INST_INT_UNSIGNED_DEFAULT OFF)  
      GLOBAL_SET(HAVE_TPETRA_INST_INT_INT_DEFAULT ON)  
    ELSE()
      IF (Tpetra_INST_INT_UNSIGNED)
        GLOBAL_SET(HAVE_TPETRA_INST_INT_LONG_LONG_DEFAULT OFF)
        GLOBAL_SET(HAVE_TPETRA_INST_INT_LONG_DEFAULT OFF)
        GLOBAL_SET(HAVE_TPETRA_INST_INT_UNSIGNED_DEFAULT ON)  
        GLOBAL_SET(HAVE_TPETRA_INST_INT_INT_DEFAULT ON)  
      ELSE()
        GLOBAL_SET(HAVE_TPETRA_INST_INT_LONG_LONG_DEFAULT OFF)
        GLOBAL_SET(HAVE_TPETRA_INST_INT_LONG_DEFAULT OFF)
        GLOBAL_SET(HAVE_TPETRA_INST_INT_UNSIGNED_DEFAULT OFF)  
        GLOBAL_SET(HAVE_TPETRA_INST_INT_INT_DEFAULT ON)  
      ENDIF()
    ENDIF()
  ENDIF()  
ENDIF()

# Enable Global Ordinals

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_INT_INT
  HAVE_TPETRA_INST_INT_INT
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes over GlobalOrdinal = int.  This option is ${HAVE_TPETRA_INST_INT_INT_DEFAULT} by default."
  ${HAVE_TPETRA_INST_INT_INT_DEFAULT}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_INT_UNSIGNED
  HAVE_TPETRA_INST_INT_UNSIGNED
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes over GlobalOrdinal = unsigned.  This option is ${HAVE_TPETRA_INST_INT_UNSIGNED_DEFAULT} by default."
  ${HAVE_TPETRA_INST_INT_UNSIGNED_DEFAULT}
  )
 
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_INT_LONG
  HAVE_TPETRA_INST_INT_LONG
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes over GlobalOrdinal = long.  This option is ${HAVE_TPETRA_INST_INT_LONG_DEFAULT} by default."
  ${HAVE_TPETRA_INST_INT_LONG_DEFAULT}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_INT_LONG_LONG
  HAVE_TPETRA_INST_INT_LONG_LONG
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes over GlobalOrdinal = long long.  This option is ${HAVE_TPETRA_INST_INT_LONG_LONG_DEFAULT} by default."
  ${HAVE_TPETRA_INST_INT_LONG_LONG_DEFAULT}
  ) 

# Append to list

IF(Tpetra_INST_INT_INT) 
  LIST(APPEND TpetraCore_ETI_SCALARS "int")
  LIST(APPEND TpetraCore_ETI_GORDS "int")
ENDIF()
IF(Tpetra_INST_INT_UNSIGNED) 
  LIST(APPEND TpetraCore_ETI_SCALARS "unsigned")
  LIST(APPEND TpetraCore_ETI_GORDS "unsigned")
ENDIF()
IF(Tpetra_INST_INT_LONG) 
  LIST(APPEND TpetraCore_ETI_SCALARS "long")
  LIST(APPEND TpetraCore_ETI_GORDS "long")
ENDIF()
IF(Tpetra_INST_INT_LONG_LONG) 
  LIST(APPEND TpetraCore_ETI_SCALARS "long long")
  LIST(APPEND TpetraCore_ETI_GORDS "long long")
ENDIF()

#-----------------------------------------------------------------------
#---------------------Scalars-------------------------------------------
#-----------------------------------------------------------------------
# The user gets what he asks for, and if nothing is asked for enable "double"
# Double is always enabled since tests explicitly depend on it.

# Default Logic

IF (${${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION})
  IF (Tpetra_INST_FLOAT) 
    GLOBAL_SET(HAVE_TPETRA_INST_FLOAT_DEFAULT ON)
    GLOBAL_SET(HAVE_TPETRA_INST_DOUBLE_DEFAULT ON)  
  ELSE()
    GLOBAL_SET(HAVE_TPETRA_INST_FLOAT_DEFAULT OFF)
    GLOBAL_SET(HAVE_TPETRA_INST_DOUBLE_DEFAULT ON)
  ENDIF()  
ENDIF()

# Enable (Floaty) Scalars

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_DOUBLE
  HAVE_TPETRA_INST_DOUBLE
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes over Scalar = double.  This option is ${HAVE_TPETRA_INST_DOUBLE_DEFAULT} by default."
  ${HAVE_TPETRA_INST_DOUBLE_DEFAULT}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_FLOAT
  HAVE_TPETRA_INST_FLOAT
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes over Scalar = float.  This option is ${HAVE_TPETRA_INST_FLOAT_DEFAULT} by default."
  ${HAVE_TPETRA_INST_FLOAT_DEFAULT}
  )

# Enable complex
 
IF (HAVE_COMPLEX_BLAS)
  IF (Tpetra_INST_DOUBLE)
    GLOBAL_SET(HAVE_Tpetra_INST_COMPLEX_DOUBLE_DEFAULT ON)
  ELSE()
    GLOBAL_SET(HAVE_Tpetra_INST_COMPLEX_DOUBLE_DEFAULT OFF)
  ENDIF()
  IF (Tpetra_INST_FLOAT)
    GLOBAL_SET(HAVE_Tpetra_INST_COMPLEX_FLOAT_DEFAULT ON)
  ELSE()
    GLOBAL_SET(HAVE_Tpetra_INST_COMPLEX_FLOAT_DEFAULT OFF)
  ENDIF()
ELSE()
  GLOBAL_SET(HAVE_Tpetra_INST_COMPLEX_DOUBLE_DEFAULT OFF)
  GLOBAL_SET(HAVE_Tpetra_INST_COMPLEX_FLOAT_DEFAULT OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_COMPLEX_FLOAT
  HAVE_TPETRA_INST_COMPLEX_FLOAT
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes with Scalar = std::complex<float>.  This option is ${HAVE_TPETRA_INST_COMPLEX_FLOAT_DEFAULT} by default."
  ${HAVE_Tpetra_INST_COMPLEX_FLOAT_DEFAULT}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_COMPLEX_DOUBLE
  HAVE_TPETRA_INST_COMPLEX_DOUBLE
  "If explicit template instantiation (ETI) is enabled, instantiate Tpetra classes with Scalar = std::complex<double>.  This option is ${HAVE_TPETRA_INST_COMPLEX_DOUBLE_DEFAULT} by default."
  ${HAVE_Tpetra_INST_COMPLEX_DOUBLE_DEFAULT}
  )

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
  
# Check for support in BLAS
  
IF(Tpetra_INST_FLOAT AND NOT HAVE_TEUCHOS_BLASFLOAT)
  MESSAGE(SEND_ERROR "You may not instantiate Tpetra classes with Scalar = float, if the BLAS library you are using does not support float (S) arithmetic.")
ENDIF()
IF(Tpetra_INST_COMPLEX_FLOAT AND (NOT HAVE_COMPLEX_BLAS OR NOT HAVE_TEUCHOS_BLASFLOAT))
  MESSAGE(SEND_ERROR "You may not instantiate Tpetra classes with Scalar = std::complex<float>, if the BLAS library you are using does not support complex (C) or float (S) arithmetic.")
ENDIF()
IF(Tpetra_INST_COMPLEX_DOUBLE AND NOT HAVE_COMPLEX_BLAS)
  MESSAGE(SEND_ERROR "You may not instantiate Tpetra classes with Scalar = std::complex<double>, if the BLAS library you are using does not support double-complex (Z) arithmetic.")
ENDIF()
IF(Tpetra_INST_QD_REAL AND Tpetra_ENABLE_Kokkos_Refactor)
  MESSAGE(SEND_ERROR "Scalar = qd_real does not currently work with the Kokkos refactor version of Tpetra.  We will fix this at some point.  Please let us know if you need this capability urgently.")
ENDIF()
IF(Tpetra_INST_DD_REAL AND Tpetra_ENABLE_Kokkos_Refactor)
  MESSAGE(SEND_ERROR "Scalar = dd_real does not currently work with the Kokkos refactor version of Tpetra.  We will fix this at some point.  Please let us know if you need this capability urgently.")
ENDIF()

# Add to list 

IF(Tpetra_INST_FLOAT) 
  LIST(APPEND TpetraCore_ETI_SCALARS "float")
ENDIF()
IF(Tpetra_INST_DOUBLE) 
  LIST(APPEND TpetraCore_ETI_SCALARS "double")
ENDIF()
IF(Tpetra_INST_COMPLEX_FLOAT) 
  LIST(APPEND TpetraCore_ETI_SCALARS "std::complex<float>")
ENDIF()
IF(Tpetra_INST_COMPLEX_DOUBLE) 
  LIST(APPEND TpetraCore_ETI_SCALARS "std::complex<double>")
ENDIF()
IF(Tpetra_INST_QD_REAL) 
  LIST(APPEND TpetraCore_ETI_SCALARS "qd_real")
ENDIF()
IF(Tpetra_INST_DD_REAL) 
  LIST(APPEND TpetraCore_ETI_SCALARS "dd_real")
ENDIF()

IF (${${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION})
  IF(Tpetra_INST_FLOAT OR Tpetra_INST_DOUBLE OR Tpetra_INST_COMPLEX_FLOAT OR Tpetra_INST_COMPLEX_DOUBLE OR Tpetra_INST_QD_REAL OR Tpetra_INST_DD_REAL)
    GLOBAL_SET(HAVE_TPETRA_SCALAR ON)
  ELSE()
    MESSAGE(FATAL_ERROR "Tpetra: No Scalar types are enabled. One Scalar type needs to be enabled to use ETI.") 
  ENDIF()
  IF(Tpetra_INST_INT_INT OR Tpetra_INST_INT_UNSIGNED OR Tpetra_INST_INT_LONG OR Tpetra_INST_INT_LONG_LONG)
    GLOBAL_SET(HAVE_TPETRA_GLOBALORDINAL ON)
  ELSE()
    MESSAGE(FATAL_ERROR "Tpetra: No GlobalOrdinal type is enabled. One GlobalOrdinal type needs to be enabled to use ETI.")   
  ENDIF() 
ENDIF()

ASSERT_DEFINED(HAVE_TPETRA_SCALAR)
ASSERT_DEFINED(HAVE_TPETRA_GLOBALORDINAL)

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

IF (Tpetra_ENABLE_Kokkos_Refactor AND ${PACKAGE_NAME}_ENABLE_RTI)
  MESSAGE(FATAL_ERROR "You enabled the Reduction/Transformation Interface (RTI) in Tpetra by setting Tpetra_ENABLE_RTI, yet the Kokkos refactor (new) version of Tpetra is enabled.  We have not yet certified RTI for use with the new version of Tpetra.  If you really need RTI right now, you may bring back the old version of Tpetra temporarily by setting Tpetra_ENABLE_Kokkos_Refactor:BOOL=OFF.  Please note that the old version of Tpetra will go away completely with the 12.0 release of Trilinos.")
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
