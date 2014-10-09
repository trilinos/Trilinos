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

TRIBITS_SET_ST_FOR_DEV_MODE(ST_FOR_DEV_PT_FOR_RELEASE)

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

IF (ST_FOR_DEV_PT_FOR_RELEASE AND HAVE_TEUCHOS_BLASFLOAT) 
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

IF (ST_FOR_DEV_PT_FOR_RELEASE AND HAVE_COMPLEX_BLAS AND HAVE_TEUCHOS_BLASFLOAT)
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

IF (ST_FOR_DEV_PT_FOR_RELEASE AND HAVE_COMPLEX_BLAS)
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

IF (ST_FOR_DEV_PT_FOR_RELEASE AND Tpetra_ENABLE_QD)
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

IF (ST_FOR_DEV_PT_FOR_RELEASE AND Tpetra_ENABLE_QD)
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
  "If explicit instantiation is enabled, instantiate over global ordinal long."
  ON
  )
IF(Tpetra_INST_INT_LONG) 
  LIST(APPEND Tpetra_ETI_GORDS   "long")
  LIST(APPEND Tpetra_ETI_SCALARS "long")
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_INT_LONG_LONG
  HAVE_TPETRA_INST_INT_LONG_LONG
  "If explicit instantiation is enabled, instantiate over global ordinal long long."
  ${Teuchos_ENABLE_LONG_LONG_INT}
  )
IF(Tpetra_INST_INT_LONG_LONG) 
  LIST(APPEND Tpetra_ETI_GORDS   "long long")
  LIST(APPEND Tpetra_ETI_SCALARS "long long")
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
  Tpetra_ENABLE_Kokkos_Refactor
  TPETRA_ENABLE_KOKKOS_REFACTOR
  "Enable the \"Kokkos refactor\" version of Tpetra.  In order to enable this, you must also enable the following Kokkos subpackages that do not build by default: Core, Compat, LinAlg, Containers, Algorithms, and MpiComm."
  OFF
  )

# The "Kokkos refactor" version of Tpetra currently requires that
# several Kokkos subpackages be enabled.  These subpackages are not
# enabled by default, because they are not currently (as of August
# 2014) marked Primary Stable.  This will change in the future, but
# for now, you must enable these subpackages explicitly.  The error
# message below explains how to do that.
IF (Tpetra_ENABLE_Kokkos_Refactor AND
    (NOT Tpetra_ENABLE_KokkosCore OR
     NOT Tpetra_ENABLE_KokkosCompat OR
     NOT Tpetra_ENABLE_KokkosLinAlg OR
     NOT Tpetra_ENABLE_KokkosContainers OR
     NOT Tpetra_ENABLE_KokkosAlgorithms OR
     NOT Tpetra_ENABLE_KokkosMpiComm))
  MESSAGE(FATAL_ERROR "Dear user: You enabled the Tpetra_ENABLE_Kokkos_Refactor option, which turns on the Kokkos Refactor version of Tpetra.  However, you did not enable all required packages in order for this to work.  The Kokkos Refactor version of Tpetra depends on some Kokkos subpackages which currently do not build by default.  For now, you must enable them explicitly in order to use the Kokkos Refactor version of Tpetra.  Please enable Kokkos subpackages Core, Compat, LinAlg, Containers, Algorithms, and MpiComm.  Despite the name \"MpiComm\", you may enable that subpackage even if not building Trilinos with MPI support.  Enable these subpackages by adding the following to your CMake configuration:\n-D Trilinos_ENABLE_KokkosAlgorithms:BOOL=ON\n-D Trilinos_ENABLE_KokkosCore:BOOL=ON\n-D Trilinos_ENABLE_KokkosCompat:BOOL=ON\n-D Trilinos_ENABLE_KokkosLinAlg:BOOL=ON\n-D Trilinos_ENABLE_KokkosContainers:BOOL=ON\n-D Trilinos_ENABLE_KokkosAlgorithms:BOOL=ON\n-D Trilinos_ENABLE_KokkosMpiComm:BOOL=ON\nIf you are running Trilinos' check-in test script, you must enable these subpackages explicitly, using arguments to --enable-packages, like this:\n--enable-packages=KokkosCore,KokkosCompat,KokkosContainers,KokkosLinAlg,KokkosAlgorithms,KokkosMpiComm,<YOUR_PACKAGES>\n(Replace <YOUR_PACKAGES> with the other packages you want to enable.)  Also, the check-in test script will not enable these subpackages with the default builds MPI_DEBUG and SERIAL_RELEASE, so you must create another configuration and use it with the --extra-builds=... command-line argument.")
ENDIF()

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
  MESSAGE( FATAL_ERROR "If Tpetra_ENABLE_Kokkos_Refactor_Map is ON, then Tpetra's Kokkos refactor code must be enabled as well.  This requires Tpetra_ENABLE_Kokkos_Refactor to be ON, and the Core, Compat, LinAlg, Containers, and MpiComm subpackages of Kokkos to be enabled." )
ENDIF ()

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
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

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
