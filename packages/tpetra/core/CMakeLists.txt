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

ASSERT_DEFINED(Kokkos_ENABLE_Cuda)
ASSERT_DEFINED(Kokkos_ENABLE_Cuda_UVM)
IF (Kokkos_ENABLE_Cuda AND NOT Kokkos_ENABLE_Cuda_UVM)
  MESSAGE (FATAL_ERROR "If CUDA is enabled in Kokkos, Tpetra requires that Kokkos' UVM support be enabled.  You may do this by setting the CMake option Kokkos_ENABLE_Cuda_UVM:BOOL=ON (WARNING: IT IS CASE SENSITIVE!) and running CMake again.\n\nDetails for developers: UVM stands for \"Unified Virtual Memory\".  It lets code running on the host processor access GPU memory.  There is a difference between CUDA's support for UVM, and Kokkos' support for UVM.  Versions of CUDA >= 6 have UVM support built in by default.  Kokkos always supports this.  In particular, Kokkos always has a memory space for UVM allocations, called Kokkos::CudaUVMSpace.  \"Turning on UVM support in Kokkos\" means two things:\n\n1. Kokkos::Cuda::memory_space (the CUDA execution space's default memory space) is Kokkos::CudaUVMSpace, rather than Kokkos::CudaSpace.\n\n2. Kokkos::DualView<T, Kokkos::Cuda> only uses a single allocation, in the Kokkos::CudaUVMSpace memory space.")
ENDIF ()

SET(Tpetra_MACHINE_XML_FILE_DIR 
    ${CMAKE_CURRENT_SOURCE_DIR}/machine_files
    CACHE INTERNAL "")

# Set variables which were implicitly set when the Kokkos packages were optional dependencies
GLOBAL_SET(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT ON)
GLOBAL_SET(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMM ON)
GLOBAL_SET(HAVE_TPETRACORE_KOKKOSCORE ON)
GLOBAL_SET(HAVE_TPETRACORE_KOKKOSALGORITHMS ON)
GLOBAL_SET(HAVE_TPETRACORE_KOKKOSCONTAINERS ON)

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
  ${PACKAGE_NAME}_ENABLE_BUGTESTS
  HAVE_TPETRA_BUGTESTS
  "Enabling testing of previous Tpetra bugs."
  ${${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE}
  )

######################################################################
# Explicit template instantiation (ETI) and test instantiation logic
######################################################################

# For Scalar types, Node types, and (LocalOrdinal, GlobalOrdinal) type
# pairs, Tpetra defines macros which say
#
#   1. Whether Tpetra instantiates (if ETI is enabled) and/or tests
#      its objects for that type (or type pair)
#   2. Whether the type (or type pair) is available
#
# If ETI is enabled, both are identical.  If not, the settings of the
# latter depend on both Kokkos settings and Trilinos' settings.
#
# The availability macros have the following pattern: HAVE_TPETRA_[TYPE]
# The instantiation macros have the pattern:          HAVE_TPETRA_INST_[TYPE]
# The CMake options have the pattern:                 Tpetra_INST_[TYPE]
#
# The availability and instantiation macros are THE SAME.
#
# The ordinal types are: INT_INT, INT_UNSIGNED, INT_LONG, INT_LONG_LONG, INT_UNSIGNED_LONG
# The node types are:    SERIALCLASSIC, SERIAL, OPENMP, PTHREAD, CUDA
# The scalar types are:  DOUBLE, FLOAT, COMPLEX_DOUBLE, COMPLEX_FLOAT
#
# A downstream package or application should use the availability
# macros to check whether something can be used/instantiated.
#
# Note on (1) above: Tpetra also adds enabled GlobalOrdinal types
# internally to the Scalar types for instantiation of MultiVector.
#
# Note on (2) above: some of the types can only be enabled if certain
# requirements are fulfilled.

SET(TpetraCore_ETI_SCALARS "")
SET(TpetraCore_ETI_SCALARS_NO_ORDS "") # exclude all ordinal types (GlobalOrdinal and int)
SET(TpetraCore_ETI_GORDS   "")
SET(TpetraCore_ETI_LORDS   "int")
SET(TpetraCore_ETI_NODES   "")

# ============================================================
# Node types
# ============================================================

# Decide what Node types to enable.  If ETI is ON, "enabling a Node
# type" means that Tpetra objects with Node as a template parameter
# will
# 
#   1. get instantiated explicitly for that Node type, and
#   2. get tested for that Node type, if their test is also
#      templated on Node type.
#
# If ETI is OFF, #1 no longer holds, but #2 still holds.
#
# We enable exactly one Node type by default, whether ETI is ON or
# OFF.  This keeps build times for tests down in the non-ETI case, and
# library sizes small in the ETI case.

GLOBAL_SET(HAVE_TPETRA_INST_SERIAL_DEFAULT OFF)
GLOBAL_SET(HAVE_TPETRA_INST_PTHREAD_DEFAULT OFF)
GLOBAL_SET(HAVE_TPETRA_INST_OPENMP_DEFAULT OFF)  
GLOBAL_SET(HAVE_TPETRA_INST_CUDA_DEFAULT OFF)  

IF (Kokkos_ENABLE_Cuda) 
  # If you bothered to build with CUDA enabled, then we assume you
  # only want CUDA.
  GLOBAL_SET(HAVE_TPETRA_INST_CUDA_DEFAULT ON)  
ELSE()
  IF (Kokkos_ENABLE_OpenMP) 
    # It takes effort (setting a nondefault compiler flag) to turn
    # on OpenMP, so we assume that if OpenMP is ON, then you must
    # want it as the Node type.
    GLOBAL_SET(HAVE_TPETRA_INST_OPENMP_DEFAULT ON)  
  ELSE()
    IF (Kokkos_ENABLE_Serial)
      GLOBAL_SET(HAVE_TPETRA_INST_SERIAL_DEFAULT ON)
    ELSE()
      IF (Kokkos_ENABLE_Pthread)
        GLOBAL_SET(HAVE_TPETRA_INST_PTHREAD_DEFAULT ON)
      ELSE()
        MESSAGE(FATAL_ERROR "Tpetra: No Kokkos execution space is enabled.  This message should only come up if you explicitly disabled all the Kokkos execution spaces.  If you don't want Tpetra to use threads, try setting Kokkos_ENABLE_Serial=ON.")
      ENDIF()  # Kokkos_ENABLE_Pthread
    ENDIF() # Kokkos_ENABLE_Serial
  ENDIF() # Kokkos_ENABLE_OpenMP
ENDIF() # Kokkos_ENABLE_Cuda

# Don't turn on both Serial and Threads by default.  Prefer Serial,
# because Threads is a last resort.  Of course, users may still
# enable Threads explicitly, by setting TPETRA_INST_PTHREAD=ON.
IF (HAVE_TPETRA_INST_SERIAL_DEFAULT AND HAVE_TPETRA_INST_PTHREAD_DEFAULT)
  GLOBAL_SET(HAVE_TPETRA_INST_PTHREAD_DEFAULT OFF)
ENDIF ()

# The Kokkos refactor version of Tpetra is ON, so the "classic" Node
# type must be OFF.  (We will enforce this below.)
GLOBAL_SET(HAVE_TPETRA_INST_SERIALCLASSIC_DEFAULT OFF)

# MESSAGE(STATUS "Tpetra: Default Node enables:")
# MESSAGE(STATUS "   KokkosSerialWrapperNode:  ${HAVE_TPETRA_INST_SERIAL_DEFAULT}")
# MESSAGE(STATUS "   KokkosThreadsWrapperNode: ${HAVE_TPETRA_INST_PTHREAD_DEFAULT}")
# MESSAGE(STATUS "   KokkosOpenMPWrapperNode:  ${HAVE_TPETRA_INST_OPENMP_DEFAULT}")
# MESSAGE(STATUS "   KokkosCudaWrapperNode:    ${HAVE_TPETRA_INST_CUDA_DEFAULT}")

# Set up all the public "Tpetra_INST_*" CMake options for enabling
# Node types.  Each enables both instantations and tests in the ETI ON
# case, and only tests in the ETI OFF case.  Nevertheless, the options
# exist whether ETI is ON or OFF.
#
# While setting up the options, check that the default Node type is
# actually enabled.  Print out a helpful STATUS message if a Kokkos
# execution space is ON, but ETI for its corresponding Node type is
# OFF.  This may help users diagnose link errors.

# Kokkos::Serial (Kokkos::Compat::KokkosSerialWrapperNode)
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_SERIAL
  HAVE_TPETRA_INST_SERIAL
  "Instantiate and/or test Tpetra classes over Node = Kokkos::Compat::KokkosSerialWrapperNode.  If ETI is OFF, enable tests for that Node type.  This option is ${HAVE_TPETRA_INST_SERIAL_DEFAULT} by default."
  ${HAVE_TPETRA_INST_SERIAL_DEFAULT}
  )
GLOBAL_SET(HAVE_TPETRA_SERIAL ${Tpetra_INST_SERIAL})
IF (Tpetra_INST_SERIAL)
  IF (NOT Kokkos_ENABLE_Serial)
    MESSAGE(FATAL_ERROR "Tpetra: The Kokkos::Serial execution space is disabled, but you enabled the corresponding Tpetra Node type by setting Tpetra_INST_SERIAL=ON.  If you want to enable instantiation and use of Kokkos::Serial in Tpetra, you must enable the Kokkos::Serial execution space by setting Kokkos_ENABLE_Serial=ON.")
  ENDIF ()
ELSE () # NOT Tpetra_INST_SERIAL
  IF (HAVE_KOKKOSCLASSIC_DEFAULTNODE_SERIALWRAPPERNODE)
    MESSAGE(FATAL_ERROR "Tpetra: Node = Kokkos::Compat::KokkosSerialWrapperNode is disabled (since Tpetra_INST_SERIAL=OFF), but you set it as the default Node type.  Try setting the CMake options Kokkos_ENABLE_Serial:BOOL=ON and Tpetra_INST_SERIAL:BOOL=ON.")
  ENDIF ()
  IF (Kokkos_ENABLE_Serial)
    MESSAGE(STATUS "NOTE: Kokkos::Serial is ON (the CMake option Kokkos_ENABLE_Serial is ON), but the corresponding Tpetra Node type is disabled.  If you want to enable instantiation and use of Kokkos::Serial in Tpetra, please also set the CMake option Tpetra_INST_SERIAL:BOOL=ON.  If you use the Kokkos::Serial Node type in Tpetra without doing this, you will get link errors!")
  ENDIF ()
ENDIF () # Tpetra_INST_SERIAL

# Kokkos::OpenMP (Kokkos::Compat::KokkosOpenMPWrapperNode)
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_OPENMP
  HAVE_TPETRA_INST_OPENMP
  "Instantiate and/or test Tpetra classes over Node = Kokkos::Compat::KokkosOpenMPWrapperNode.  If ETI is OFF, enable tests for that Node type.  This option is ${HAVE_TPETRA_INST_OPENMP_DEFAULT} by default."
  ${HAVE_TPETRA_INST_OPENMP_DEFAULT}
  )
GLOBAL_SET(HAVE_TPETRA_OPENMP ${Tpetra_INST_OPENMP})
IF (Tpetra_INST_OPENMP)
  IF (NOT Kokkos_ENABLE_OpenMP)
    MESSAGE(FATAL_ERROR "Tpetra: The Kokkos::OpenMP execution space is disabled, but you enabled the corresponding Tpetra Node type by setting Tpetra_INST_OPENMP=ON.  If you want to enable instantiation and use of Kokkos::OpenMP in Tpetra, you must enable the Kokkos::OpenMP execution space by setting Kokkos_ENABLE_OpenMP=ON.")
  ENDIF ()
ELSE () # NOT Tpetra_INST_OPENMP
  IF (HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPWRAPPERNODE)
    MESSAGE(FATAL_ERROR "Tpetra: Node = Kokkos::Compat::KokkosOpenMPWrapperNode is disabled (since Tpetra_INST_SERIAL=OFF), but you set it as the default Node type.  Try setting the CMake option Tpetra_INST_SERIAL:BOOL=ON.  If that does not work, also add Trilinos_ENABLE_OpenMP:BOOL=ON and possibly also Kokkos_ENABLE_OpenMP:BOOL=ON.")
  ENDIF()
  IF (Kokkos_ENABLE_OpenMP)
    MESSAGE(STATUS "NOTE: Kokkos::OpenMP is ON (the CMake option Kokkos_ENABLE_OpenMP is ON), but the corresponding Tpetra Node type is disabled.  If you want to enable instantiation and use of Kokkos::OpenMP in Tpetra, please also set the CMake option Tpetra_INST_OPENMP:BOOL=ON.  If you use the Kokkos::OpenMP Node type in Tpetra without doing this, you will get link errors!")
  ENDIF ()
ENDIF ()

# Kokkos::Threads (Kokkos::Compat::KokkosThreadsWrapperNode)
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_PTHREAD
  HAVE_TPETRA_INST_PTHREAD
  "Instantiate and/or test Tpetra classes over Node = Kokkos::Compat::KokkosThreadsWrapperNode.  If ETI is OFF, enable tests for that Node type.  This option is ${HAVE_TPETRA_INST_PTHREAD_DEFAULT} by default."
  ${HAVE_TPETRA_INST_PTHREAD_DEFAULT}
  )
GLOBAL_SET(HAVE_TPETRA_PTHREAD ${Tpetra_INST_PTHREAD})
IF (Tpetra_INST_PTHREAD)
  IF (NOT Kokkos_ENABLE_Pthread)
    MESSAGE(FATAL_ERROR "Tpetra: The Kokkos::Threads execution space is disabled, but you enabled the corresponding Tpetra Node type by setting Tpetra_INST_PTHREAD=ON [sic; it is really PTHREAD, not PTHREADS].  If you want to enable instantiation and use of Kokkos::Threads in Tpetra, you must enable the Kokkos::Threads execution space by setting Kokkos_ENABLE_Pthread=ON [again sic].")
  ENDIF ()
  IF (Tpetra_INST_OPENMP)
    MESSAGE(WARNING "We do NOT recommend enabling both the Kokkos::Serial and Kokkos::Threads Node types in Tpetra in the same build!  The threads that the Kokkos::Threads execution space creates will interfere with OpenMP's threads, causing performance problems.  You are seeing this warning because both CMake options Tpetra_INST_OPENMP and Tpetra_INST_PTHREAD are ON, which is not the default.")
  ENDIF ()
ELSE () # NOT Tpetra_INST_PTHREAD
  IF (HAVE_KOKKOSCLASSIC_DEFAULTNODE_THREADSWRAPPERNODE) 
    MESSAGE(FATAL_ERROR "Tpetra: Node = Kokkos::Compat::KokkosThreadsWrapperNode is disabled (since Tpetra_INST_PTHREAD=OFF), but you set it as the default Node type.  Try setting the CMake options Kokkos_ENABLE_Pthread:BOOL=ON and Tpetra_INST_PTHREAD:BOOL=ON.  (Yes, it is 'Pthread', not 'Pthreads'.  Sorry, that wasn't my idea.)")
  ENDIF()
  IF (Kokkos_ENABLE_Pthread)
    MESSAGE(STATUS "NOTE: Kokkos::Threads is ON (the CMake option Kokkos_ENABLE_Pthread is ON), but the corresponding Tpetra Node type is disabled.  Note that unlike with other Kokkos execution spaces, Tpetra disables Kokkos::Threads by default, as long as some other Kokkos execution space is enabled.  If you want to enable instantiation and use of Kokkos::Threads in Tpetra, please also set the CMake option Tpetra_INST_PTHREAD:BOOL=ON.  If you use the Kokkos::Threads version of Tpetra without doing this, you will get link errors!  We do NOT recommend using Kokkos::Threads if Kokkos::OpenMP is enabled.")
  ENDIF ()
ENDIF () # Tpetra_INST_PTHREAD

# Kokkos::Cuda (Kokkos::Compat::KokkosCudaWrapperNode)
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_CUDA
  HAVE_TPETRA_INST_CUDA
  "Instantiate and/or test Tpetra classes over Node = Kokkos::Compat::KokkosCudaWrapperNode.  If ETI is OFF, enable tests for that Node type.  This option is ${HAVE_TPETRA_INST_CUDA_DEFAULT} by default."
  ${HAVE_TPETRA_INST_CUDA_DEFAULT}
  )
GLOBAL_SET(HAVE_TPETRA_CUDA ${Tpetra_INST_CUDA})
IF (Tpetra_INST_CUDA)
  IF (NOT Kokkos_ENABLE_Cuda)
    MESSAGE(FATAL_ERROR "Tpetra: The Kokkos::Cuda execution space is disabled, but you enabled the corresponding Tpetra Node type by setting Tpetra_INST_CUDA=ON.  If you want to enable instantiation and use of Kokkos::Cuda in Tpetra, you must enable the Kokkos::Cuda execution space by setting Kokkos_ENABLE_Cuda=ON.")
  ENDIF ()
ELSE () # NOT Tpetra_INST_CUDA
  IF (HAVE_KOKKOSCLASSIC_DEFAULTNODE_CUDAWRAPPERNODE) 
    MESSAGE(FATAL_ERROR "Tpetra: Node = Kokkos::Compat::KokkosCudaWrapperNode is disabled (since Tpetra_INST_CUDA=OFF), but you set it as the default Node type.  Try setting the CMake options Kokkos_ENABLE_Cuda:BOOL=ON and Tpetra_INST_CUDA:BOOL=ON.  If you are building with a CUDA-capable compiler and Kokkos can detect that, then you are unlikely to see this message, since both Kokkos and Tpetra enable CUDA support by default in that case.")
  ENDIF()
  IF (Kokkos_ENABLE_Cuda)
    MESSAGE(STATUS "NOTE: Kokkos::Cuda is ON (the CMake option Kokkos_ENABLE_Cuda is ON), but the corresponding Tpetra Node type is disabled.  If you want to enable instantiation and use of Kokkos::Cuda in Tpetra, please also set the CMake option Tpetra_INST_CUDA:BOOL=ON.  If you use the Kokkos::Cuda version of Tpetra without doing this, you will get link errors!")
  ENDIF ()
ENDIF () # Tpetra_INST_CUDA

#
# We've survived the gauntlet of tests above, so it's now safe to
# construct the list of enabled Node types.
#
IF(Tpetra_INST_SERIAL) 
  LIST(APPEND TpetraCore_ETI_NODES "Kokkos::Compat::KokkosSerialWrapperNode")
ENDIF()
IF(Tpetra_INST_PTHREAD) 
  LIST(APPEND TpetraCore_ETI_NODES "Kokkos::Compat::KokkosThreadsWrapperNode")
ENDIF()
IF(Tpetra_INST_OPENMP) 
  LIST(APPEND TpetraCore_ETI_NODES "Kokkos::Compat::KokkosOpenMPWrapperNode")
ENDIF()
IF(Tpetra_INST_CUDA) 
  LIST(APPEND TpetraCore_ETI_NODES "Kokkos::Compat::KokkosCudaWrapperNode")
ENDIF()

# Tell users what Nodes are enabled.
MESSAGE(STATUS "Tpetra Node availability (ON means available): ")
MESSAGE(STATUS "   KokkosSerialWrapperNode:  ${HAVE_TPETRA_SERIAL}")
MESSAGE(STATUS "   KokkosThreadsWrapperNode: ${HAVE_TPETRA_PTHREAD}")
MESSAGE(STATUS "   KokkosOpenMPWrapperNode:  ${HAVE_TPETRA_OPENMP}")
MESSAGE(STATUS "   KokkosCudaWrapperNode:    ${HAVE_TPETRA_CUDA}")

# ============================================================
# GlobalOrdinal types
# ============================================================

# Decide what GlobalOrdinal types to enable by default.  If ETI is ON,
# "enabling a GlobalOrdinal type" means that Tpetra objects with
# GlobalOrdinal as a template parameter will
# 
#   1. get instantiated explicitly for that GlobalOrdinal type, and
#   2. get tested for that GlobalOrdinal type, if their test is also
#      templated on GlobalOrdinal type.
#
# If ETI is OFF, #1 no longer holds, but #2 still holds.
#
# Currently, we always enable int, because MueLu requires it for
# examples and tests (this is Bug 6358).  Furthermore, if Teuchos
# enables long long support (which it now does by default with C++11)
# and CUDA is OFF (CUDA doesn't seem to like GO=long long), we also
# enable long long (which C++11 / C99 guarantees to be >= 64 bits).
# We currently only enable other types if the user asks for them.
#
# Note that long is 64 bits on some platforms, and 32 bits on other
# platforms.  On the latter, you may use long long instead, which
# C++11 / C99 guarantees to be at least 64 bits long.
#
# Tpetra has support for GlobalOrdial = unsigned and unsigned long,
# but only reluctantly.  Both are OFF by default.  GO = unsigned helps
# catch potential errors when GlobalOrdinal is unsigned and smaller
# than global_size_t (e.g., using the wrong
# Teuchos::OrdinalTraits<>::invalid().

# GO = int is enabled by default, but only because of Bug 6358.
GLOBAL_SET(HAVE_TPETRA_INST_INT_INT_DEFAULT ON)

# GO = long long is enabled by default if and only if CUDA is OFF.
# However, if the user enabled GO = long EXPLICITLY, we disable GO =
# long long by default, because we don't want to enable multiple
# "long" GO types by default.
IF (Kokkos_ENABLE_Cuda)
  GLOBAL_SET(HAVE_TPETRA_INST_INT_LONG_LONG_DEFAULT OFF)
ELSE ()
  # Tpetra currently needs Teuchos to support long long.
  # This may change in the future, though.
  IF (DEFINED Teuchos_ENABLE_LONG_LONG_INT)
    GLOBAL_SET(HAVE_TPETRA_INST_INT_LONG_LONG_DEFAULT ${Teuchos_ENABLE_LONG_LONG_INT})
  ELSE ()
    GLOBAL_SET(HAVE_TPETRA_INST_INT_LONG_LONG_DEFAULT ON)
  ENDIF ()
ENDIF ()
IF (DEFINED Tpetra_INST_INT_LONG) # if the user set it explicitly (either ON or OFF)
  IF (Tpetra_INST_INT_LONG) # if the user set it to ON explicitly
    # If the user explicitly set Tpetra_INST_INT_LONG=ON, then turn
    # off GO = long long by default.  We have to go through the
    # two-level IF statement because the Tpetra_INST_INT_LONG option
    # doesn't exist yet; we have not yet encountered the
    # TRIBITS_ADD_OPTION_AND_DEFINE statement (below) that defines it
    # and gives it a default value.  (Note the circular dependency
    # between the default value of Tpetra_INST_INT_LONG and
    # Tpetra_INST_INT_LONG_LONG.  This is deliberate.  The point is
    # that only one of them is defined by default at a time.
    GLOBAL_SET(HAVE_TPETRA_INST_INT_LONG_LONG_DEFAULT OFF)
  ENDIF ()
ENDIF ()

# If GO = long long is OFF by default, turn ON GO = long by default.
# If GO = long long is ON by default, turn OFF GO = long by default.
# The point is to have EXACTLY one potentially 64-bit GO type enabled
# by default.  I say "potentially" because the C++ standard only
# guarantees long to be at least 32 bits.
IF (HAVE_TPETRA_INST_INT_LONG_LONG_DEFAULT)
  GLOBAL_SET(HAVE_TPETRA_INST_INT_LONG_DEFAULT OFF)
ELSE ()
  GLOBAL_SET(HAVE_TPETRA_INST_INT_LONG_DEFAULT ON)
ENDIF ()

# GO = unsigned is NOT enabled by default.
GLOBAL_SET(HAVE_TPETRA_INST_INT_UNSIGNED_DEFAULT OFF)
# GO = unsigned long is NOT enabled by default.
GLOBAL_SET(HAVE_TPETRA_INST_INT_UNSIGNED_LONG_DEFAULT OFF)

#
# Define options for enabling various GlobalOrdinal types
#

TRIBITS_ADD_OPTION_AND_DEFINE(  
  Tpetra_INST_INT_INT
  HAVE_TPETRA_INST_INT_INT
  "Enable GlobalOrdinal = int in Tpetra.  This option is ${HAVE_TPETRA_INST_INT_INT_DEFAULT} by default."
  ${HAVE_TPETRA_INST_INT_INT_DEFAULT}
  )
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_INT_UNSIGNED
  HAVE_TPETRA_INST_INT_UNSIGNED
  "Enable GlobalOrdinal = unsigned in Tpetra.  This option is ${HAVE_TPETRA_INST_INT_UNSIGNED_DEFAULT} by default.  Please don't set this option.  We only maintain it for backwards compatibility."
  ${HAVE_TPETRA_INST_INT_UNSIGNED_DEFAULT}
  )
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_INT_LONG_LONG
  HAVE_TPETRA_INST_INT_LONG_LONG
  "Enable GlobalOrdinal = long long in Tpetra.  This option is ${HAVE_TPETRA_INST_INT_LONG_LONG_DEFAULT} by default."
  ${HAVE_TPETRA_INST_INT_LONG_LONG_DEFAULT}
  ) 
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_INT_LONG
  HAVE_TPETRA_INST_INT_LONG
  "Enable GlobalOrdinal = long in Tpetra.  This option is ${HAVE_TPETRA_INST_INT_LONG_DEFAULT} by default."
  ${HAVE_TPETRA_INST_INT_LONG_DEFAULT}
  )
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_INT_UNSIGNED_LONG
  HAVE_TPETRA_INST_INT_UNSIGNED_LONG
  "Enable GlobalOrdinal = unsigned long in Tpetra.  This option is ${HAVE_TPETRA_INST_INT_UNSIGNED_LONG_DEFAULT} by default.  We do NOT recommend using this type.  We are NOT responsible for ensuring that it builds or passes tests.  If you want a 64-bit type, we prefer that you use long long (Tpetra_INST_INT_LONG_LONG=ON), or long (Tpetra_INST_INT_LONG=ON) if you know that it is 64 bits.  If you are interested in using unsigned long, Roger Pawlowski is the responsible party; please contact him."
  ${HAVE_TPETRA_INST_INT_UNSIGNED_LONG_DEFAULT}
  )

IF (DEFINED Teuchos_ENABLE_LONG_LONG_INT)
  IF (Tpetra_INST_INT_LONG_LONG AND NOT Teuchos_ENABLE_LONG_LONG_INT)
    MESSAGE (FATAL_ERROR "Tpetra: GlobalOrdinal = long long is enabled in Tpetra but is not enabled in Teuchos. Set Teuchos_ENABLE_LONG_LONG_INT=ON or set Tpetra_INST_INT_LONG_LONG=OFF.  Note that Tpetra requires C++11, and Teuchos_ENABLE_LONG_LONG_INT is ON by default if the compiler supports C++11.  Thus, if you are seeing this message, you probably set Teuchos_ENABLE_LONG_LONG_INT=OFF explicitly.  Please don't do that.  If you don't want to enable GlobalOrdinal = long long in Tpetra, set Tpetra_INST_INT_LONG_LONG=OFF and leave Teuchos alone.")
  ENDIF ()
ENDIF ()

# Warn on use of options that we don't like.
IF (Tpetra_INST_INT_UNSIGNED)
  MESSAGE (WARNING "Tpetra: You set Tpetra_INST_INT_UNSIGNED.  This means that you want to use GlobalOrdinal = unsigned with Tpetra objects.  Please don't do that.  I only maintain this option for backwards compatibility.  While I have gone through some effort to make unsigned GlobalOrdinal types work in Tpetra, they are not tested and thus I don't want to promise that they still work.")
ENDIF ()
IF (Tpetra_INST_INT_UNSIGNED_LONG)
  MESSAGE (WARNING "Tpetra: You set Tpetra_INST_INT_UNSIGNED_LONG.  This means that you want to use GlobalOrdinal = unsigned long with Tpetra objects.  We do NOT recommend using this type.  We are NOT responsible for ensuring that it builds or passes tests.  If you want a 64-bit type, we prefer that you use long long (Tpetra_INST_INT_LONG_LONG=ON), or long (Tpetra_INST_INT_LONG=ON) if you know that it is 64 bits.  If you are interested in using unsigned long, Roger Pawlowski is the responsible party; please contact him.")
ENDIF ()

# Set availability macros for various GlobalOrdinal types
GLOBAL_SET(HAVE_TPETRA_INT_INT ${Tpetra_INST_INT_INT})
GLOBAL_SET(HAVE_TPETRA_INT_LONG ${Tpetra_INST_INT_LONG})
GLOBAL_SET(HAVE_TPETRA_INT_LONG_LONG ${Tpetra_INST_INT_LONG_LONG})
GLOBAL_SET(HAVE_TPETRA_INT_UNSIGNED ${Tpetra_INST_INT_UNSIGNED})
GLOBAL_SET(HAVE_TPETRA_INT_UNSIGNED_LONG ${Tpetra_INST_INT_UNSIGNED_LONG})

# Check if no GlobalOrdinal types are enabled.  Print a "warning"
# (helpful message) in that case.  (Don't actually call it a WARNING;
# that might break builds.)
IF (Tpetra_INST_INT_INT OR Tpetra_INST_INT_UNSIGNED OR Tpetra_INST_INT_LONG OR Tpetra_INST_INT_LONG_LONG OR Tpetra_INST_INT_UNSIGNED_LONG)
  GLOBAL_SET (HAVE_TPETRA_GLOBALORDINAL ON)
ELSE()
  # mfh 18 Aug 2015: CASL wants to set GlobalOrdinal = unsigned long.
  # Tpetra doesn't have an option for this (for good reason: I don't
  # want to promise support for unsigned GlobalOrdinal types --
  # "unsigned" is already a bad idea).  Thus, CASL sets it using
  # TriBITS' ETI support.  As a result, their build disables ALL
  # GlobalOrdinal types that Tpetra claims to support; they explicitly
  # set the four Tpetra_INST_INT_* options above to OFF.  Tpetra's
  # CMakeLists.txt files get processed before its ETI logic, so the
  # code here doesn't see the GlobalOrdinal type that CASL adds.  This
  # is why we don't error out if none of the above options are ON.
  GLOBAL_SET (HAVE_TPETRA_GLOBALORDINAL OFF)
  MESSAGE (STATUS "Tpetra: NOTE: No GlobalOrdinal type is enabled.  This means that you, the user, explicitly disabled all GlobalOrdinal types that Tpetra enables by default.  If you know what you are doing, this might be OK.  If you don't know what you are doing, please don't do this.")
ENDIF()

# Append GlobalOrdinal types to instantiation /test list.  Add them to
# the list of Scalar types as well, since Tpetra needs Scalar =
# GlobalOrdinal for some communication routines.
IF(Tpetra_INST_INT_INT) 
  LIST(APPEND TpetraCore_ETI_SCALARS "int")
  LIST(APPEND TpetraCore_ETI_GORDS "int")
ENDIF()
IF(Tpetra_INST_INT_LONG) 
  LIST(APPEND TpetraCore_ETI_SCALARS "long")
  LIST(APPEND TpetraCore_ETI_GORDS "long")
ENDIF()
IF(Tpetra_INST_INT_LONG_LONG) 
  LIST(APPEND TpetraCore_ETI_SCALARS "long long")
  LIST(APPEND TpetraCore_ETI_GORDS "long long")
ENDIF()
IF(Tpetra_INST_INT_UNSIGNED) 
  LIST(APPEND TpetraCore_ETI_SCALARS "unsigned")
  LIST(APPEND TpetraCore_ETI_GORDS "unsigned")
ENDIF()
IF(Tpetra_INST_INT_UNSIGNED_LONG) 
  LIST(APPEND TpetraCore_ETI_SCALARS "unsigned long")
  LIST(APPEND TpetraCore_ETI_GORDS "unsigned long")
ENDIF()

# ============================================================
# Scalar types
# ============================================================

# The warning makes sense whether or not ETI is enabled.
IF (${PACKAGE_NAME}_ENABLE_quadmath AND ${PACKAGE_NAME}_ENABLE_CUDA)
  MESSAGE(WARNING "You have enabled the 'quadmath' TPL for __float128 support, but have also enabled CUDA.  __float128 does not work in CUDA code.  To simplify Tpetra's CMake logic, we have thus turned off Tpetra's ETI (explicit template instantiation) and tests for Scalar = __float128.  If you really want to use __float128 in Tpetra, AND really want to enable CUDA, please talk to the Tpetra developers.")
ENDIF ()

# Scalar = double is enabled by default.
GLOBAL_SET(Tpetra_INST_DOUBLE_DEFAULT ON)

# Set up the Tpetra_INST_DOUBLE option: whether to instantiate / test
# Scalar = double.
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_DOUBLE
  HAVE_TPETRA_INST_DOUBLE
  "Instantiate and/or test Tpetra classes over Scalar = double.  This option is ${Tpetra_INST_DOUBLE_DEFAULT} by default."
  ${Tpetra_INST_DOUBLE_DEFAULT}
  )
ASSERT_DEFINED (Tpetra_INST_DOUBLE)
GLOBAL_SET(HAVE_TPETRA_DOUBLE ${Tpetra_INST_DOUBLE})

# Scalar = float is NOT enabled by default.
GLOBAL_SET(Tpetra_INST_FLOAT_DEFAULT OFF)

# Set up the Tpetra_INST_FLOAT option: whether to instantiate / test
# Scalar = float.
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_FLOAT
  HAVE_TPETRA_INST_FLOAT
  "Instantiate and/or test classes over Scalar = float.  This option is ${Tpetra_INST_FLOAT_DEFAULT} by default."
  ${Tpetra_INST_FLOAT_DEFAULT}
  )
ASSERT_DEFINED (Tpetra_INST_FLOAT)
GLOBAL_SET(HAVE_TPETRA_FLOAT ${Tpetra_INST_FLOAT})

# Make sure that it's OK to instantiate/test over Scalar = float.
# Tpetra requires BLAS support for this case.  It presumes that a BLAS
# could exist without float (S) support.
ASSERT_DEFINED (HAVE_TEUCHOS_BLASFLOAT)
IF (Tpetra_INST_FLOAT AND NOT HAVE_TEUCHOS_BLASFLOAT)
  MESSAGE(FATAL_ERROR "Tpetra: Tpetra_INST_FLOAT is ON (meaning that you want to instantiate and/or test Tpetra with Scalar = float), but HAVE_TEUCHOS_BLASFLOAT is OFF.  This means that you are linking with a BLAS library that lacks float (S) support.  Tpetra needs a BLAS implementation that supports float.")
ENDIF ()

# Decide whether to instantiate / test Scalar = std::complex<float> or
# std::complex<double> by default.  Both require a complex-arithmetic
# BLAS.  Furthermore, std::complex<T> requires the corresponding
# Scalar=T instantiation to be ON.
#
# In a debug build, we enable Scalar = complex, in order that this
# capability be tested.  In a release build, we disable this case, in
# order to reduce build time.  Users may always enable it explicitly.
ASSERT_DEFINED (Tpetra_ENABLE_DEBUG)
# mfh 18 Aug 2015: Teuchos currently doesn't even define
# HAVE_COMPLEX_BLAS, in the case where one would expect it to be OFF.
# Thanks to Andrew Bradley for pointing this out.
IF (DEFINED HAVE_COMPLEX_BLAS)
  IF (HAVE_COMPLEX_BLAS)
    IF (Tpetra_INST_DOUBLE)
      GLOBAL_SET(Tpetra_INST_COMPLEX_DOUBLE_DEFAULT ${Tpetra_ENABLE_DEBUG})
    ELSE()
      GLOBAL_SET(Tpetra_INST_COMPLEX_DOUBLE_DEFAULT OFF)
    ENDIF()
    IF (Tpetra_INST_FLOAT)
      GLOBAL_SET(Tpetra_INST_COMPLEX_FLOAT_DEFAULT ${Tpetra_ENABLE_DEBUG})
    ELSE()
      GLOBAL_SET(Tpetra_INST_COMPLEX_FLOAT_DEFAULT OFF)
    ENDIF()
  ELSE()
    GLOBAL_SET(Tpetra_INST_COMPLEX_DOUBLE_DEFAULT OFF)
    GLOBAL_SET(Tpetra_INST_COMPLEX_FLOAT_DEFAULT OFF)
  ENDIF() # HAVE_COMPLEX_BLAS
ELSE () # NOT DEFINED HAVE_COMPLEX_BLAS (consider OFF)
  GLOBAL_SET(Tpetra_INST_COMPLEX_DOUBLE_DEFAULT OFF)
  GLOBAL_SET(Tpetra_INST_COMPLEX_FLOAT_DEFAULT OFF)
ENDIF () # DEFINED HAVE_COMPLEX_BLAS

# Set up the Tpetra_INST_COMPLEX_DOUBLE option: whether to instantiate
# / test Scalar = std::complex<double>.
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_COMPLEX_DOUBLE
  HAVE_TPETRA_INST_COMPLEX_DOUBLE
  "Instantiate and/or test Tpetra classes with Scalar = std::complex<double>.  This option is ${HAVE_TPETRA_INST_COMPLEX_DOUBLE_DEFAULT} by default."
  ${Tpetra_INST_COMPLEX_DOUBLE_DEFAULT}
  )
ASSERT_DEFINED (Tpetra_INST_COMPLEX_DOUBLE)
GLOBAL_SET (HAVE_TPETRA_COMPLEX_DOUBLE ${Tpetra_INST_COMPLEX_DOUBLE})

# Enabling Scalar = std::complex<T> requires that the corresponding
# Scalar=T instantiation be enabled as well.  
#
# FIXME (mfh 17 Aug 2015) It's not clear to me why this must be the
# case.
IF (Tpetra_INST_COMPLEX_DOUBLE AND NOT Tpetra_INST_DOUBLE)
  MESSAGE (FATAL_ERROR "Tpetra: Tpetra_INST_COMPLEX_DOUBLE is ON, meaning that you want to want to instantiate and/or test Tpetra classes with Scalar = std::complex<double>.  However, Tpetra_INST_DOUBLE is OFF.  Enabling Scalar = std::complex<T> requires that the corresponding Scalar=T real instantiation be enabled as well.")
ENDIF ()

# Tpetra requires that the BLAS work for Scalar, as long as Scalar is
# one of the four types (S, D, C, Z) that the BLAS promises to
# implement.  Check whether the BLAS works for complex arithmetic.

# mfh 18 Aug 2015: Teuchos currently doesn't even define
# HAVE_COMPLEX_BLAS, in the case where one would expect it to be OFF.
# Thanks to Andrew Bradley for pointing this out.
IF (Tpetra_INST_COMPLEX_DOUBLE)
  IF (NOT DEFINED HAVE_COMPLEX_BLAS)
    MESSAGE(FATAL_ERROR "Tpetra: Tpetra_INST_COMPLEX_DOUBLE is ON (meaning that you want to want to instantiate and/or test Tpetra classes with Scalar = std::complex<double>), but HAVE_COMPLEX_BLAS is not defined.  This means that you are linking with a BLAS library that lacks complex arithmetic (Z) support.  Tpetra needs a BLAS implementation that supports complex arithmetic.")
  ELSEIF (NOT HAVE_COMPLEX_BLAS)
    MESSAGE(FATAL_ERROR "Tpetra: Tpetra_INST_COMPLEX_DOUBLE is ON (meaning that you want to want to instantiate and/or test Tpetra classes with Scalar = std::complex<double>), but HAVE_COMPLEX_BLAS is OFF.  This means that you are linking with a BLAS library that lacks complex arithmetic (Z) support.  Tpetra needs a BLAS implementation that supports complex arithmetic.")
  ENDIF ()
ENDIF ()

# Set up the Tpetra_INST_COMPLEX_FLOAT option: whether to instantiate
# / test Scalar = std::complex<float>.
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_COMPLEX_FLOAT
  HAVE_TPETRA_INST_COMPLEX_FLOAT
  "Instantiate and/or test Tpetra classes with Scalar = std::complex<float>.  This option is ${HAVE_TPETRA_INST_COMPLEX_FLOAT_DEFAULT} by default."
  ${Tpetra_INST_COMPLEX_FLOAT_DEFAULT}
  )
ASSERT_DEFINED (Tpetra_INST_COMPLEX_FLOAT)
GLOBAL_SET (HAVE_TPETRA_COMPLEX_FLOAT ${Tpetra_INST_COMPLEX_FLOAT})

# Enabling Scalar = std::complex<T> requires that the corresponding
# Scalar=T instantiation be enabled as well.
#
# FIXME (mfh 17 Aug 2015) It's not clear to me why this must be the
# case.
IF (Tpetra_INST_COMPLEX_FLOAT AND NOT Tpetra_INST_FLOAT)
  MESSAGE (FATAL_ERROR "Tpetra: Tpetra_INST_COMPLEX_FLOAT is ON, meaning that you want to want to instantiate and/or test Tpetra classes with Scalar = std::complex<float>.  However, Tpetra_INST_FLOAT is OFF.  Enabling Scalar = std::complex<T> requires that the corresponding Scalar=T real instantiation be enabled as well.")
ENDIF ()

# Tpetra requires that the BLAS work for Scalar, as long as Scalar is
# one of the four types (S, D, C, Z) that the BLAS promises to
# implement.  Check whether the BLAS works for complex float (C)
# arithmetic.
IF (Tpetra_INST_COMPLEX_FLOAT)
  IF (NOT DEFINED HAVE_COMPLEX_BLAS)
    MESSAGE(FATAL_ERROR "Tpetra: Tpetra_INST_COMPLEX_FLOAT is ON (meaning that you want to want to instantiate and/or test Tpetra classes with Scalar = std::complex<float>), but HAVE_COMPLEX_BLAS is not defined.  This means that you are linking with a BLAS library that lacks complex arithmetic (Z) support.  Tpetra needs a BLAS implementation that supports complex arithmetic.")
  ELSEIF (NOT HAVE_COMPLEX_BLAS)
    MESSAGE(FATAL_ERROR "Tpetra: Tpetra_INST_COMPLEX_FLOAT is ON (meaning that you want to want to instantiate and/or test Tpetra classes with Scalar = std::complex<float>), but HAVE_COMPLEX_BLAS is OFF.  This means that you are linking with a BLAS library that lacks complex arithmetic (Z) support.  Tpetra needs a BLAS implementation that supports complex arithmetic.")
  ENDIF ()
ENDIF ()

# Optionally enable explicit template instantiation (ETI) and tests
# for Scalar = __float128.  __float128 is a GCC-specific C++ language
# extension.  It requires building with GCC, enabling the 'quadmath'
# TPL, and setting build flags -std=gnu+11 (not c++11, else warnings)
# and -fext-numeric-literals (else build errors).  Be sure also to
# turn off -ansi -pedantic; __float128 is a language extension, so it
# is NOT ANSI C++.
#
# __float128 won't work with Kokkos::Cuda.  We could try to be clever
# about turning off Tpetra's ETI and tests for __float128 with the
# Kokkos::Cuda execution space, but instead, we simply disable
# __float128 ETI and tests if CUDA is enabled.
ASSERT_DEFINED (${PACKAGE_NAME}_ENABLE_quadmath)
ASSERT_DEFINED (${PACKAGE_NAME}_ENABLE_CUDA)
IF (${PACKAGE_NAME}_ENABLE_quadmath AND NOT ${PACKAGE_NAME}_ENABLE_CUDA)
  # FIXME (mfh 19 Mar 2015) It doesn't currently build.
  # Once it does, set this to ON.
  SET (Tpetra_INST_FLOAT128_DEFAULT OFF)
ELSE ()
  SET (Tpetra_INST_FLOAT128_DEFAULT OFF)
ENDIF ()

# Set up the Tpetra_INST_FLOAT128 option: whether to instantiate /
# test Scalar = std::complex<float>.  It is OFF by default.  Users may
# shut this off explicitly by setting Tpetra_INST_FLOAT128=OFF.
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_FLOAT128
  HAVE_TPETRA_INST_FLOAT128
  "Instantiate and/or test Tpetra classes with Scalar = __float128.  This option is ${HAVE_TPETRA_INST_FLOAT128_DEFAULT} by default.  __float128 is a GCC-specific C++ language extension.  It requires building with GCC, enabling the 'quadmath' TPL, and setting build flags -std=gnu+11 (not c++11, else warnings) and -fext-numeric-literals (else build errors).  Be sure also to turn off -ansi -pedantic; __float128 is a language extension, so it is NOT ANSI C++."
  ${Tpetra_INST_FLOAT128_DEFAULT}
  )
ASSERT_DEFINED (Tpetra_INST_FLOAT128)
GLOBAL_SET (HAVE_TPETRA_FLOAT128 ${Tpetra_INST_FLOAT128})

IF (Tpetra_INST_FLOAT128 AND NOT ${PACKAGE_NAME}_ENABLE_quadmath)
  MESSAGE (FATAL_ERROR "You may not enable ETI (explicit template instantiation) or tests in Tpetra for Scalar = __float128, unless you enable the 'quadmath' TPL.  Please also note that __float128 is a GCC-specific C++ language extension.  It requires building with GCC, enabling the 'quadmath' TPL, and setting build flags -std=gnu+11 (not c++11, else warnings) and -fext-numeric-literals (else build errors).  Be sure also to turn off -ansi -pedantic; __float128 is a language extension, so it is NOT ANSI C++.")
ENDIF ()

# Decide whether to instantiate / test Scalar = qd_real ("quad
# double", provided by the QD TPL) by default.  It requires the QD TPL
# and does not work with the Kokkos refactor version of Tpetra (due to
# the missing Kokkos::atomic_* and volatile overloaded functions /
# methods).
ASSERT_DEFINED (${PACKAGE_NAME}_ENABLE_QD)
ASSERT_DEFINED (Tpetra_ENABLE_Kokkos_Refactor)
IF (${PACKAGE_NAME}_ENABLE_QD AND NOT Tpetra_ENABLE_Kokkos_Refactor)
  GLOBAL_SET(HAVE_TPETRA_INST_QD_REAL_DEFAULT ON)
ELSE()
  GLOBAL_SET(HAVE_TPETRA_INST_QD_REAL_DEFAULT OFF)
ENDIF()

# Set up the Tpetra_INST_QD_REAL option: whether to instantiate / test
# Scalar = qd_real.
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_QD_REAL
  HAVE_TPETRA_INST_QD_REAL
  "Instantiate and/or test Tpetra classes over Scalar = qd_real (requries TPL QD).  This option is ${HAVE_TPETRA_INST_QD_REAL_DEFAULT} by default."
  ${HAVE_TPETRA_INST_QD_REAL_DEFAULT}  
  )
ASSERT_DEFINED (Tpetra_INST_QD_REAL)
GLOBAL_SET (HAVE_TPETRA_QD_REAL ${Tpetra_INST_QD_REAL})

IF (Tpetra_INST_QD_REAL AND Tpetra_ENABLE_Kokkos_Refactor)
  MESSAGE (SEND_ERROR "Scalar = qd_real does not currently work with the Kokkos refactor version of Tpetra.  We will fix this at some point.  Please let us know if you need this capability urgently.")
ENDIF ()

# Decide whether to instantiate / test Scalar = dd_real ("double
# double", provided by the QD TPL) by default.  It requires the QD TPL
# and does not work with the Kokkos refactor version of Tpetra (due to
# the missing Kokkos::atomic_* and volatile overloaded functions /
# methods).
IF (${PACKAGE_NAME}_ENABLE_QD AND NOT Tpetra_ENABLE_Kokkos_Refactor)
  GLOBAL_SET(HAVE_TPETRA_INST_DD_REAL_DEFAULT ON)
ELSE()
  GLOBAL_SET(HAVE_TPETRA_INST_DD_REAL_DEFAULT OFF)
ENDIF()

# Set up the Tpetra_INST_DD_REAL option: whether to instantiate / test
# Scalar = dd_real.
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_INST_DD_REAL
  HAVE_TPETRA_INST_DD_REAL
  "Instantiate and/or test Tpetra classes over Scalar = dd_real (requries TPL QD).  This option is ${HAVE_TPETRA_INST_DD_REAL_DEFAULT} by default."
  ${HAVE_TPETRA_INST_DD_REAL_DEFAULT}  
  )
ASSERT_DEFINED (Tpetra_INST_DD_REAL)
GLOBAL_SET (HAVE_TPETRA_DD_REAL ${Tpetra_INST_DD_REAL})

IF (Tpetra_INST_DD_REAL AND Tpetra_ENABLE_Kokkos_Refactor)
  MESSAGE (SEND_ERROR "Scalar = dd_real does not currently work with the Kokkos refactor version of Tpetra.  We will fix this at some point.  Please let us know if you need this capability urgently.")
ENDIF ()

#
# Build list of Scalar types over which to instantiate / test.
#
IF(Tpetra_INST_FLOAT) 
  LIST(APPEND TpetraCore_ETI_SCALARS "float")
  LIST(APPEND TpetraCore_ETI_SCALARS_NO_ORDS "float")
  GLOBAL_SET (HAVE_TPETRA_SCALAR ON)
ENDIF()
IF(Tpetra_INST_DOUBLE) 
  LIST(APPEND TpetraCore_ETI_SCALARS "double")
  LIST(APPEND TpetraCore_ETI_SCALARS_NO_ORDS "double")
  GLOBAL_SET (HAVE_TPETRA_SCALAR ON)
ENDIF()
IF(Tpetra_INST_COMPLEX_FLOAT) 
  LIST(APPEND TpetraCore_ETI_SCALARS "std::complex<float>")
  LIST(APPEND TpetraCore_ETI_SCALARS_NO_ORDS "std::complex<float>")
  GLOBAL_SET (HAVE_TPETRA_SCALAR ON)
ENDIF()
IF(Tpetra_INST_COMPLEX_DOUBLE) 
  LIST(APPEND TpetraCore_ETI_SCALARS "std::complex<double>")
  LIST(APPEND TpetraCore_ETI_SCALARS_NO_ORDS "std::complex<double>")
  GLOBAL_SET (HAVE_TPETRA_SCALAR ON)
ENDIF()
IF(Tpetra_INST_FLOAT128) 
  LIST(APPEND TpetraCore_ETI_SCALARS "__float128")
  LIST(APPEND TpetraCore_ETI_SCALARS_NO_ORDS "__float128")
  GLOBAL_SET (HAVE_TPETRA_SCALAR ON)
ENDIF()
IF(Tpetra_INST_QD_REAL) 
  LIST(APPEND TpetraCore_ETI_SCALARS "qd_real")
  LIST(APPEND TpetraCore_ETI_SCALARS_NO_ORDS "qd_real")
  GLOBAL_SET (HAVE_TPETRA_SCALAR ON)
ENDIF()
IF(Tpetra_INST_DD_REAL) 
  LIST(APPEND TpetraCore_ETI_SCALARS "dd_real")
  LIST(APPEND TpetraCore_ETI_SCALARS_NO_ORDS "dd_real")
  GLOBAL_SET (HAVE_TPETRA_SCALAR ON)
ENDIF()

# Make sure that at one Scalar type got activated
IF (NOT HAVE_TPETRA_SCALAR)
  MESSAGE (FATAL_ERROR "You must enable at least one Scalar type in Tpetra.  If you see this message, you probably set nondefault CMake options to disable all of Tpetra's Scalar types.  Please don't do that.")
ENDIF ()

# Make sure that Tpetra enables at least one Scalar type.
ASSERT_DEFINED(HAVE_TPETRA_SCALAR)

#
# Other Tpetra options
#

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
  "This option has been DEPRECATED and its value will be IGNORED.  We keep it only for backwards compatibility."
  OFF
  )
IF (Tpetra_ENABLE_Kokkos_Refactor_Map)
  MESSAGE(FATAL_ERROR "The Tpetra_ENABLE_Kokkos_Refactor_Map option has been DEPRECATED and must no longer be set.  It turns out that we took a different path with the refactor of Tpetra::Map, that didn't require two separate versions of the class.  Thus, there is no need to ask for the Kokkos refactor version of Map; there was only ever one version of Map.")
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
  MESSAGE(FATAL_ERROR "You enabled the Reduction/Transformation Interface (RTI) in Tpetra by setting ${PACKAGE_NAME}_ENABLE_RTI, yet the Kokkos refactor (new) version of Tpetra is enabled.  We have not yet certified RTI for use with the new version of Tpetra.  If you really need RTI right now, you may bring back the old version of Tpetra temporarily by setting Tpetra_ENABLE_Kokkos_Refactor:BOOL=OFF.  Please note that the old version of Tpetra will go away completely with the 12.0 release of Trilinos.")
ENDIF ()

#
# Add libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(inout)
ADD_SUBDIRECTORY(ext)
IF (${PACKAGE_NAME}_ENABLE_RTI)
  ADD_SUBDIRECTORY(rti)
ENDIF()
 
TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

GLOBAL_SET(TpetraCore_ETI_SCALARS ${TpetraCore_ETI_SCALARS})
GLOBAL_SET(TpetraCore_ETI_SCALARS_NO_ORDS ${TpetraCore_ETI_SCALARS_NO_ORDS})
GLOBAL_SET(TpetraCore_ETI_GORDS   ${TpetraCore_ETI_GORDS})
GLOBAL_SET(TpetraCore_ETI_LORDS   ${TpetraCore_ETI_LORDS})
GLOBAL_SET(TpetraCore_ETI_NODES   ${TpetraCore_ETI_NODES})

#
# Do standard subpackage postprocessing
#
TRIBITS_SUBPACKAGE_POSTPROCESS()
