############################################################
# Define the package
############################################################

# NOTE: If the TBB (Intel's Threading Building Blocks) TPL is enabled,
# it's probably not possible to build Tpetra with shadowing warnings.
#TRIBITS_PACKAGE( Tpetra ENABLE_SHADOWING_WARNINGS CLEANED )
TRIBITS_PACKAGE( Tpetra CLEANED )

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

TRIBITS_SET_ST_FOR_DEV_MODE(ST_FOR_DEV_PT_FOR_RELEASE)

# Kokkos and Tpetra require C++11 as of the 12.0 release.
#
# I left a special case for CUDA because I'm not sure whether TriBITS
# counts C++11 as enabled in that case; please feel free to revise if
# you know better.
ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX11)
ASSERT_DEFINED(TPL_ENABLE_CUDA)
IF (NOT ${PROJECT_NAME}_ENABLE_CXX11 AND NOT TPL_ENABLE_CUDA)
  MESSAGE(FATAL_ERROR "As of Trilinos 12.0, C++11 is REQUIRED when building Kokkos and its downstream packages, which include Tpetra.  You MUST enable C++11 when building these packages.  Trilinos should enable this by default, if your C++11 compiler supports C++11 (that is, if it's not too old).  If you set Trilinos_ENABLE_CXX11=OFF, then you are disabling C++11.  DON'T DO THAT if you want to build Kokkos or its downstream packages, including Tpetra.  If you aren't happy about that, then you need to disable Kokkos explicitly, by setting Trilinos_ENABLE_Kokkos=OFF.  This WILL disable Tpetra and its downstream packages.")
ENDIF ()

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_ENABLE_Kokkos_Refactor
  TPETRA_ENABLE_KOKKOS_REFACTOR
  "Enable the \"Kokkos refactor\" version of Tpetra.  DO NOT SET THIS EXPLICITLY!  It is ON by default and that is the only supported value."
  ON
  )
IF (NOT Tpetra_ENABLE_Kokkos_Refactor)
  MESSAGE(FATAL_ERROR "The 'classic' version of Tpetra has been deprecated for a long time, and is not longer available.  If you don't want to see this message, please stop setting Tpetra_ENABLE_Kokkos_Refactor to OFF.  (It is ON by default; you should not set it explicitly.)")
ENDIF ()

######################################################################
# Explicit template instantiation (ETI) and test instantiation logic
######################################################################

# For
#   - Scalar types,
#   - (LocalOrdinal, GlobalOrdinal) type pairs, and
#   - Node types,
#
# Tpetra defines macros which say
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
# The node types are:    SERIAL, OPENMP, PTHREAD, CUDA
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
#
# Tpetra also keeps a list of the Kokkos Device types (specializations
# of Kokkos::Device) over which it does instantiations and/or tests.
# This is the Tpetra_ETI_DEVICES CMake variable.  Eventually, we will
# deprecate and remove the above Node instantiation list, and just use
# the list of enabled Devices.  We use devices rather than execution
# spaces, because a Kokkos::Device is an (execution space, memory
# space) pair.  This means that it is more general.  For example,
# Device<Cuda, CudaSpace> != Device<Cuda, CudaUVMSpace>, yet some
# users may want both represented.  The particular use case I have in
# mind here is Device<OpenMP, $MEM_SPACE> where $MEM_SPACE could be
# either HostSpace or the HBM (high-bandwidth memory) memory space.
#

SET(${PACKAGE_NAME}_ETI_SCALARS "")
# Exclude all ordinal types (GlobalOrdinal and int)
SET(${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS "")
SET(${PACKAGE_NAME}_ETI_LORDS "")
SET(${PACKAGE_NAME}_ETI_GORDS "")
SET(${PACKAGE_NAME}_ETI_DEVICES "")
SET(${PACKAGE_NAME}_ETI_NODES "")

# ============================================================
# Kokkos execution space types
# ============================================================

# Decide what Kokkos execution spaces to use by default.  If ETI is
# ON, "using an execution space ExecutionSpace" means that Tpetra
# objects with Node =
# Kokkos::Compat::KokkosDeviceWrapperNode<ExecutionSpace> will
#
#   1. get instantiated explicitly for that Node type, and
#   2. get tested for that Node type, if their test is also
#      templated on Node type.
#
# If ETI is OFF, #1 no longer holds, but #2 still holds.
#
# Tpetra uses exactly one Kokkos execution space by default, whether
# ETI is ON or OFF.  This keeps build times for tests down in the
# non-ETI case, and library sizes small in the ETI case.

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
  IF (HAVE_TPETRA_DEFAULTNODE_SERIALWRAPPERNODE)
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
  IF (HAVE_TPETRA_DEFAULTNODE_OPENMPWRAPPERNODE)
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
  IF (HAVE_TPETRA_DEFAULTNODE_THREADSWRAPPERNODE) 
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
  IF (HAVE_TPETRA_DEFAULTNODE_CUDAWRAPPERNODE) 
    MESSAGE(FATAL_ERROR "Tpetra: Node = Kokkos::Compat::KokkosCudaWrapperNode is disabled (since Tpetra_INST_CUDA=OFF), but you set it as the default Node type.  Try setting the CMake options Kokkos_ENABLE_Cuda:BOOL=ON and Tpetra_INST_CUDA:BOOL=ON.  If you are building with a CUDA-capable compiler and Kokkos can detect that, then you are unlikely to see this message, since both Kokkos and Tpetra enable CUDA support by default in that case.")
  ENDIF()
  IF (Kokkos_ENABLE_Cuda)
    MESSAGE(STATUS "NOTE: Kokkos::Cuda is ON (the CMake option Kokkos_ENABLE_Cuda is ON), but the corresponding Tpetra Node type is disabled.  If you want to enable instantiation and use of Kokkos::Cuda in Tpetra, please also set the CMake option Tpetra_INST_CUDA:BOOL=ON.  If you use the Kokkos::Cuda version of Tpetra without doing this, you will get link errors!")
  ENDIF ()
ENDIF () # Tpetra_INST_CUDA

#
# Construct the list of enabled Device and Node types.
#
IF(Tpetra_INST_SERIAL) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_DEVICES "Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>")
  LIST(APPEND ${PACKAGE_NAME}_ETI_NODES "Kokkos::Compat::KokkosSerialWrapperNode")
ENDIF()
IF(Tpetra_INST_PTHREAD) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_DEVICES "Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>")
  LIST(APPEND ${PACKAGE_NAME}_ETI_NODES "Kokkos::Compat::KokkosThreadsWrapperNode")
ENDIF()
IF(Tpetra_INST_OPENMP) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_DEVICES "Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>")
  LIST(APPEND ${PACKAGE_NAME}_ETI_NODES "Kokkos::Compat::KokkosOpenMPWrapperNode")
ENDIF()
IF(Tpetra_INST_CUDA) 
  # TpetraCore currently requires UVM, so it doesn't make sense to
  # instantiate over Kokkos::CudaSpace as well as
  # Kokkos::CudaUVMSpace.
  LIST(APPEND ${PACKAGE_NAME}_ETI_DEVICES "Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>")
  LIST(APPEND ${PACKAGE_NAME}_ETI_NODES "Kokkos::Compat::KokkosCudaWrapperNode")
ENDIF()

# Tell users what Nodes are enabled.
MESSAGE(STATUS "Tpetra execution space availability (ON means available): ")
MESSAGE(STATUS "  - Serial:  ${HAVE_TPETRA_SERIAL}")
MESSAGE(STATUS "  - Threads: ${HAVE_TPETRA_PTHREAD}")
MESSAGE(STATUS "  - OpenMP:  ${HAVE_TPETRA_OPENMP}")
MESSAGE(STATUS "  - Cuda:    ${HAVE_TPETRA_CUDA}")

# Fix Github Issue #190 by making sure that users enabled at least one
# Node type.
IF(NOT Tpetra_INST_SERIAL AND NOT Tpetra_INST_PTHREAD AND NOT Tpetra_INST_OPENMP AND NOT Tpetra_INST_CUDA)
  MESSAGE(FATAL_ERROR "It appears that you have disabled all of Tpetra's Node types.  This is a bad idea; please don't do this.  This may have happened either if all Kokkos execution spaces got disabled, or if you explicitly disabled some Tpetra Node types.  To fix this, look in your CMake configuration to see if you set any of the following variables to OFF explicitly:  Kokkos_ENABLE_Cuda, Kokkos_ENABLE_OpenMP, Kokkos_ENABLE_Serial, Kokkos_ENABLE_Pthread, Tpetra_INST_SERIAL, Tpetra_INST_PTHREAD, Tpetra_INST_OPENMP, or Tpetra_INST_CUDA.")
ENDIF()

# ============================================================
# LocalOrdinal types
# ============================================================

# FIXME (mfh 16 Oct 2015) Tpetra currently only supports LocalOrdinal
# = int.  At some point, the list of allowed LocalOrdinal types will
# expand to include 64-bit signed integer types.

LIST (APPEND ${PACKAGE_NAME}_ETI_LORDS "int")

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
# Tpetra has support for GlobalOrdinal = unsigned and unsigned long,
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

#
# Construct the list of enabled GlobalOrdinal types.  Add its entries
# to the list of Scalar types as well, since Tpetra needs Scalar =
# GlobalOrdinal for some communication routines.
#
IF(Tpetra_INST_INT_INT) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS "int")
  LIST(APPEND ${PACKAGE_NAME}_ETI_GORDS "int")
ENDIF()
IF(Tpetra_INST_INT_LONG) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS "long")
  LIST(APPEND ${PACKAGE_NAME}_ETI_GORDS "long")
ENDIF()
IF(Tpetra_INST_INT_LONG_LONG) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS "long long")
  LIST(APPEND ${PACKAGE_NAME}_ETI_GORDS "long long")
ENDIF()
IF(Tpetra_INST_INT_UNSIGNED) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS "unsigned")
  LIST(APPEND ${PACKAGE_NAME}_ETI_GORDS "unsigned")
ENDIF()
IF(Tpetra_INST_INT_UNSIGNED_LONG) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS "unsigned long")
  LIST(APPEND ${PACKAGE_NAME}_ETI_GORDS "unsigned long")
ENDIF()

# ============================================================
# Scalar types
# ============================================================

# Tpetra decides what Scalar types to enable, based in part on which
# TPLs are enabled.
#
# mfh 16 Oct 2015: In TriBITS, a package depends on its subpackages.
# Thus, it's OK for Tpetra the package to check a variable defined in
# its subpackage TpetraCore.  Nevertheless, we cautiously check
# whether the variable in question is defined before using it.

IF (NOT DEFINED ${PACKAGE_NAME}_ENABLE_quadmath)
  IF (DEFINED TpetraCore_ENABLE_quadmath)
    SET (${PACKAGE_NAME}_ENABLE_quadmath "${TpetraCore_ENABLE_quadmath}")
  ELSEIF (DEFINED TPL_ENABLE_quadmath)
    SET (${PACKAGE_NAME}_ENABLE_quadmath "${TPL_ENABLE_quadmath}")
  ELSE ()
    SET (${PACKAGE_NAME}_ENABLE_quadmath OFF)
  ENDIF ()
ENDIF ()
ASSERT_DEFINED (${PACKAGE_NAME}_ENABLE_quadmath)

IF (NOT DEFINED ${PACKAGE_NAME}_ENABLE_CUDA)
  IF (DEFINED TpetraCore_ENABLE_CUDA)
    SET (${PACKAGE_NAME}_ENABLE_CUDA "${TpetraCore_ENABLE_CUDA}")
  ELSEIF (DEFINED TPL_ENABLE_CUDA)
    SET (${PACKAGE_NAME}_ENABLE_CUDA "${TPL_ENABLE_CUDA}")
  ELSE ()
    SET (${PACKAGE_NAME}_ENABLE_CUDA OFF)
  ENDIF ()
ENDIF ()
ASSERT_DEFINED (${PACKAGE_NAME}_ENABLE_CUDA)

IF (NOT DEFINED ${PACKAGE_NAME}_ENABLE_QD)
  IF (DEFINED TpetraCore_ENABLE_QD)
    SET (${PACKAGE_NAME}_ENABLE_QD "${TpetraCore_ENABLE_QD}")
  ELSEIF (DEFINED TPL_ENABLE_QD)
    SET (${PACKAGE_NAME}_ENABLE_QD "${TPL_ENABLE_QD}")
  ELSE ()
    SET (${PACKAGE_NAME}_ENABLE_QD OFF)
  ENDIF ()
ENDIF ()
ASSERT_DEFINED (${PACKAGE_NAME}_ENABLE_QD)

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
# could exist without float (S) support.  Furthermore, it presumes
# that Teuchos supports float (this is mainly necessary for correct
# linking when ETI is ON, with downstream packages like Thyra, which
# get used in surprising places, e.g., in Ifpack2).
IF (Tpetra_INST_FLOAT)
  IF (NOT DEFINED Teuchos_ENABLE_FLOAT)
    MESSAGE(FATAL_ERROR "Tpetra: Tpetra_INST_FLOAT is ON (meaning that you want to instantiate and/or test Tpetra with Scalar = float), but Teuchos_ENABLE_FLOAT is not defined.  Tpetra needs Teuchos to support float.  Please set Teuchos_ENABLE_FLOAT:BOOL=ON, reconfigure, and rebuild.")
  ELSEIF (NOT Teuchos_ENABLE_FLOAT)
    MESSAGE(FATAL_ERROR "Tpetra: Tpetra_INST_FLOAT is ON (meaning that you want to instantiate and/or test Tpetra with Scalar = float), but Teuchos_ENABLE_FLOAT is OFF.  Tpetra needs Teuchos to support float.  Please set Teuchos_ENABLE_FLOAT:BOOL=ON, reconfigure, and rebuild.")
  ENDIF ()

  IF (NOT DEFINED HAVE_TEUCHOS_BLASFLOAT)
    MESSAGE(FATAL_ERROR "Tpetra: Tpetra_INST_FLOAT is ON (meaning that you want to instantiate and/or test Tpetra with Scalar = float), but HAVE_TEUCHOS_BLASFLOAT is not defined.  This means that you are linking with a BLAS library that lacks float (S) support.  Tpetra needs a BLAS implementation that supports float.")
  ELSEIF (NOT HAVE_TEUCHOS_BLASFLOAT)
    MESSAGE(FATAL_ERROR "Tpetra: Tpetra_INST_FLOAT is ON (meaning that you want to instantiate and/or test Tpetra with Scalar = float), but HAVE_TEUCHOS_BLASFLOAT is OFF.  This means that you are linking with a BLAS library that lacks float (S) support.  Tpetra needs a BLAS implementation that supports float.")
  ENDIF ()
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

IF (Tpetra_INST_COMPLEX_DOUBLE)
  IF (NOT DEFINED Teuchos_ENABLE_COMPLEX)
    MESSAGE(FATAL_ERROR "Tpetra: Tpetra_INST_COMPLEX_DOUBLE is ON (meaning that you want to want to instantiate and/or test Tpetra classes with Scalar = std::complex<double>), but Teuchos_ENABLE_COMPLEX is not defined.  Please set Teuchos_ENABLE_COMPLEX:BOOL=ON, reconfigure, and rebuild.")
  ELSEIF (NOT Teuchos_ENABLE_COMPLEX)
    MESSAGE(FATAL_ERROR "Tpetra: Tpetra_INST_COMPLEX_DOUBLE is ON (meaning that you want to want to instantiate and/or test Tpetra classes with Scalar = std::complex<double>), but Teuchos_ENABLE_COMPLEX is OFF.  Please set Teuchos_ENABLE_COMPLEX:BOOL=ON, reconfigure, and rebuild.")
  ENDIF ()

  # mfh 18 Aug 2015: Teuchos currently doesn't even define
  # HAVE_COMPLEX_BLAS, in the case where one would expect it to be
  # OFF.  Thanks to Andrew Bradley for pointing this out.  This is why
  # we start with a check for whether it is NOT DEFINED; otherwise,
  # testing "IF (NOT HAVE_COMPLEX_BLAS)" would be a CMake syntax
  # error.
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
# arithmetic.  Also, check whether Teuchos supports complex arithmetic
# at all.
IF (Tpetra_INST_COMPLEX_FLOAT)
  IF (NOT DEFINED Teuchos_ENABLE_COMPLEX)
    MESSAGE(FATAL_ERROR "Tpetra: Tpetra_INST_COMPLEX_FLOAT is ON (meaning that you want to want to instantiate and/or test Tpetra classes with Scalar = std::complex<float>), but Teuchos_ENABLE_COMPLEX is not defined.  Please set Teuchos_ENABLE_COMPLEX:BOOL=ON, reconfigure, and rebuild.")
  ELSEIF (NOT Teuchos_ENABLE_COMPLEX)
    MESSAGE(FATAL_ERROR "Tpetra: Tpetra_INST_COMPLEX_FLOAT is ON (meaning that you want to want to instantiate and/or test Tpetra classes with Scalar = std::complex<float>), but Teuchos_ENABLE_COMPLEX is OFF.  Please set Teuchos_ENABLE_COMPLEX:BOOL=ON, reconfigure, and rebuild.")
  ENDIF ()

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
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS "float")
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS "float")
  GLOBAL_SET (HAVE_TPETRA_SCALAR ON)
ENDIF()
IF(Tpetra_INST_DOUBLE) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS "double")
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS "double")
  GLOBAL_SET (HAVE_TPETRA_SCALAR ON)
ENDIF()
IF(Tpetra_INST_COMPLEX_FLOAT) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS "std::complex<float>")
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS "std::complex<float>")
  GLOBAL_SET (HAVE_TPETRA_SCALAR ON)
ENDIF()
IF(Tpetra_INST_COMPLEX_DOUBLE) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS "std::complex<double>")
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS "std::complex<double>")
  GLOBAL_SET (HAVE_TPETRA_SCALAR ON)
ENDIF()
IF(Tpetra_INST_FLOAT128) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS "__float128")
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS "__float128")
  GLOBAL_SET (HAVE_TPETRA_SCALAR ON)
ENDIF()
IF(Tpetra_INST_QD_REAL) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS "qd_real")
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS "qd_real")
  GLOBAL_SET (HAVE_TPETRA_SCALAR ON)
ENDIF()
IF(Tpetra_INST_DD_REAL) 
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS "dd_real")
  LIST(APPEND ${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS "dd_real")
  GLOBAL_SET (HAVE_TPETRA_SCALAR ON)
ENDIF()

# Make sure that at one Scalar type got activated
IF (NOT HAVE_TPETRA_SCALAR)
  MESSAGE (FATAL_ERROR "You must enable at least one Scalar type in Tpetra.  If you see this message, you probably set nondefault CMake options to disable all of Tpetra's Scalar types.  Please don't do that.")
ENDIF ()

# Make sure that Tpetra enables at least one Scalar type.
ASSERT_DEFINED(HAVE_TPETRA_SCALAR)

GLOBAL_SET(${PACKAGE_NAME}_ETI_SCALARS ${${PACKAGE_NAME}_ETI_SCALARS})
GLOBAL_SET(${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS ${${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS})
GLOBAL_SET(${PACKAGE_NAME}_ETI_LORDS   ${${PACKAGE_NAME}_ETI_LORDS})
GLOBAL_SET(${PACKAGE_NAME}_ETI_GORDS   ${${PACKAGE_NAME}_ETI_GORDS})
GLOBAL_SET(${PACKAGE_NAME}_ETI_DEVICES ${${PACKAGE_NAME}_ETI_DEVICES})
GLOBAL_SET(${PACKAGE_NAME}_ETI_NODES   ${${PACKAGE_NAME}_ETI_NODES})


############################################################
# Process subpackages
############################################################

TRIBITS_PROCESS_SUBPACKAGES()

TRIBITS_PACKAGE_DEF()

#
# Identify what files to exclude when making a source tarball.
#
TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_EXCLUDE_FILES(
  core/doc/CodingGuidelines
  core/doc/TpetraDesign
  )

#
# Do standard package postprocessing
#
TRIBITS_PACKAGE_POSTPROCESS()
