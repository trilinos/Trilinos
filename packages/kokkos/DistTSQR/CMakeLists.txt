# Make sure we added this option to the build
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TSQR)

IF(${PACKAGE_NAME}_ENABLE_TSQR)

  SET(HEADERS "")
  SET(SOURCES "")
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

  # Non-generated header files
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  # "Tsqr_Mpi*" files will only have nontrivial contents if HAVE_MPI
  # is #defined via Teuchos_ConfigDefs.hpp.
  APPEND_SET(HEADERS
    Tsqr_DistTsqrHelper.hpp
    Tsqr_DistTsqr.hpp
    Tsqr_DistTsqrRB.hpp
    Tsqr_FullTsqrTest.hpp
    Tsqr_generateStack.hpp
    Tsqr_GlobalTimeStats.hpp
    Tsqr_GlobalVerify.hpp
    Tsqr.hpp
    Tsqr_MessengerBase.hpp 
    Tsqr_Mgs.hpp 
    Tsqr_MgsTest.hpp 
    Tsqr_MpiDatatype.hpp 
    Tsqr_MpiMessenger.hpp 
    Tsqr_ParTest.hpp 
    Tsqr_printGlobalMatrix.hpp 
    Tsqr_Random_GlobalMatrix.hpp
    Tsqr_RMessenger.hpp
    Tsqr_TestSetup.hpp
    Tsqr_TestUtils.hpp
    Tsqr_TeuchosMessenger.hpp
    Tsqr_TrivialMessenger.hpp
    Tsqr_TsqrTest.hpp
    )
  # "Tsqr_Mpi*" files will only have nontrivial contents if HAVE_MPI
  # is #defined via Teuchos_ConfigDefs.hpp.
  APPEND_SET(SOURCES
    Tsqr_GlobalTimeStats.cpp
    Tsqr_MpiDatatype.cpp 
    )
  
  ASSERT_DEFINED(TSQR_ENABLE_Intel_TBB)
  # Headers for TBB-enabled shared-memory parallel TSQR.
  IF (TSQR_ENABLE_Intel_TBB)
    APPEND_SET(HEADERS 
      TbbTsqr_TbbMgs.hpp
      )
  ENDIF ()

  PACKAGE_ADD_LIBRARY(
    kokkosdisttsqr
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS kokkos kokkosnodetsqr
    )

  # Our tests test for both performance and numerical accuracy.  They
  # also function as general benchmarks, given the right command-line
  # arguments.
  PACKAGE_ADD_TEST_DIRECTORIES(test)

ENDIF () # ${PACKAGE_NAME}_ENABLE_TSQR

