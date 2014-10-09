# Make sure we added this option to the build
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TSQR)

IF(${PACKAGE_NAME}_ENABLE_TSQR)

  SET(HEADERS "")
  SET(SOURCES "")
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

  # Add TSQR-specific CMake files to the CMake module path.  Save the
  # old module path so we can restore it afterwards.
  set (CMAKE_MODULE_PATH_SAVE ${CMAKE_MODULE_PATH})
  set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

  # Generate a header file (Tsqr_Config.hpp) with the appropriate
  # #defines.  Tsqr_Config.hpp.in is in the cmake/ subdirectory, which
  # was added to the CMAKE_MODULE_PATH above.  Putting the header in
  # CMAKE_CURRENT_BINARY_DIR will ensure correct installation of the
  # generated header file.
  configure_file ("Tsqr_Config.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/Tsqr_Config.hpp")
  # Teuchos/src/CMakeLists.txt shows the canonical way to add a
  # generated header file (which doesn't start out in the source tree)
  # to HEADERS.  Putting the header in CMAKE_CURRENT_BINARY_DIR will
  # ensure correct installation of the generated header file.
  APPEND_SET(HEADERS
    ${CMAKE_CURRENT_BINARY_DIR}/Tsqr_Config.hpp)

  # Non-generated header files
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  APPEND_SET(HEADERS
    Tsqr_ApplyType.hpp
    Tsqr_CacheBlocker.hpp
    Tsqr_CacheBlockingStrategy.hpp
    Tsqr_CombineDefault.hpp
    Tsqr_Combine.hpp
    Tsqr_CombineBenchmark.hpp
    Tsqr_CombineBenchmarker.hpp
    Tsqr_CombineNative.hpp
    Tsqr_CombineTest.hpp
    Tsqr_ConfigDefs.hpp
    Tsqr_KokkosNodeTsqr.hpp
    Tsqr_KokkosNodeTsqrTest.hpp
    Tsqr_LocalVerify.hpp
    Tsqr_Matrix.hpp
    Tsqr_MatView.hpp
    Tsqr_nodeTestProblem.hpp
    Tsqr_NodeTsqr.hpp
    Tsqr_NodeTsqrFactory.hpp
    Tsqr_Random_MatrixGenerator.hpp
    Tsqr_Random_NormalGenerator.hpp
    Tsqr_ScalarTraits.hpp
    Tsqr_SeqTest.hpp
    Tsqr_SequentialCholeskyQR.hpp
    Tsqr_SequentialTsqr.hpp
    Tsqr_StatTimeMonitor.hpp
    Tsqr_TbbTest.hpp
    Tsqr_TimeStats.hpp
    Tsqr_TrivialTimer.hpp
    Tsqr_Util.hpp
    Tsqr_verifyTimerConcept.hpp
    )
  # Headers for TBB-enabled shared-memory parallel TSQR.
  IF (KokkosClassic_ENABLE_TSQR_Intel_TBB)
    APPEND_SET(HEADERS
      TbbTsqr_ApplyTask.hpp
      TbbTsqr_CacheBlockTask.hpp
      TbbTsqr_ExplicitQTask.hpp
      TbbTsqr_FactorTask.hpp
      TbbTsqr_FillWithZerosTask.hpp
      TbbTsqr.hpp
      TbbTsqr_Partitioner.hpp
      TbbTsqr_RevealRankTask.hpp
      TbbTsqr_TbbParallelTsqr.hpp
      TbbTsqr_TbbRecursiveTsqr_Def.hpp
      TbbTsqr_TbbRecursiveTsqr.hpp
      TbbTsqr_UnCacheBlockTask.hpp
      )
  ENDIF ()
  # Non-generated source files.  Does not include source files
  # corresponding to complex arithmetic routines; those are below.  We
  # have to separate them out because Trilinos does not build with
  # complex arithmetic support by default.
  APPEND_SET(SOURCES
    Tsqr_ApplyType.cpp 
    Tsqr_CombineTest.cpp
    Tsqr_SeqTest.cpp
    Tsqr_StatTimeMonitor.cpp
    Tsqr_TimeStats.cpp
    Tsqr_TrivialTimer.cpp
    )

  TRIBITS_ADD_LIBRARY(
    kokkosnodetsqr
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS kokkos
    )

  # Restore original CMake module path
  set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH_SAVE})

  # Our tests test for both performance and numerical accuracy.  They
  # also function as general benchmarks, given the right command-line
  # arguments.
  TRIBITS_ADD_TEST_DIRECTORIES(test)

ENDIF () # ${PACKAGE_NAME}_ENABLE_TSQR
