IF (Tacho_ENABLE_INT_INT)
  SET(TACHO_USE_INT_INT ON)
  SET(TACHO_USE_INT_SIZE_T OFF)
ELSE()
  SET(TACHO_USE_INT_INT OFF)
  SET(TACHO_USE_INT_SIZE_T ON)
ENDIF()

SET(TACHO_HAVE_MKL ${TPL_ENABLE_MKL})
SET(TACHO_HAVE_LAPACK ${TPL_ENABLE_LAPACK})
SET(TACHO_HAVE_BLAS ${TPL_ENABLE_BLAS})

SET(TACHO_HAVE_METIS ${TPL_ENABLE_METIS})
SET(TACHO_HAVE_METIS_MT ${TPL_ENABLE_MTMETIS})
SET(TACHO_HAVE_SCOTCH ${TPL_ENABLE_Scotch})
SET(TACHO_HAVE_TRILINOS_SS ${Trilinos_ENABLE_TrilinosSS})
SET(TACHO_HAVE_SUITESPARSE ${TPL_ENABLE_Cholmod})
                       
SET(TACHO_HAVE_VTUNE ${TPL_ENABLE_VTune})

#-----------------------------------------------------------------------------

TRIBITS_CONFIGURE_FILE(Tacho_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/impl)
#-----------------------------------------------------------------------------

SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#-----------------------------------------------------------------------------

SET(HEADERS_PUBLIC "")
SET(HEADERS_PRIVATE "")
SET(SOURCES "")

FILE(GLOB HEADERS_PUBLIC *.hpp impl/*.hpp)
LIST( APPEND HEADERS_PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/Tacho_config.h )


APPEND_SET(SOURCES
  impl/Tacho_Util.cpp
  impl/Tacho_Blas_External.cpp
  impl/Tacho_Lapack_External.cpp
  impl/Tacho_GraphTools_Metis.cpp
  impl/Tacho_SymbolicTools.cpp
  eti/Tacho_Solver_ETI_ChseLevTaskScheduler_complex_double_Cuda.cpp
  eti/Tacho_Solver_ETI_ChseLevTaskScheduler_complex_double_OpenMP.cpp
  eti/Tacho_Solver_ETI_ChseLevTaskScheduler_complex_double_Serial.cpp
  eti/Tacho_Solver_ETI_ChseLevTaskScheduler_complex_float_Cuda.cpp
  eti/Tacho_Solver_ETI_ChseLevTaskScheduler_complex_float_OpenMP.cpp
  eti/Tacho_Solver_ETI_ChseLevTaskScheduler_complex_float_Serial.cpp
  eti/Tacho_Solver_ETI_ChseLevTaskScheduler_double_Cuda.cpp
  eti/Tacho_Solver_ETI_ChseLevTaskScheduler_double_OpenMP.cpp
  eti/Tacho_Solver_ETI_ChseLevTaskScheduler_double_Serial.cpp
  eti/Tacho_Solver_ETI_ChseLevTaskScheduler_float_Cuda.cpp
  eti/Tacho_Solver_ETI_ChseLevTaskScheduler_float_OpenMP.cpp
  eti/Tacho_Solver_ETI_ChseLevTaskScheduler_float_Serial.cpp
  eti/Tacho_Solver_ETI_TaskScheduler_complex_double_Cuda.cpp
  eti/Tacho_Solver_ETI_TaskScheduler_complex_double_OpenMP.cpp
  eti/Tacho_Solver_ETI_TaskScheduler_complex_double_Serial.cpp
  eti/Tacho_Solver_ETI_TaskScheduler_complex_float_Cuda.cpp
  eti/Tacho_Solver_ETI_TaskScheduler_complex_float_OpenMP.cpp
  eti/Tacho_Solver_ETI_TaskScheduler_complex_float_Serial.cpp
  eti/Tacho_Solver_ETI_TaskScheduler_double_Cuda.cpp
  eti/Tacho_Solver_ETI_TaskScheduler_double_OpenMP.cpp
  eti/Tacho_Solver_ETI_TaskScheduler_double_Serial.cpp
  eti/Tacho_Solver_ETI_TaskScheduler_float_Cuda.cpp
  eti/Tacho_Solver_ETI_TaskScheduler_float_OpenMP.cpp
  eti/Tacho_Solver_ETI_TaskScheduler_float_Serial.cpp
  eti/Tacho_Solver_ETI_TaskSchedulerMultiple_complex_double_Cuda.cpp
  eti/Tacho_Solver_ETI_TaskSchedulerMultiple_complex_double_OpenMP.cpp
  eti/Tacho_Solver_ETI_TaskSchedulerMultiple_complex_double_Serial.cpp
  eti/Tacho_Solver_ETI_TaskSchedulerMultiple_complex_float_Cuda.cpp
  eti/Tacho_Solver_ETI_TaskSchedulerMultiple_complex_float_OpenMP.cpp
  eti/Tacho_Solver_ETI_TaskSchedulerMultiple_complex_float_Serial.cpp
  eti/Tacho_Solver_ETI_TaskSchedulerMultiple_double_Cuda.cpp
  eti/Tacho_Solver_ETI_TaskSchedulerMultiple_double_OpenMP.cpp
  eti/Tacho_Solver_ETI_TaskSchedulerMultiple_double_Serial.cpp
  eti/Tacho_Solver_ETI_TaskSchedulerMultiple_float_Cuda.cpp
  eti/Tacho_Solver_ETI_TaskSchedulerMultiple_float_OpenMP.cpp
  eti/Tacho_Solver_ETI_TaskSchedulerMultiple_float_Serial.cpp
)

#-----------------------------------------------------------------------------

TRIBITS_ADD_LIBRARY(
    tacho
    HEADERS ${HEADERS_PUBLIC}
    NOINSTALLHEADERS ${HEADERS_PRIVATE}
    SOURCES ${SOURCES}
    DEPLIBS
    )
