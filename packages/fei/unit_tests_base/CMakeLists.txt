
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../base)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(SOURCES "")

APPEND_SET(SOURCES
fei_UBase_CSRMat_CSVec.cpp
fei_UBase_CommMap.cpp
fei_UBase_CommUtils.cpp
fei_UBase_DirBC.cpp
fei_UBase_DofMapper.cpp
fei_UBase_LinearDecomposition.cpp
fei_UBase_MatrixTraits.cpp
fei_UBase_Params.cpp
fei_UBase_Reducer.cpp
fei_UBase_ReverseMapper.cpp
fei_UBase_exchangeCommMapData.cpp
fei_UBase_set_dof_mappings.cpp
fei_UBase_set_shared_ids.cpp
fei_UBase_utils.cpp
fei_ubase.cpp
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    fei_ubase
    SOURCES ${SOURCES}
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
   STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_TEST(
    fei_ubase
    NAME fei_ubase
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 2
   STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_TEST(
    fei_ubase
    NAME fei_ubase
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 3
   STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_TEST(
    fei_ubase
    NAME fei_ubase
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 4
   STANDARD_PASS_OUTPUT
  )

