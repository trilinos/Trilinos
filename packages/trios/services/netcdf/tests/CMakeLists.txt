INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddTest)

# Need to include the binary dir so we can find config.h
include_directories(${CMAKE_BINARY_DIR})

# Need to include the source dir so we can find aggregation.h
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)
# Need to include the binary dir so we can find netcdf_args.h
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../src)

LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../src)

cmake_policy(VERSION 2.8)


######################  UNIT TESTS  ##################################################

TRIBITS_ADD_EXECUTABLE(
  netcdf-remote-test
  SOURCES netcdf_test.c
  DEPLIBS ${Trios_NETCDF_LIBRARIES} ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  netcdf-test
  SOURCES netcdf_test.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  aggregation_test
  SOURCES aggregation_test.cpp ../src/aggregation.cpp ../src/netcdf_debug.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)


######################  PROGRAMS  ##################################################

TRIBITS_ADD_EXECUTABLE(
  simple_xy_wr
  SOURCES simple_xy_wr.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  simple_xy_wr_remote
  SOURCES simple_xy_wr.c
  DEPLIBS ${Trios_NETCDF_LIBRARIES} ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  simple_xy_rd
  SOURCES simple_xy_rd.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  simple_xy_rd_remote
  SOURCES simple_xy_rd.c
  DEPLIBS ${Trios_NETCDF_LIBRARIES} ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  ncmpi_xy_wr
  SOURCES ncmpi_xy_wr.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  ncmpi_xy_wr_remote
  SOURCES ncmpi_xy_wr.c
  DEPLIBS ${Trios_NETCDF_LIBRARIES} ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  ncmpi_xy_rd
  SOURCES ncmpi_xy_rd.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  ncmpi_xy_rd_remote
  SOURCES ncmpi_xy_rd.c
  DEPLIBS ${Trios_NETCDF_LIBRARIES} ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)


######### PRES_TEMP_4D EXAMPLE ####################

TRIBITS_ADD_EXECUTABLE(
  pres_temp_4D_wr
  SOURCES pres_temp_4D_wr.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  pres_temp_4D_wr_remote
  SOURCES pres_temp_4D_wr.c
  DEPLIBS ${Trios_NETCDF_LIBRARIES} ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  pres_temp_4D_rd
  SOURCES pres_temp_4D_rd.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  pres_temp_4D_rd_remote
  SOURCES pres_temp_4D_rd.c
  DEPLIBS ${Trios_NETCDF_LIBRARIES} ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  ncmpi_pres_temp_4D_wr
  SOURCES ncmpi_pres_temp_4D_wr.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  ncmpi_pres_temp_4D_wr_remote
  SOURCES ncmpi_pres_temp_4D_wr.c
  DEPLIBS ${Trios_NETCDF_LIBRARIES} ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  ncmpi_pres_temp_4D_rd
  SOURCES ncmpi_pres_temp_4D_rd.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  ncmpi_pres_temp_4D_rd_remote
  SOURCES ncmpi_pres_temp_4D_rd.c
  DEPLIBS ${Trios_NETCDF_LIBRARIES} ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)


######### SVC_PRES_TEMP EXAMPLE ####################

TRIBITS_ADD_EXECUTABLE(
  sfc_pres_temp_wr
  SOURCES sfc_pres_temp_wr.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  sfc_pres_temp_wr_remote
  SOURCES sfc_pres_temp_wr.c
  DEPLIBS ${Trios_NETCDF_LIBRARIES} ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  sfc_pres_temp_rd
  SOURCES sfc_pres_temp_rd.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  sfc_pres_temp_rd_remote
  SOURCES sfc_pres_temp_rd.c
  DEPLIBS ${Trios_NETCDF_LIBRARIES} ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

######### Integer tests from parallel-netcdf distribution ####################

TRIBITS_ADD_EXECUTABLE(
  test_read_coll
  SOURCES test_read_coll.c testutils.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  test_read_coll_remote
  SOURCES test_read_coll.c testutils.c
  DEPLIBS ${Trios_NETCDF_LIBRARIES} ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  test_read_indep
  SOURCES test_read_indep.c testutils.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  test_read_indep_remote
  SOURCES test_read_indep.c testutils.c
  DEPLIBS ${Trios_NETCDF_LIBRARIES} ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  test_write_coll
  SOURCES test_write_coll.c testutils.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  test_write_coll_remote
  SOURCES test_write_coll.c testutils.c
  DEPLIBS ${Trios_NETCDF_LIBRARIES} ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  test_write_indep
  SOURCES test_write_indep.c testutils.c
  DEPLIBS ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE(
  test_write_indep_remote
  SOURCES test_write_indep.c testutils.c
  DEPLIBS ${Trios_NETCDF_LIBRARIES} ${DEPLIBS}
  COMM serial mpi
  NOEXEPREFIX
)
