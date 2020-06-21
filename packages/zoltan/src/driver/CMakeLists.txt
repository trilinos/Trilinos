
SET(COMMON_DRIVER_SRC
  dr_err.c
  dr_input.c
  dr_par_util.c
  dr_util.c
  dr_elem.c
  dr_setfixed.c
  dr_output.c
  dr_chaco_io.c
  dr_random_io.c
  dr_elem_util.c
  dr_eval.c
  dr_gnuplot.c
  dr_hg_io.c
  dr_hg_readfile.c
  dr_mm_readfile.c
  dr_mmio.c
  ${PACKAGE_SOURCE_DIR}/src/ch/ch_dist_graph.c
  ${PACKAGE_SOURCE_DIR}/src/ch/ch_input_assign.c
  ${PACKAGE_SOURCE_DIR}/src/ch/ch_input_graph.c
  ${PACKAGE_SOURCE_DIR}/src/ch/ch_init_dist.c
  ${PACKAGE_SOURCE_DIR}/src/ch/ch_input_geom.c
  ${PACKAGE_SOURCE_DIR}/src/ch/ch_input_read.c
  )

IF (${PACKAGE_NAME}_ENABLE_Zlib)

  APPEND_SET(COMMON_DRIVER_SRC
  dr_compress.c
    )

ENDIF()


SET(COMMON_DRIVER_HDRS
  dr_compress_const.h
  dr_hg_readfile.h
  dr_const.h
  dr_externs.h
  dr_input_const.h
  dr_util_const.h
  dr_err_const.h
  dr_par_util_const.h
  dr_elem_const.h
  dr_elem_util_const.h
  dr_maps_const.h
  dr_eval_const.h
  dr_loadbal_const.h
  dr_mmio.h
  dr_output_const.h
  dr_param_file.h
  dr_dd.h
  ${PACKAGE_SOURCE_DIR}/src/ch/ch_input_const.h
  ${PACKAGE_SOURCE_DIR}/src/ch/ch_init_dist_const.h
  )

SET(ZDRIVE_SOURCES
  ${COMMON_DRIVER_SRC}
  dr_main.c
  dr_exoII_io.c
  dr_loadbal.c
  dr_migrate.c
  dr_maps.c
  dr_dd.c
  dr_param_file.c
  order_test.c
  )

SET(ZCPPDRIVE_SOURCES
  ${COMMON_DRIVER_SRC}
  dr_mainCPP.cpp
  dr_exoII_ioCPP.cpp
  dr_loadbalCPP.cpp
  dr_migrateCPP.cpp
  dr_mapsCPP.cpp
  dr_ddCPP.cpp
  dr_param_fileCPP.cpp
  )


TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/src/ch)

TRIBITS_ADD_EXECUTABLE(
  zdrive
  NOEXEPREFIX
  SOURCES ${COMMON_DRIVER_HDRS} ${ZDRIVE_SOURCES}
  COMM serial mpi
  CATEGORIES BASIC PERFORMANCE
  )

IF (${PACKAGE_NAME}_ENABLE_CPPDRIVER)
  TRIBITS_ADD_EXECUTABLE(
    zCPPdrive
    NOEXEPREFIX
    SOURCES ${COMMON_DRIVER_HDRS} ${ZCPPDRIVE_SOURCES}
    COMM serial mpi
    )

ENDIF()

