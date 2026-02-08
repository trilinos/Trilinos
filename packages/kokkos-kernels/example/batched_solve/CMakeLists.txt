kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR})
kokkoskernels_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

kokkoskernels_add_executable(static_pivoting SOURCES static_pivoting.cpp)

kokkoskernels_add_executable(team_GMRES SOURCES team_GMRES.cpp)

kokkoskernels_add_executable(serial_pttrs SOURCES serial_pttrs.cpp)

kokkoskernels_add_executable(serial_pbtrs SOURCES serial_pbtrs.cpp)

kokkoskernels_add_executable(serial_gbtrs SOURCES serial_gbtrs.cpp)

kokkoskernels_add_executable(serial_getrs SOURCES serial_getrs.cpp)

kokkoskernels_add_executable(serial_ger SOURCES serial_ger.cpp)

kokkoskernels_add_executable(serial_syr SOURCES serial_syr.cpp)

kokkoskernels_add_executable(serial_tbsv SOURCES serial_tbsv.cpp)

kokkoskernels_add_executable(serial_iamax SOURCES serial_iamax.cpp)

kokkoskernels_add_executable(serial_laswp SOURCES serial_laswp.cpp)

kokkoskernels_add_executable(serial_lacgv SOURCES serial_lacgv.cpp)

kokkoskernels_add_executable(serial_trsv SOURCES serial_trsv.cpp)

kokkoskernels_add_executable(serial_trsm SOURCES serial_trsm.cpp)

kokkoskernels_add_executable(serial_axpy SOURCES serial_axpy.cpp)
