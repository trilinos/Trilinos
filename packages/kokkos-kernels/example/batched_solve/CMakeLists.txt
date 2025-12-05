KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

KOKKOSKERNELS_ADD_EXECUTABLE(
  static_pivoting
  SOURCES static_pivoting.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  team_GMRES
  SOURCES team_GMRES.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  serial_pttrs
  SOURCES serial_pttrs.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  serial_pbtrs
  SOURCES serial_pbtrs.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  serial_gbtrs
  SOURCES serial_gbtrs.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  serial_getrs
  SOURCES serial_getrs.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  serial_ger
  SOURCES serial_ger.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  serial_syr
  SOURCES serial_syr.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  serial_tbsv
  SOURCES serial_tbsv.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  serial_iamax
  SOURCES serial_iamax.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  serial_laswp
  SOURCES serial_laswp.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  serial_lacgv
  SOURCES serial_lacgv.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  serial_trsv
  SOURCES serial_trsv.cpp
  )
