INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)

TRIBITS_SUBPACKAGE(Aprepro)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_BINARY_DIR}"
)

ADD_DEFINITIONS(-DNO_EXODUSII)

SET(HEADERS "")
SET(SOURCES ap_builtin.c  ap_stat.c  y.tab.c  ap_util.c  ap_init.c
	    ap_hash.c     ap_liballoc.c  ap_units.c)

TRIBITS_ADD_LIBRARY(aprepro_app_lib SOURCES ${SOURCES} DEPLIBS suplib)
TRIBITS_ADD_EXECUTABLE(aprepro NOEXEPREFIX NOEXESUFFIX SOURCES ap_aprepro.c DEPLIBS aprepro_app_lib suplib m COMM serial mpi)
install_executable(aprepro)

TRIBITS_SUBPACKAGE_POSTPROCESS()

