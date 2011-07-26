INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)

SUBPACKAGE(Aprepro)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

ADD_DEFINITIONS(-DNO_EXODUSII)

SET(HEADERS "")
SET(SOURCES builtin.c  main.c     stat.c    getopt.c  getopt1.c
	    y.tab.c    util.c     init.c    hash.c    liballoc.c
	    add_to_log.c  units.c getline.c getdelim.c)

PACKAGE_ADD_EXECUTABLE(aprepro NOEXEPREFIX NOEXESUFFIX SOURCES ${SOURCES} DEPLIBS m)
install_executable(aprepro)

