INCLUDE(PackageAddExecutable)
SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

ADD_DEFINITIONS(-DNO_EXODUSII)

SET(HEADERS "")
SET(SOURCES builtin.c  main.c     stat.c    getopt.c  getopt1.c
	    y.tab.c    util.c     init.c    hash.c    liballoc.c
	    add_to_log.c  units.c getline.c getdelim.c)

# If I use the PACKAGE_ADD_EXECUTABLE, it adds all libraries in this
# package to the link line for aprepro.  This includes aprepro_lib
# and we get duplicat symbols.  Use the ADD_EXECUTABLE rule instead.
# PACKAGE_ADD_EXECUTABLE(aprepro NOEXEPREFIX NOEXESUFFIX SOURCES ${SOURCES} )

ADD_EXECUTABLE(aprepro ${SOURCES})
install_executable(aprepro)

