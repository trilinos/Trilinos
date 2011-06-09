INCLUDE(PackageAddExecutable)
INCLUDE(CheckIncludeFile)

IF (Trilinos_ENABLE_Zoltan)
   INCLUDE_DIRECTORIES(${Zoltan_INCLUDE_DIRS})
   SET(ZOLTAN_DEP zoltan simpi)
   ADD_DEFINITIONS(-DUSE_ZOLTAN)
ENDIF()

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS "")
SET(SOURCES elb_allo.c       elb_elem.c       elb_err.c        elb_exo_util.c
	    elb_graph.c      elb_groups.c     elb_inp.c        elb_loadbal.c
	    elb_main.c       elb_output.c     elb_time.c       elb_util.c
	    add_to_log.c     getopt.c         md_getsubopt.c)

CHECK_INCLUDE_FILE( sys/resource.h HAS_GETRUSAGE )
PACKAGE_ADD_EXECUTABLE(nem_slice NOEXEPREFIX NOEXESUFFIX SOURCES ${SOURCES} DEPLIBS ${ZOLTAN_DEP} exodus nemesis COMM serial)
install_executable(nem_slice)

