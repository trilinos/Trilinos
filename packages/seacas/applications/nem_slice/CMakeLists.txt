INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsLibraryMacros)
INCLUDE(CheckIncludeFile)

TRIBITS_SUBPACKAGE(Nemslice)

IF (${PACKAGE_NAME}_ENABLE_Zoltan)
   INCLUDE_DIRECTORIES(${Zoltan_INCLUDE_DIRS})
   SET(ZOLTAN_DEP zoltan)
   ADD_DEFINITIONS(-DUSE_ZOLTAN)
ENDIF()

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS "")
SET(SOURCES elb_allo.c       elb_elem.c       elb_err.c        elb_exo_util.c
	    elb_graph.c      elb_groups.c     elb_inp.c        elb_loadbal.c
	    elb_output.c     elb_time.c       elb_util.c
	    md_getsubopt.c)

CHECK_INCLUDE_FILE( sys/resource.h HAS_GETRUSAGE )

TRIBITS_ADD_LIBRARY(nem_slice_app_lib SOURCES ${SOURCES} DEPLIBS ${ZOLTAN_DEP} exodus nemesis chaco suplib)
TRIBITS_ADD_EXECUTABLE(nem_slice NOEXEPREFIX NOEXESUFFIX SOURCES elb_main.c DEPLIBS ${ZOLTAN_DEP} nem_slice_app_lib COMM serial mpi)
install_executable(nem_slice)

TRIBITS_SUBPACKAGE_POSTPROCESS()

