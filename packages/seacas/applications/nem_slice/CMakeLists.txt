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
SET(SOURCES elb_allo.C       elb_elem.C       elb_err.C        elb_exo_util.C
	    elb_graph.C      elb_groups.C     elb_inp.C        elb_loadbal.C
	    elb_output.C     elb_time.C       elb_util.C
	    md_getsubopt.C)

CHECK_INCLUDE_FILE( sys/resource.h HAS_GETRUSAGE )

TRIBITS_ADD_LIBRARY(nem_slice_app_lib SOURCES ${SOURCES} DEPLIBS ${ZOLTAN_DEP} exodus nemesis chaco suplib)
TRIBITS_ADD_EXECUTABLE(nem_slice NOEXEPREFIX NOEXESUFFIX SOURCES elb_main.C DEPLIBS ${ZOLTAN_DEP} nem_slice_app_lib COMM serial mpi)
install_executable(nem_slice)

TRIBITS_SUBPACKAGE_POSTPROCESS()

