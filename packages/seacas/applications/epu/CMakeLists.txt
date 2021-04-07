# This will build one or two versions of epu.
# In a MPI-enabled build, it will build 'epu' which is a serial application and 'pepu' which is a parallel-enabled epu
# In a serial, non-MPI-enabled build, it will only build 'epu'.
# On some platforms, both applications may perform the same if run in serial; however, some platforms require
# parallel applications to be run on a compute-partition.  We used to just build a parallel-enabled epu
# when MPI was enabled, but this caused issues for users on those platforms which would not allow a
# parallel application to be run on a non-compute partition, so now we build two versions...

TRIBITS_SUBPACKAGE(Epu)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/EP_ExodusFile.C
		    ${CMAKE_CURRENT_SOURCE_DIR}/EP_Internals.C
		    ${CMAKE_CURRENT_SOURCE_DIR}/EP_ParallelDisks.C
		    ${CMAKE_CURRENT_SOURCE_DIR}/EP_SystemInterface.C)

#disabling on PGI due to an odd compiler issue with the compiler generated assembly.
#this error only seems to manifest on pgi 11.1, however, the version of the compiler
#is not given for pgi compilers unfortunately so we can't be more specific
IF(NOT __LINUX_COMPILER_PGI)
  TRIBITS_ADD_EXECUTABLE(
    epu
    NOEXEPREFIX
    NOEXESUFFIX
    SOURCES ${SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/epu.C
    INSTALLABLE
    )

IF (TPL_ENABLE_MPI)
  TRIBITS_ADD_EXECUTABLE(
    pepu
    NOEXEPREFIX
    NOEXESUFFIX
    SOURCES ${SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/pepu.C
    COMM mpi
    INSTALLABLE
    )
ENDIF()
ENDIF()

TRIBITS_SUBPACKAGE_POSTPROCESS()
