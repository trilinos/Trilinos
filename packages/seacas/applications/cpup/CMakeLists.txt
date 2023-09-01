# This will build one or two versions of epu.
# In a MPI-enabled build, it will build 'epu' which is a serial application and 'pepu' which is a parallel-enabled epu
# In a serial, non-MPI-enabled build, it will only build 'epu'.
# On some platforms, both applications may perform the same if run in serial; however, some platforms require
# parallel applications to be run on a compute-partition.  We used to just build a parallel-enabled epu
# when MPI was enabled, but this caused issues for users on those platforms which would not allow a
# parallel application to be run on a non-compute partition, so now we build two versions...

TRIBITS_SUBPACKAGE(Cpup)

TRIBITS_INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/CP_SystemInterface.C)

ASSERT_DEFINED(TPL_ENABLE_CGNS)
IF (TPL_ENABLE_CGNS)
  TRIBITS_ADD_EXECUTABLE(
    cpup
    NOEXEPREFIX
    NOEXESUFFIX
    SOURCES ${SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/cpup.C
    INSTALLABLE
    )
ENDIF()

TRIBITS_SUBPACKAGE_POSTPROCESS()
