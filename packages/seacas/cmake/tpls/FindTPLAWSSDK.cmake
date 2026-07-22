# Copyright(C) 2025 National Technology & Engineering Solutions
# of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
# NTESS, the U.S. Government retains certain rights in this software.
#
# See packages/seacas/LICENSE for details

MESSAGE("-- Using FIND_PACKAGE(AWSSDK ...) ...")

FIND_PACKAGE(AWSSDK REQUIRED COMPONENTS s3 transfer)

IF (AWSSDK_FOUND)
  # For compatibility with TriBITS:
  SET(DOCSTR "List of semi-colon separated paths to look for the TPL AWSSDK")
  # AWSSDK gives back a list of library names not absolute pathnames.
  # Downstream projects won't link correctly.  Rewrite each library
  # name to include the library installation path.
  SET(AWSSDK_TARGETS ${AWSSDK_LIBRARIES})
  FOREACH(TARGET IN LISTS AWSSDK_TARGETS)
    MESSAGE(STATUS "Adding ${TARGET} to AWSSDK_LIB_PATHNAMES as ${AWSSDK_LIB_DIR}/lib${TARGET}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    LIST(APPEND AWSSDK_LIB_PATHNAMES "${AWSSDK_LIB_DIR}/lib${TARGET}${CMAKE_SHARED_LIBRARY_SUFFIX}")
  ENDFOREACH()
  SET(TPL_AWSSDK_LIBRARIES ${AWSSDK_LIB_PATHNAMES} CACHE PATH ${DOCSTR})
  SET(TPL_AWSSDK_INCLUDE_DIRS ${AWSSDK_INCLUDE_DIRS} CACHE PATH ${DOCSTR})
  SET(TPL_AWSSDK_LIBRARY_DIRS ${AWSSDK_LIB_DIR} CACHE PATH ${DOCSTR})
ENDIF()

#
# Third, call TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()
#
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES(
  AWSSDK
  REQUIRED_HEADERS aws/core/Aws.h aws/s3/S3Client.h aws/transfer/TransferManager.h
  REQUIRED_LIBS_NAMES aws-cpp-sdk-core aws-cpp-sdk-s3 aws-cpp-sdk-transfer
  )

# NOTE: If FIND_PACKAGE(AWSSDK ...) was called and successfully found AWSSDK, then
# TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES() will use the already-set
# variables TPL_AWSSDK_INCLUDE_DIRS and TPL_AWSSDK_LIBRARIES and then print them
# out (and set some other standard variables as well).  This is the final
# "hook" into the TriBITS TPL system.
