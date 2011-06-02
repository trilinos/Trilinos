
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

SET(HAVE_EPETRA ON)
SET(HAVE_ISORROPIA_TPETRA OFF)
SET(HAVE_ISORROPIA_ZOLTAN ON)
SET(HAVE_EPETRAEXT ${PACKAGE_NAME}_ENABLE_EpetraExt)

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS)
SET(SOURCES)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  Isorropia_CostDescriber.hpp
  Isorropia_Operator.hpp
  Isorropia_Partitioner.hpp
  Isorropia_Colorer.hpp
  Isorropia_Orderer.hpp
  Isorropia_Redistributor.hpp
  Isorropia_Utils.hpp
  Isorropia_Exception.hpp
  Isorropia_ConfigDefs.hpp
  Isorropia_Version.hpp
  Isorropia_LevelScheduler.hpp
  )

SET(SOURCES ${SOURCES}
  Isorropia_Utils.cpp
  Isorropia_Exception.cpp
  )

#
# src/epetra
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/epetra)
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../utils)

SET(HEADERS ${HEADERS}
  epetra/Isorropia_EpetraCostDescriber.hpp
  epetra/Isorropia_EpetraPartitioner.hpp
  epetra/Isorropia_EpetraOperator.hpp
  epetra/Isorropia_EpetraColorer.hpp
  epetra/Isorropia_EpetraOrderer.hpp
  epetra/Isorropia_EpetraProber.hpp
  epetra/Isorropia_Epetra.hpp
  epetra/Isorropia_EpetraRedistributor.hpp
  epetra/Isorropia_EpetraLibrary.hpp
  epetra/Isorropia_EpetraLevelScheduler.hpp
  epetra/Isorropia_EpetraMatcher.hpp
  )

SET(SOURCES ${SOURCES}
  epetra/Isorropia_EpetraCostDescriber.cpp
  epetra/Isorropia_EpetraOperator.cpp
  epetra/Isorropia_EpetraPartitioner.cpp
  epetra/Isorropia_EpetraColorer.cpp
  epetra/Isorropia_EpetraOrderer.cpp
  epetra/Isorropia_EpetraProber.cpp
  epetra/Isorropia_Epetra.cpp
  epetra/Isorropia_EpetraRedistributor.cpp
  epetra/Isorropia_EpetraLibrary.cpp
  epetra/Isorropia_EpetraLevelScheduler.cpp
  epetra/Isorropia_EpetraMatcher.cpp
  )

#
# src/tpetra
#

IF (${PACKAGE_NAME}_ENABLE_TPETRA)
#   SET(HAVE_ISORROPIA_TPETRA ON)

   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/tpetra)

   SET(HEADERS ${HEADERS}
       tpetra/Isorropia_TpetraCostDescriber.hpp
       tpetra/Isorropia_TpetraOperator.hpp
       tpetra/Isorropia_TpetraPartitioner.hpp
       tpetra/Isorropia_Tpetra.hpp
       tpetra/Isorropia_TpetraRedistributor.hpp
       tpetra/Isorropia_TpetraLibrary.hpp
       tpetra/Isorropia_TpetraLevelScheduler.hpp
      )

   SET(SOURCES ${SOURCES}
       tpetra/Isorropia_TpetraCostDescriber.cpp
       tpetra/Isorropia_TpetraOperator.cpp
       tpetra/Isorropia_TpetraPartitioner.cpp
       tpetra/Isorropia_Tpetra.cpp
       tpetra/Isorropia_TpetraRedistributor.cpp
       tpetra/Isorropia_TpetraLibrary.cpp
       tpetra/Isorropia_TpetraLevelScheduler.cpp
      )

   # Zoltan
   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/tpetra/zoltanLib)

   SET(HEADERS ${HEADERS}
       tpetra/zoltanLib/Isorropia_TpetraQueryObject.hpp
       tpetra/Isorropia_TpetraZoltanLib.hpp
      )

   SET(SOURCES ${SOURCES}
       tpetra/zoltanLib/Isorropia_TpetraQueryObject.cpp
       tpetra/Isorropia_TpetraZoltanLib.cpp
      )



ENDIF()




# zoltan

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/epetra/zoltanLib)

SET(HEADERS ${HEADERS}
  epetra/zoltanLib/QueryObject.hpp
  epetra/Isorropia_EpetraZoltanLib.hpp
  #epetra/zoltanLib/Isorropia_Zoltan_Repartition.hpp
  )

SET(SOURCES ${SOURCES}
  epetra/zoltanLib/QueryObject.cpp
  epetra/Isorropia_EpetraZoltanLib.cpp
  #epetra/zoltanLib/Isorropia_Zoltan_Repartition.cpp
  )



#
# 2D partitioning
#

IF (${PACKAGE_NAME}_ENABLE_2DPARTITION)

  APPEND_SET(HEADERS
    Isorropia_Partitioner2D.hpp
    epetra/Isorropia_EpetraPartitioner2D.hpp
    )

  APPEND_SET(SOURCES
    epetra/Isorropia_EpetraPartitioner2D.cpp
    )

ENDIF()



#
# C) Define the targets for package's library(s)
#
PACKAGE_ADD_LIBRARY(
  isorropia
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#PACKAGE_ADD_LIBRARY(
#  isorropia
#  HEADERS ${HEADERS}
#  SOURCES ${SOURCES}
#  DEPLIBS isorropia_test_utils
#  )
