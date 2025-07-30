

#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(NOINSTALLHEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src-rbgen)

# Abstract headers
APPEND_SET(HEADERS
    RBGen_BurkardtFileIOHandler.h
    RBGen_EpetraMVFileIOFactory.h
    RBGen_EpetraMVMethodFactory.h
    RBGen_EpetraMVPreprocessorFactory.h
    RBGen_LapackPOD.h
    RBGen_MSPreprocessor.h
    RBGen_Utils.h
    )

APPEND_SET(SOURCES
    RBGen_BurkardtFileIOHandler.cpp
    RBGen_EpetraMVFileIOFactory.cpp
    RBGen_EpetraMVMethodFactory.cpp
    RBGen_EpetraMVPreprocessorFactory.cpp
    RBGen_LapackPOD.cpp
    RBGen_MSPreprocessor.cpp
    RBGen_Utils.cpp
    )

# Anasazi-based methods
APPEND_SET(HEADERS
    RBGen_AnasaziPOD.h
    RBGen_IncSVDPOD.h
    RBGen_ISVDMultiCD.h
    RBGen_ISVD_MultiCDUDV.h
    RBGen_ISVDMultiSDA.h
    RBGen_ISVD_MultiSDAUDV.h
    RBGen_ISVDMultiSDB.h
    RBGen_ISVD_MultiSDBUDV.h
    RBGen_ISVDSingle.h
    RBGen_ISVD_SingleUDV.h
    RBGen_ISVDUDV.h
    RBGen_StSVD_RTR.h
    )

APPEND_SET(SOURCES
    RBGen_AnasaziPOD.cpp
    RBGen_IncSVDPOD.cpp
    RBGen_ISVDMultiCD.cpp
    RBGen_ISVD_MultiCDUDV.cpp
    RBGen_ISVDMultiSDA.cpp
    RBGen_ISVD_MultiSDAUDV.cpp
    RBGen_ISVDMultiSDB.cpp
    RBGen_ISVD_MultiSDBUDV.cpp
    RBGen_ISVDSingle.cpp
    RBGen_ISVD_SingleUDV.cpp
    RBGen_ISVDUDV.cpp
    RBGen_StSVD_RTR.cpp
    )



#
# C) Define the targets for package's library(s)
#
IF (${PACKAGE_NAME}_ENABLE_RBGen)
  TRIBITS_ADD_LIBRARY(
    rbgenepetra
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS rbgen anasaziepetra anasazi
    )
ENDIF()

