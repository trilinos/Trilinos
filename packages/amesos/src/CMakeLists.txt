

#
# A) Package-specific configuration options
#

#SET(HAVE_BOOL ON)
#SET(HAVE_STL ON)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# support files
#

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SOURCES ${SOURCES}
  Amesos_Status.cpp
  Amesos_Control.cpp
  Amesos_Support.cpp
  )

#
# DSCPACK files
#

IF (${PACKAGE_NAME}_ENABLE_DSCPACK)
  SET(HEADERS ${HEADERS}
    Amesos_Dscpack.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Dscpack.cpp
    )
ENDIF()

#
# UMFPACK files
#

IF (${PACKAGE_NAME}_ENABLE_UMFPACK)
  SET(HEADERS ${HEADERS}
    Amesos_Umfpack.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Umfpack.cpp
    )
ENDIF()

#
# MUMPS files
#

IF (${PACKAGE_NAME}_ENABLE_MUMPS)
   
   IF(TPL_ENABLE_MPI)


     SET(HEADERS ${HEADERS}
     Amesos_Mumps.h
     )

     SET(SOURCES ${SOURCES}
     Amesos_Mumps.cpp
     )
  ELSE()

     MESSAGE(WARNING "\n\n***Amesos with MUMPS requires MPI, MUMPS will not be available***\n\n")
   
   ENDIF()     

ENDIF()

#
# PARDISO files
#

IF (${PACKAGE_NAME}_ENABLE_PARDISO OR ${PACKAGE_NAME}_ENABLE_PARDISO_MKL)
  SET(HEADERS ${HEADERS}
    Amesos_Pardiso.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Pardiso.cpp
    )
ENDIF()

#
# MKL CSS files
#

IF (${PACKAGE_NAME}_ENABLE_CSS_MKL)
  SET(HEADERS ${HEADERS}
    Amesos_CssMKL.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_CssMKL.cpp
    )
ENDIF()

#
# CSparse files
#

IF (${PACKAGE_NAME}_ENABLE_CSparse OR ${PACKAGE_NAME}_ENABLE_CSparse)
  SET(HEADERS ${HEADERS}
    Amesos_CSparse.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_CSparse.cpp
    )
ENDIF()

#
# TAUCS files
#

IF (${PACKAGE_NAME}_ENABLE_TAUCS)
  SET(HEADERS ${HEADERS}
    Amesos_Taucs.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Taucs.cpp
    )
ENDIF()

#
# KLU files
#

IF (${PACKAGE_NAME}_ENABLE_KLU)
    SET(HEADERS ${HEADERS}
   Amesos_Klu.h
    )

    SET(SOURCES ${SOURCES}
    Amesos_Klu.cpp
    )

ENDIF()

#
# MC64 files
#

IF (${PACKAGE_NAME}_ENABLE_MC64)
  SET(HEADERS ${HEADERS}
    Amesos_MC64.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_MC64.cpp
    )
ENDIF()

#
# SCALAPACK files
#

IF (${PACKAGE_NAME}_ENABLE_SCALAPACK)
  SET(HEADERS ${HEADERS}
    Amesos_Scalapack.h
    Amesos_SCALAPACK_wrappers.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Scalapack.cpp
    )
ENDIF()

#
# FACTORY files
#

SET(HEADERS ${HEADERS}
  Amesos.h
  )

SET(SOURCES ${SOURCES}
  Amesos.cpp
  )

#
# SUPERLU files
#

IF (${PACKAGE_NAME}_ENABLE_SuperLU)
  SET(HEADERS ${HEADERS}
    Amesos_Superlu.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Superlu.cpp
    )
ENDIF()

#
# SUPERLUDIST files
#

IF (${PACKAGE_NAME}_ENABLE_SuperLUDist)
  SET(HEADERS ${HEADERS}
    Amesos_Superludist.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Superludist.cpp
    )
ENDIF()

#
# COMP files
#

SET(HEADERS ${HEADERS}
  Amesos_BaseSolver.h
  Amesos_ComponentBaseSolver.h
  Amesos_ConfigDefs.h
  CrsMatrixTranspose.h
  Amesos_NoCopiable.h
  Amesos_Utils.h
  Amesos_Time.h
  Amesos_Control.h
  Amesos_Support.h
  Amesos_Status.h
  Amesos_Scaling.h
  Amesos_Reordering.h
  )

SET(SOURCES ${SOURCES}
  CrsMatrixTranspose.cpp  
  )

IF (${PACKAGE_NAME}_ENABLE_LAPACK)
  SET(HEADERS ${HEADERS}
    Amesos_Lapack.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Lapack.cpp
    )
ENDIF()



##############################

#
# src
#

SET(HEADERS ${HEADERS}
  Amesos_BaseSolver.h
  Amesos_ComponentBaseSolver.h
  Amesos_TestRowMatrix.h
  )

SET(SOURCES ${SOURCES}

  )

IF (TPL_ENABLE_MPI)

  SET(SOURCES ${SOURCES}

    )

ENDIF()






#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  amesos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
