# $Header$ 

# Setup direct dependencies
SET(DIRECT_DEPS Teuchos Epetra)

# Package-specific configuration options
IF(TEUCHOS_ENABLE_COMPLEX)
  SET(HAVE_TEUCHOS_COMPLEX TRUE)
ENDIF(TEUCHOS_ENABLE_COMPLEX)

# Create <PACKAGE>_config.h
CONFIGURE_FILE(
  ${${PROJECT_NAME}_SOURCE_DIR}/cmake/${PROJECT_NAME}_config.h.in 
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_config.h
)
INSTALL(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_config.h
  DESTINATION ${TRILINOS_INSTALL_LIB_INCLUDE_DIR}
  COMPONENT ${PROJECT_NAME}
  )

# Setup directories based on current directory
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Setup directories based on dependency directories
FOREACH(DIRECT_DEP ${DIRECT_DEPS})
  INCLUDE_DIRECTORIES(${${DIRECT_DEP}_INCLUDE_DIRS})
  LINK_DIRECTORIES(${${DIRECT_DEP}_LIBRARY_DIRS})
ENDFOREACH()

# Set headers
SET(HEADERS
	AnasaziBasicEigenproblem.hpp 
	AnasaziBasicOrthoManager.hpp 
	AnasaziBasicOutputManager.hpp 
	AnasaziBasicSort.hpp 
	AnasaziBlockDavidson.hpp 
	AnasaziBlockDavidsonSolMgr.hpp 
	AnasaziBlockKrylovSchur.hpp 
	AnasaziBlockKrylovSchurSolMgr.hpp 
	AnasaziConfigDefs.hpp 
	AnasaziEigenproblem.hpp 
	AnasaziEigensolver.hpp 
	AnasaziEigensolverDecl.hpp 
	AnasaziEpetraAdapter.hpp 
	AnasaziLOBPCG.hpp 
	AnasaziLOBPCGSolMgr.hpp 
	AnasaziHelperTraits.hpp 
	AnasaziMatOrthoManager.hpp 
	AnasaziSolverUtils.hpp 
	AnasaziMultiVec.hpp 
	AnasaziMultiVecTraits.hpp 
	AnasaziOperator.hpp 
	AnasaziOperatorTraits.hpp 
	AnasaziOrthoManager.hpp 
	AnasaziOutputManager.hpp 
	AnasaziSimpleLOBPCGSolMgr.hpp 
	AnasaziSolverManager.hpp 
	AnasaziSortManager.hpp 
	AnasaziStatusTest.hpp 
	AnasaziStatusTestDecl.hpp 
	AnasaziStatusTestCombo.hpp 
	AnasaziStatusTestMaxIters.hpp 
	AnasaziStatusTestWithOrdering.hpp 
	AnasaziStatusTestOutput.hpp 
	AnasaziStatusTestResNorm.hpp 
	AnasaziSVQBOrthoManager.hpp 
	AnasaziTypes.hpp
  )

# Set sources
SET(SOURCES
	AnasaziEpetraAdapter.cpp 
	AnasaziVersion.cpp
  )

# Target library setup
ADD_LIBRARY(anasazi ${HEADERS} ${SOURCES})
LIST(APPEND TARGET_LIBS anasazi)
LIST(APPEND LIB_DEPS epetra teuchos)
TARGET_LINK_LIBRARIES(anasazi ${LIB_DEPS})

# MPI setup
IF(Trilinos_ENABLE_MPI)
  LIST(APPEND LIB_DEPS ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_DIRS})
  TARGET_LINK_LIBRARIES(anasazi ${LIB_DEPS})
ENDIF(Trilinos_ENABLE_MPI)

INSTALL(
  TARGETS anasazi
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT ${PROJECT_NAME}
  )

INSTALL(
  FILES ${HEADERS}
  DESTINATION ${TRILINOS_INSTALL_INCLUDE_DIR}
  COMPONENT ${PROJECT_NAME}
  )

# Define current package dependency variables
GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT INCLUDE_DIRECTORIES)
SET(${PROJECT_NAME}_INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_DIRS} ${INCLUDE_DIRS_CURRENT} CACHE INTERNAL "${PROJECT_NAME} include directories")
GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT LINK_DIRECTORIES)
SET(${PROJECT_NAME}_LIBRARY_DIRS ${${PROJECT_NAME}_LIBRARY_DIRS}  ${LIBRARY_DIRS_CURRENT} CACHE INTERNAL "${PROJECT_NAME} library directories")
SET(${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_LIBRARIES} ${TARGET_LIBS} ${LIB_DEPS} CACHE INTERNAL "${PROJECT_NAME} libraries")
