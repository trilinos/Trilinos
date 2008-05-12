# $Header$ 

IF(TEUCHOS_ENABLE_COMPLEX)
  SET(HAVE_TEUCHOS_COMPLEX TRUE)
ENDIF(TEUCHOS_ENABLE_COMPLEX)

CONFIGURE_FILE(${${PROJECT_NAME}_SOURCE_DIR}/cmake/${PROJECT_NAME}_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_config.h)
INSTALL(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_config.h
  DESTINATION ${TRILINOS_INSTALL_LIB_INCLUDE_DIR}
  COMPONENT ${PROJECT_NAME}
  )

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

SET(SOURCES
	AnasaziEpetraAdapter.cpp 
	AnasaziVersion.cpp
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${Teuchos_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${Teuchos_BINARY_DIR}/src)
INCLUDE_DIRECTORIES(${Epetra_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${Epetra_BINARY_DIR}/src)
LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

ADD_LIBRARY(anasazi ${HEADERS} ${SOURCES})
LIST(APPEND LIBRARIES anasazi)

TARGET_LINK_LIBRARIES(anasazi epetra teuchos)

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

# Info for ${PROJECT_NAME}Config.cmake
GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT INCLUDE_DIRECTORIES)
SET(${PROJECT_NAME}_INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_DIRS} ${INCLUDE_DIRS_CURRENT} CACHE INTERNAL "${PROJECT_NAME} include directories")
GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT LINK_DIRECTORIES)
SET(${PROJECT_NAME}_LIBRARY_DIRS ${${PROJECT_NAME}_LIBRARY_DIRS}  ${LIBRARY_DIRS_CURRENT} CACHE INTERNAL "${PROJECT_NAME} library directories")
SET(${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_LIBRARIES} ${LIBRARIES} CACHE INTERNAL "${PROJECT_NAME} libraries")
