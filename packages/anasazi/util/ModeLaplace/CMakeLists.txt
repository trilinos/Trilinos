# $Header$

# Setup direct dependencies
SET(DIRECT_DEPS Anasazi)

SET(CURRENT_PACKAGE ANASAZI_MODELAPLACE)
CONFIGURE_FILE("${TRILINOS_HOME_DIR}/packages/Trilinos_DLLExportMacro.h.in"
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_ModeLaplace_DLLExportMacro.h)

# Setup directories based on current directory
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Setup directories based on dependency directories
FOREACH(DIRECT_DEP ${DIRECT_DEPS})
  INCLUDE_DIRECTORIES(${${DIRECT_DEP}_INCLUDE_DIRS})
  LINK_DIRECTORIES(${${DIRECT_DEP}_LIBRARY_DIRS})
ENDFOREACH(DIRECT_DEP ${DIRECT_DEPS})

# Set headers
SET(HEADERS
  BlockPCGSolver.h 
	CheckingTools.h 
	ModeLaplace1DQ1.h 
	ModeLaplace1DQ2.h 
	ModeLaplace2DQ1.h 
	ModeLaplace2DQ2.h 
	ModeLaplace3DQ1.h 
	ModeLaplace3DQ2.h 
	ModalProblem.h 
	ModeLaplace.h 
	SortingTools.h
  )

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_ModeLaplace_DLLExportMacro.h
  )

# Set headers
SET(SOURCES
	BlockPCGSolver.cpp 
	CheckingTools.cpp 
	ModeLaplace1DQ1.cpp 
	ModeLaplace1DQ2.cpp 
	ModeLaplace2DQ1.cpp 
	ModeLaplace2DQ2.cpp 
	ModeLaplace3DQ1.cpp 
	ModeLaplace3DQ2.cpp 
	SortingTools.cpp
  )

# Target library setup
ADD_LIBRARY(anasazi-mode-laplace ${HEADERS} ${SOURCES} -DANASAZI_MODELAPLACE_LIB_EXPORTS_MODE)
LIST(APPEND TARGET_LIBS anasazi-mode-laplace)
LIST(APPEND LIB_DEPS epetra teuchos)
TARGET_LINK_LIBRARIES(anasazi-mode-laplace ${LIB_DEPS})

# Define current package dependency variables
GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT INCLUDE_DIRECTORIES)
SET(${PROJECT_NAME}_INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_DIRS} ${INCLUDE_DIRS_CURRENT} CACHE INTERNAL "${PROJECT_NAME} include directories")
GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT LINK_DIRECTORIES)
SET(${PROJECT_NAME}_LIBRARY_DIRS ${${PROJECT_NAME}_LIBRARY_DIRS} ${LIBRARY_DIRS_CURRENT} CACHE INTERNAL "${PROJECT_NAME} library directories")
SET(${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_LIBRARIES} ${TARGET_LIBS} ${LIB_DEPS} CACHE INTERNAL "${PROJECT_NAME} libraries")
