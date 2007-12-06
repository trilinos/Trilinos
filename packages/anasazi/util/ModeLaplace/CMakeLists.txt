	
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

INCLUDE_DIRECTORIES(${Anasazi_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${Anasazi_BINARY_DIR}/src)
INCLUDE_DIRECTORIES(${Teuchos_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${Teuchos_BINARY_DIR}/src)
INCLUDE_DIRECTORIES(${Epetra_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${Epetra_BINARY_DIR}/src)

ADD_LIBRARY(anasazi-mode-laplace ${HEADERS} ${SOURCES})

TARGET_LINK_LIBRARIES(anasazi-mode-laplace epetra teuchos)

