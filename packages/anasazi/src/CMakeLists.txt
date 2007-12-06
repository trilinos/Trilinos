SET(HAVE_IOSTREAM TRUE)
SET(HAVE_ITERATOR TRUE)
SET(HAVE_NUMERIC TRUE)
SET(HAVE_SSTREAM TRUE)
SET(HAVE_STDEXCEPT TRUE)
SET(HAVE_VECTOR TRUE)

CONFIGURE_FILE(${Anasazi_SOURCE_DIR}/cmake/Anasazi_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/Anasazi_config.h)

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

ADD_LIBRARY(anasazi ${HEADERS} ${SOURCES})

TARGET_LINK_LIBRARIES(anasazi epetra teuchos)

