
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

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
  Ifpack_ConfigDefs.h
  Ifpack_Version.h
  Ifpack_IlukGraph.h
  Ifpack_CrsRiluk.h
  Ifpack_ILU.h
  Ifpack_IKLU.h
  Ifpack_IKLU_Utils.h
  Ifpack_ILUT.h
  Ifpack_OverlapGraph.h
  Ifpack_ScalingType.h
  Ifpack_CrsIct.h
  Ifpack_ICT.h
  Ifpack_IC.h
  Ifpack_IC_Utils.h
  icrout_cholesky_mex.h
  ifp_parameters.h	
  Ifpack_Preconditioner.h
  Ifpack_CondestType.h
  Ifpack_Condest.h
  Ifpack.h
  Ifpack_DiagPreconditioner.h
  Ifpack_PointRelaxation.h
  Ifpack_BlockRelaxation.h
  Ifpack_Container.h
  Ifpack_SparseContainer.h
  Ifpack_DenseContainer.h
  Ifpack_AdditiveSchwarz.h
  Ifpack_LocalFilter.h
  Ifpack_NodeFilter.h
  Ifpack_DropFilter.h
  Ifpack_SparsityFilter.h
  Ifpack_ReorderFilter.h
  Ifpack_SingletonFilter.h
  Ifpack_DiagonalFilter.h
  Ifpack_Graph.h
  Ifpack_Graph_Epetra_RowMatrix.h
  Ifpack_Graph_Epetra_CrsGraph.h
  Ifpack_Partitioner.h
  Ifpack_OverlappingPartitioner.h
  Ifpack_LinearPartitioner.h
  Ifpack_GreedyPartitioner.h
  Ifpack_METISPartitioner.h
  Ifpack_EquationPartitioner.h
  Ifpack_UserPartitioner.h
  Ifpack_Utils.h
  Ifpack_HashTable.h
  Ifpack_Reordering.h
  Ifpack_RCMReordering.h
  Ifpack_SPARSKIT.h
  Ifpack_Chebyshev.h
  Ifpack_METISReordering.h
  Ifpack_OverlappingRowMatrix.h
  Ifpack_ValidParameters.h
  )

SET(SOURCES ${SOURCES}
  Ifpack_IlukGraph.cpp
  Ifpack_CrsRiluk.cpp
  Ifpack_ILU.cpp
  Ifpack_IKLU.cpp
  Ifpack_IKLU_Utils.cpp
  Ifpack_ILUT.cpp
  Ifpack_CrsIct.cpp
  Ifpack_ICT.cpp
  Ifpack_IC.cpp
  Ifpack_IC_Utils.cpp
  Ifpack_OverlapGraph.cpp
  icrout_cholesky_mex.c
  icrout_quicksort.c
  ifp_parameters.cpp
  Ifpack_DiagPreconditioner.cpp
  Ifpack_PointRelaxation.cpp
  Ifpack_Condest.cpp
  Ifpack.cpp
  Ifpack_DenseContainer.cpp
  Ifpack_LocalFilter.cpp
  Ifpack_NodeFilter.cpp
  Ifpack_DropFilter.cpp
  Ifpack_SparsityFilter.cpp
  Ifpack_ReorderFilter.cpp
  Ifpack_SingletonFilter.cpp
  Ifpack_DiagonalFilter.cpp
  Ifpack_Graph_Epetra_RowMatrix.cpp
  Ifpack_Graph_Epetra_CrsGraph.cpp
  Ifpack_Utils.cpp
  Ifpack_Chebyshev.cpp
  Ifpack_OverlappingPartitioner.cpp
  Ifpack_LinearPartitioner.cpp
  Ifpack_GreedyPartitioner.cpp
  Ifpack_METISPartitioner.cpp
  Ifpack_EquationPartitioner.cpp
  Ifpack_UserPartitioner.cpp
  Ifpack_RCMReordering.cpp
  Ifpack_METISReordering.cpp
  Ifpack_OverlappingRowMatrix.cpp
  Ifpack_ValidParameters.cpp
  )

IF (${PACKAGE_NAME}_ENABLE_Amesos)

  INCLUDE_DIRECTORIES(  
    ${PACKAGE_SOURCE_DIR}/../amesos/src/SuiteSparse/AMD/Include
    ${PACKAGE_SOURCE_DIR}/../amesos/src/SuiteSparse/UFconfig
    )

  Set(HEADERS ${HEADERS}
    Ifpack_AMDReordering.h
    Ifpack_Amesos.h
    )
  
  SET(SOURCES ${SOURCES}
    Ifpack_AMDReordering.cpp
    Ifpack_Amesos.cpp
    )

ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  ifpack
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
