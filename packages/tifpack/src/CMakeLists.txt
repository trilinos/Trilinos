
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

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
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  Tifpack_ConfigDefs.hpp
  Tifpack_Version.hpp
  Tifpack_IlukGraph.hpp
  Tifpack_CrsRiluk.hpp
  Tifpack_ILU.hpp
  Tifpack_IKLU.hpp
  Tifpack_IKLU_Utils.hpp
  Tifpack_ILUT.hpp
  Tifpack_OverlapGraph.hpp
  Tifpack_ScalingType.hpp
  Tifpack_CrsIct.hpp
  Tifpack_ICT.hpp
  Tifpack_IC.hpp
  Tifpack_IC_Utils.hpp
  icrout_cholesky_mex.hpp
  ifp_parameters.hpp	
  Tifpack_Preconditioner.hpp
  Tifpack_CondestType.hpp
  Tifpack_Condest.hpp
  Tifpack.hpp
  Tifpack_DiagPreconditioner.hpp
  Tifpack_PointRelaxation.hpp
  Tifpack_BlockRelaxation.hpp
  Tifpack_Container.hpp
  Tifpack_SparseContainer.hpp
  Tifpack_DenseContainer.hpp
  Tifpack_AdditiveSchwarz.hpp
  Tifpack_LocalFilter.hpp
  Tifpack_NodeFilter.hpp
  Tifpack_DropFilter.hpp
  Tifpack_SparsityFilter.hpp
  Tifpack_ReorderFilter.hpp
  Tifpack_SingletonFilter.hpp
  Tifpack_DiagonalFilter.hpp
  Tifpack_Graph.hpp
  Tifpack_Graph_Tpetra_RowMatrix.hpp
  Tifpack_Graph_Tpetra_CrsGraph.hpp
  Tifpack_Partitioner.hpp
  Tifpack_OverlappingPartitioner.hpp
  Tifpack_LinearPartitioner.hpp
  Tifpack_GreedyPartitioner.hpp
  Tifpack_METISPartitioner.hpp
  Tifpack_EquationPartitioner.hpp
  Tifpack_UserPartitioner.hpp
  Tifpack_Utils.hpp
  Tifpack_HashTable.hpp
  Tifpack_Reordering.hpp
  Tifpack_RCMReordering.hpp
  Tifpack_SPARSKIT.hpp
  Tifpack_Chebyshev.hpp
  Tifpack_METISReordering.hpp
  Tifpack_OverlappingRowMatrix.hpp
  Tifpack_ValidParameters.hpp
  )

SET(SOURCES ${SOURCES}
  Tifpack_IlukGraph.cpp
  Tifpack_CrsRiluk.cpp
  Tifpack_ILU.cpp
  Tifpack_IKLU.cpp
  Tifpack_IKLU_Utils.cpp
  Tifpack_ILUT.cpp
  Tifpack_CrsIct.cpp
  Tifpack_ICT.cpp
  Tifpack_IC.cpp
  Tifpack_IC_Utils.cpp
  Tifpack_OverlapGraph.cpp
  icrout_cholesky_mex.c
  icrout_quicksort.c
  ifp_parameters.cpp
  Tifpack_DiagPreconditioner.cpp
  Tifpack_PointRelaxation.cpp
  Tifpack_Condest.cpp
  Tifpack.cpp
  Tifpack_DenseContainer.cpp
  Tifpack_LocalFilter.cpp
  Tifpack_NodeFilter.cpp
  Tifpack_DropFilter.cpp
  Tifpack_SparsityFilter.cpp
  Tifpack_ReorderFilter.cpp
  Tifpack_SingletonFilter.cpp
  Tifpack_DiagonalFilter.cpp
  Tifpack_Graph_Tpetra_RowMatrix.cpp
  Tifpack_Graph_Tpetra_CrsGraph.cpp
  Tifpack_Utils.cpp
  Tifpack_Chebyshev.cpp
  Tifpack_OverlappingPartitioner.cpp
  Tifpack_LinearPartitioner.cpp
  Tifpack_GreedyPartitioner.cpp
  Tifpack_METISPartitioner.cpp
  Tifpack_EquationPartitioner.cpp
  Tifpack_UserPartitioner.cpp
  Tifpack_RCMReordering.cpp
  Tifpack_METISReordering.cpp
  Tifpack_OverlappingRowMatrix.cpp
  Tifpack_ValidParameters.cpp
  )

IF (${PACKAGE_NAME}_ENABLE_Amesos)

  INCLUDE_DIRECTORIES(  
    ${PACKAGE_SOURCE_DIR}/../amesos/src/SuiteSparse/AMD/Include
    ${PACKAGE_SOURCE_DIR}/../amesos/src/SuiteSparse/UFconfig
    )

  Set(HEADERS ${HEADERS}
    Tifpack_AMDReordering.hpp
    Tifpack_Amesos.hpp
    )
  
  SET(SOURCES ${SOURCES}
    Tifpack_AMDReordering.cpp
    Tifpack_Amesos.cpp
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
