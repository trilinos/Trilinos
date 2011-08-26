
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE( Ifpack2 )

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Trilinos
  HAVE_${PACKAGE_NAME}_Trilinos
  "Building as a Trilinos package"
  ON
  )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_DEBUG
  HAVE_${PACKAGE_NAME}_DEBUG
  "Enable runtime debug checking."
  ${Trilinos_ENABLE_DEBUG}
  )

PACKAGE_ADD_EXPLICIT_INSTANTIATION_OPTION()

#Adding a define for when the TPL QD is enabled.
#This replaces an option that was reusing an existing name.
SET(HAVE_${PACKAGE_NAME}_QD "${PACKAGE_NAME}_ENABLE_QD")

#Adding a define for when MPI is enabled.
#This replaces an option that was reusing an existing name.
SET(HAVE_${PACKAGE_NAME}_MPI "${${PACKAGE_NAME}_ENABLE_MPI}")

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)

PACKAGE_ARCH_EXCLUDE_FILES(
  src/Tifpack_AMDReordering.cpp
  src/Tifpack_AMDReordering.hpp
  src/Tifpack_AdditiveSchwarz.hpp
  src/Tifpack_Amesos.cpp
  src/Tifpack_Amesos.hpp
  src/Tifpack_BlockRelaxation.hpp
  src/Tifpack_ConstructLevelFillGraph.hpp
  src/Tifpack_Container.hpp
  src/Tifpack_CrsGraph.hpp
  src/Tifpack_DenseContainer.cpp
  src/Tifpack_DenseContainer.hpp
  src/Tifpack_DiagPreconditioner.cpp
  src/Tifpack_DiagPreconditioner.hpp
  src/Tifpack_DiagonalFilter.cpp
  src/Tifpack_DiagonalFilter.hpp
  src/Tifpack_DropFilter.cpp
  src/Tifpack_DropFilter.hpp
  src/Tifpack_EquationPartitioner.cpp
  src/Tifpack_EquationPartitioner.hpp
  src/Tifpack_Graph.hpp
  src/Tifpack_Graph_Tpetra_CrsGraph.cpp
  src/Tifpack_Graph_Tpetra_CrsGraph.hpp
  src/Tifpack_Graph_Tpetra_RowMatrix.cpp
  src/Tifpack_Graph_Tpetra_RowMatrix.hpp
  src/Tifpack_GreedyPartitioner.cpp
  src/Tifpack_GreedyPartitioner.hpp
  src/Tifpack_HashTable.cpp
  src/Tifpack_HashTable.hpp
  src/Tifpack_IC.cpp
  src/Tifpack_IC.hpp
  src/Tifpack_ICT.cpp
  src/Tifpack_ICT.hpp
  src/Tifpack_IC_Utils.cpp
  src/Tifpack_IC_Utils.hpp
  src/Tifpack_IKLU.cpp
  src/Tifpack_IKLU.hpp
  src/Tifpack_IKLU_Utils.cpp
  src/Tifpack_IKLU_Utils.hpp
  src/Tifpack_ILU.cpp
  src/Tifpack_ILU.hpp
  src/Tifpack_LinearPartitioner.cpp
  src/Tifpack_LinearPartitioner.hpp
  src/Tifpack_LocalFilter.cpp
  src/Tifpack_LocalFilter.hpp
  src/Tifpack_METISPartitioner.cpp
  src/Tifpack_METISPartitioner.hpp
  src/Tifpack_METISReordering.cpp
  src/Tifpack_METISReordering.hpp
  src/Tifpack_NodeFilter.cpp
  src/Tifpack_NodeFilter.hpp
  src/Tifpack_OverlapFactor.cpp
  src/Tifpack_OverlapFactorObject.hpp
  src/Tifpack_OverlapSolveObject.cpp
  src/Tifpack_OverlapSolveObject.hpp
  src/Tifpack_OverlappingPartitioner.cpp
  src/Tifpack_OverlappingPartitioner.hpp
  src/Tifpack_OverlappingRowMatrix.cpp
  src/Tifpack_OverlappingRowMatrix.hpp
  src/Tifpack_Partitioner.hpp
  src/Tifpack_PerturbedMatrix.hpp
  src/Tifpack_RCMReordering.cpp
  src/Tifpack_RCMReordering.hpp
  src/Tifpack_ReorderFilter.cpp
  src/Tifpack_ReorderFilter.hpp
  src/Tifpack_Reordering.hpp
  src/Tifpack_SPARSKIT.cpp
  src/Tifpack_SPARSKIT.hpp
  src/Tifpack_SingletonFilter.cpp
  src/Tifpack_SingletonFilter.hpp
  src/Tifpack_SparseContainer.hpp
  src/Tifpack_SparsityFilter.cpp
  src/Tifpack_SparsityFilter.hpp
  src/Tifpack_UserPartitioner.cpp
  src/Tifpack_UserPartitioner.hpp
  src/Tifpack_Utils.cpp
  src/Tifpack_Utils.hpp
  )

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()


