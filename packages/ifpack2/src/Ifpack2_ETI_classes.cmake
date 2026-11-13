# Set the list of Ifpack2 classes for which we want to do ETI.
#
# This list lives in its own file (rather than inline in
# src/CMakeLists.txt) so that downstream packages such as Stokhos can
# INCLUDE it directly.  That is required when Ifpack2 is not configured
# in the same CMake project as the downstream package (e.g. when
# Ifpack2 is pre-installed and consumed as an external package), in
# which case Ifpack2's src/CMakeLists.txt is never processed.
#
# IFPACK2_ETI_CLASSES and IFPACK2_ETI_LO_GO_CLASSES are global
# variables so that Stokhos can do ETI.
#
# Note: refer to Ifpack2_* option variables explicitly (not
# ${PACKAGE_NAME}_*) since this file may be included from another
# package's CMakeLists.txt.

# Set the list of Ifpack2 classes templated on <Scalar, LO, GO,
# Node> for which we want to do ETI using this system.  These
# classes usually operate on sparse matrices (instances of
# Tpetra::CrsMatrix or Tpetra::RowMatrix, which also take these
# template parameters).
GLOBAL_SET(
  IFPACK2_ETI_CLASSES
  AdditiveSchwarz
  BlockRelaxation
  BorderedOperator
  Chebyshev
  Container
  BandedContainer
  DatabaseSchwarz
  DenseContainer
  Diagonal
  DiagonalFilter
  DropFilter
  Factory
  Hiptmair
  IdentitySolver
  ILUT
  LinePartitioner
  LocalFilter
  LocalSparseTriangularSolver
  MDF
  OverlappingRowMatrix
  Relaxation
  ReorderFilter
  RILUK
  SingletonFilter
  SparseContainer
  SparsityFilter
  ContainerFactory
  TriDiContainer
  Details::AdditiveSchwarzFilter
  Details::Chebyshev
  Details::ChebyshevKernel
  Details::DenseSolver
  Details::Factory
  Details::InverseDiagonalKernel
  Details::LinearSolver
  Details::OneLevelFactory
  Details::ScaledDampedResidual
  Details::TriDiSolver
  Experimental::RBILUK
  Details::LinearSolverFactory
)

IF(Ifpack2_ENABLE_Amesos2)
  APPEND_GLOBAL_SET(IFPACK2_ETI_CLASSES Details::Amesos2Wrapper)
ENDIF()

IF(Ifpack2_ENABLE_Experimental_KokkosKernels_Features)
  APPEND_GLOBAL_SET(IFPACK2_ETI_CLASSES BlockTriDiContainer)
  APPEND_GLOBAL_SET(IFPACK2_ETI_CLASSES BlockComputeResidualVector)
  APPEND_GLOBAL_SET(IFPACK2_ETI_CLASSES BlockComputeResidualAndSolve)
ENDIF()

IF(Ifpack2_ENABLE_ShyLU_NodeFastILU)
  APPEND_GLOBAL_SET(
    IFPACK2_ETI_CLASSES
    Details::FastILU_Base
    Details::Filu
    Details::Fic
    Details::Fildl
  )
ENDIF()

IF(Ifpack2_ENABLE_HYPRE)
  APPEND_GLOBAL_SET(IFPACK2_ETI_CLASSES Hypre)
ENDIF()

# Set the list of Ifpack2 classes templated on <LO, GO, Node> for
# which we want to do ETI using this system.  These classes usually
# operate on sparse graphs (instances of Tpetra::CrsGraph or
# Tpetra::RowGraph, which also take these template parameters).
GLOBAL_SET(
  IFPACK2_ETI_LO_GO_CLASSES
  LinearPartitioner
  OverlappingPartitioner
  Details::OverlappingRowGraph
)
