// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_CREATEOVERLAPGRAPH_HPP
#define IFPACK2_CREATEOVERLAPGRAPH_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Import.hpp"
#include "Teuchos_RefCountPtr.hpp"


namespace Ifpack2 {

/// \brief Construct an overlapped graph for use with Ifpack2 preconditioners.
/// \tparam GraphType A specialization of Tpetra::RowGraph or Tpetra::CrsGraph.
///
/// \warning This function is an implementation detail of Ifpack2.
///   Its interface may change or it may go away at any time.
///
/// \note This method has only been tested with GraphType =
///   Tpetra::CrsGraph.  It should also work with GraphType =
///   Tpetra::RowGraph, but I have not tested this case.
///
/// \param inputGraph [in] The input graph.  We assume that its row
///   Map is nonoverlapping.
///
/// \param overlapLevel [in] The level of overlap.  Zero means no
///   overlap, in which case this function just returns the original
///   \c inputGraph.
template<class GraphType>
Teuchos::RCP<const GraphType>
createOverlapGraph (const Teuchos::RCP<const GraphType>& inputGraph,
                    const int overlapLevel)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Tpetra::Map<typename GraphType::local_ordinal_type,
                      typename GraphType::global_ordinal_type,
                      typename GraphType::node_type> map_type;
  typedef Tpetra::Import<typename GraphType::local_ordinal_type,
                         typename GraphType::global_ordinal_type,
                         typename GraphType::node_type> import_type;
  TEUCHOS_TEST_FOR_EXCEPTION(
    overlapLevel < 0, std::invalid_argument,
    "Ifpack2::createOverlapGraph: overlapLevel must be >= 0, "
    "but you specified overlapLevel = " << overlapLevel << ".");

  const int numProcs = inputGraph->getMap ()->getComm ()->getSize ();
  if (overlapLevel == 0 || numProcs < 2) {
    return inputGraph;
  }

  RCP<const map_type> overlapRowMap = inputGraph->getRowMap ();
  RCP<const map_type> domainMap = inputGraph->getDomainMap ();
  RCP<const map_type> rangeMap = inputGraph->getRangeMap ();

  RCP<GraphType> overlapGraph;
  RCP<const GraphType> oldGraph;
  RCP<const map_type> oldRowMap;
  for (int level = 0; level < overlapLevel; ++level) {
    oldGraph = overlapGraph;
    oldRowMap = overlapRowMap;

    RCP<const import_type> overlapImporter;
    if (level == 0) {
      overlapImporter = inputGraph->getImporter ();
    } else {
      overlapImporter = oldGraph->getImporter ();
    }

    overlapRowMap = overlapImporter->getTargetMap ();
    if (level < overlapLevel - 1) {
      overlapGraph = rcp (new GraphType (overlapRowMap, 0));
    }
    else {
      // On last iteration, we want to filter out all columns except those that
      // correspond to rows in the graph.  This ensures that our graph is square
      overlapGraph = rcp (new GraphType (overlapRowMap, overlapRowMap, 0));
    }

    overlapGraph->doImport (*inputGraph, *overlapImporter, Tpetra::INSERT);
    overlapGraph->fillComplete (domainMap, rangeMap);
  }

  return overlapGraph;
}

/// \brief Construct an overlapped matrix for use with Ifpack2 preconditioners.
/// \tparam MatrixType A specialization of Tpetra::CrsMatrix.
///
/// \param inputMatrix [in] The input matrix.  We assume that its row
///   Map is nonoverlapping.
///
/// \param overlapLevel [in] The level of overlap.  Zero means no
///   overlap, in which case this function just returns the original
///   \c inputMatrix.
template<class MatrixType>
Teuchos::RCP<const MatrixType>
createOverlapMatrix (const Teuchos::RCP<const MatrixType>& inputMatrix,
                     const int overlapLevel)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef typename MatrixType::map_type map_type;
  typedef Tpetra::Import<typename MatrixType::local_ordinal_type,
    typename MatrixType::global_ordinal_type,
    typename MatrixType::node_type> import_type;

  TEUCHOS_TEST_FOR_EXCEPTION(
    overlapLevel < 0, std::invalid_argument,
    "Ifpack2::createOverlapMatrix: overlapLevel must be >= 0, "
    "but you specified overlapLevel = " << overlapLevel << ".");

  const int numProcs = inputMatrix->getMap ()->getComm ()->getSize ();
  if (overlapLevel == 0 || numProcs < 2) {
    return inputMatrix;
  }

  RCP<const map_type> overlapRowMap = inputMatrix->getRowMap ();
  RCP<const map_type> domainMap = inputMatrix->getDomainMap ();
  RCP<const map_type> rangeMap = inputMatrix->getRangeMap ();

  RCP<MatrixType> overlapMatrix;
  RCP<const MatrixType> oldMatrix;
  RCP<const map_type> oldRowMap;
  for (int level = 0; level < overlapLevel; ++level) {
    oldMatrix = overlapMatrix;
    oldRowMap = overlapRowMap;

    RCP<const import_type> overlapImporter;
    if (level == 0) {
      overlapImporter = inputMatrix->getGraph ()->getImporter ();
    } else {
      overlapImporter = oldMatrix->getGraph ()->getImporter ();
    }

    overlapRowMap = overlapImporter->getTargetMap ();
    if (level < overlapLevel - 1) {
      overlapMatrix = rcp (new MatrixType (overlapRowMap, 0));
    }
    else {
      // On last iteration, we want to filter out all columns except those that
      // correspond to rows in the matrix.  This assures that our matrix is square
      overlapMatrix = rcp (new MatrixType (overlapRowMap, overlapRowMap, 0));
    }

    overlapMatrix->doImport (*inputMatrix, *overlapImporter, Tpetra::INSERT);
    overlapMatrix->fillComplete (domainMap, rangeMap);
  }

  return overlapMatrix;
}

} // namespace Ifpack2

#endif // IFPACK2_CREATEOVERLAPGRAPH_HPP
