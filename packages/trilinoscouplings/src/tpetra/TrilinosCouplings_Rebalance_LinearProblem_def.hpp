// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TRILINOSCOUPLINGS_REBALANCE_LINEARPROBLEM_DEF_HPP
#define TRILINOSCOUPLINGS_REBALANCE_LINEARPROBLEM_DEF_HPP

/// \file TrilinosCouplings_Rebalance_LinearProblem_def.hpp
/// \brief Definition of the TrilinosCouplings::Rebalance_LinearProblem class

#include <TrilinosCouplings_Rebalance_LinearProblem_decl.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>

namespace TrilinosCouplings {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Rebalance_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
Rebalance_LinearProblem( Teuchos::RCP< Teuchos::ParameterList > paramListForZoltan2PartitioningProblem )
  : Tpetra::ViewTransform< Tpetra::LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >()
  , paramListForZoltan2PartitioningProblem_( paramListForZoltan2PartitioningProblem )
{
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Rebalance_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~Rebalance_LinearProblem()
{
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename Rebalance_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
Rebalance_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::operator()( OriginalType const & origProblem )
{
  using cg_t  = Tpetra::CrsGraph     <        LocalOrdinal, GlobalOrdinal, Node>;
  using mv_t  = Tpetra::MultiVector  <Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using cm_t  = Tpetra::CrsMatrix    <Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using lp_t  = Tpetra::LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  // Save original object
  this->origObj_ = origProblem;

  Teuchos::RCP<cm_t>       origMatrix = Teuchos::rcp<cm_t>( dynamic_cast<cm_t *>(origProblem->getMatrix().get()), false );
  Teuchos::RCP<cg_t const> origGraph  = Teuchos::rcp<cg_t const>( dynamic_cast<cg_t const *>(origMatrix->getGraph().get()), false );
  Teuchos::RCP<mv_t>       origLHS    = origProblem->getLHS();
  Teuchos::RCP<mv_t>       origRHS    = origProblem->getRHS();

  // ****************************************************************
  // Instantiate a partitioning problem and solve it.
  // ****************************************************************
  using GraphAdapter_t = Zoltan2::XpetraCrsGraphAdapter<cg_t>;
  GraphAdapter_t graphAdapter(origGraph);
  Zoltan2::PartitioningProblem<GraphAdapter_t> partitioningProblem(&graphAdapter, paramListForZoltan2PartitioningProblem_.get());
  partitioningProblem.solve();

  // ****************************************************************
  // Rebalance the graph
  // ****************************************************************
  Teuchos::RCP<cg_t> newGraph( Teuchos::null );
  graphAdapter.applyPartitioningSolution(*origGraph, newGraph, partitioningProblem.getSolution());

  // ****************************************************************
  // Rebalance the matrix
  // ****************************************************************
  using MatrixAdapter_t = Zoltan2::XpetraCrsMatrixAdapter<cm_t>;
  MatrixAdapter_t matrixAdapter(origMatrix);

  Teuchos::RCP<cm_t> newMatrix( Teuchos::null );
  matrixAdapter.applyPartitioningSolution(*origMatrix, newMatrix, partitioningProblem.getSolution());
  
  // ****************************************************************
  // Rebalance the lhs vector
  // ****************************************************************
  using MultiVectorAdapter_t = Zoltan2::XpetraMultiVectorAdapter<mv_t>;
  MultiVectorAdapter_t lhsAdapter(origLHS);

  Teuchos::RCP<mv_t> newLHS( Teuchos::null );
  lhsAdapter.applyPartitioningSolution(*origLHS, newLHS, partitioningProblem.getSolution());

  // ****************************************************************
  // Rebalance the rhs vector
  // ****************************************************************
  MultiVectorAdapter_t rhsAdapter(origRHS);

  Teuchos::RCP<mv_t> newRHS( Teuchos::null );
  rhsAdapter.applyPartitioningSolution(*origRHS, newRHS, partitioningProblem.getSolution());

  this->newObj_ = Teuchos::rcp<lp_t>( new lp_t(newMatrix, newLHS, newRHS) );

  return this->newObj_;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the TrilinosCouplings namespace!
//

#define TRILINOSCOUPLINGS_REBALANCE_LINEARPROBLEM_INSTANT(SCALAR,LO,GO,NODE) \
  template class Rebalance_LinearProblem< SCALAR , LO , GO , NODE >;

} // namespace TrilinosCouplings

#endif // TRILINOSCOUPLINGS_REBALANCE_LINEARPROBLEM_DEF_HPP
