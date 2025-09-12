// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_REBALANCE_LINEARPROBLEM_DEF_HPP
#define TPETRA_REBALANCE_LINEARPROBLEM_DEF_HPP

/// \file Tpetra_Rebalance_LinearProblem_def.hpp
/// \brief Definition of the Tpetra::Rebalance_LinearProblem class
///
/// If you want to use Tpetra::Rebalance_LinearProblem, include
/// "Tpetra_Rebalance_LinearProblem.hpp", a file which CMake generates
/// and installs for you.
///
/// If you only want the declaration of Tpetra::Rebalance_LinearProblem,
/// include "Tpetra_Rebalance_LinearProblem_decl.hpp".

#include <Tpetra_Rebalance_LinearProblem_decl.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Rebalance_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
Rebalance_LinearProblem( Teuchos::RCP< Teuchos::ParameterList > paramListForZoltan2PartitioningProblem )
  : ViewTransform< LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >()
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
  using mv_t  = MultiVector  <Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using cm_t  = CrsMatrix    <Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using lp_t  = LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  // Save original object
  this->origObj_ = origProblem;

  Teuchos::RCP<cm_t> origMatrix = Teuchos::rcp<cm_t>( dynamic_cast<cm_t *>(origProblem->getMatrix().get()), false );
  Teuchos::RCP<mv_t> origLHS    = origProblem->getLHS();
  Teuchos::RCP<mv_t> origRHS    = origProblem->getRHS();

  // ****************************************************************
  // Instantiate a partitioning problem and solve it.
  // ****************************************************************
  using MatrixAdapter_t = Zoltan2::XpetraCrsMatrixAdapter<cm_t>;
  MatrixAdapter_t matrixAdapter(origMatrix);

  Zoltan2::PartitioningProblem<MatrixAdapter_t> partitioningProblem(&matrixAdapter, paramListForZoltan2PartitioningProblem_.get());
  partitioningProblem.solve();

  // ****************************************************************
  // Rebalance the matrix
  // ****************************************************************
  Teuchos::RCP<cm_t> newMatrix( Teuchos::null );
  matrixAdapter.applyPartitioningSolution(*origMatrix,
                                          newMatrix, // The rebalanced one
                                          partitioningProblem.getSolution());

  // ****************************************************************
  // Rebalance the lhs vector
  // ****************************************************************
  using MultiVectorAdapter_t = Zoltan2::XpetraMultiVectorAdapter<mv_t>;
  MultiVectorAdapter_t lhsAdapter(origLHS);

  Teuchos::RCP<mv_t> newLHS( Teuchos::null );
  lhsAdapter.applyPartitioningSolution(*origLHS,
                                       newLHS, // The rebalanced one
                                       partitioningProblem.getSolution());

  // ****************************************************************
  // Rebalance the rhs vector
  // ****************************************************************
  MultiVectorAdapter_t rhsAdapter(origRHS);

  Teuchos::RCP<mv_t> newRHS( Teuchos::null );
  rhsAdapter.applyPartitioningSolution(*origRHS,
                                       newRHS, // The rebalanced one
                                       partitioningProblem.getSolution());

  this->newObj_ = Teuchos::rcp<lp_t>( new lp_t(newMatrix, newLHS, newRHS) );

  return this->newObj_;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_REBALANCELINEARPROBLEM_INSTANT(SCALAR,LO,GO,NODE) \
  template class Rebalance_LinearProblem< SCALAR , LO , GO , NODE >;

} // namespace Tpetra

#endif // TPETRA_REBALANCE_LINEARPROBLEM_DEF_HPP

