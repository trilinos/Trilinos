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

#include "Teuchos_Assert.hpp"
#include <TrilinosCouplings_Rebalance_LinearProblem_decl.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningHelpers.hpp>

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
  using map_t = Tpetra::Map          <        LocalOrdinal, GlobalOrdinal, Node>;

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

  // Build importer
  {
    Teuchos::ArrayRCP<GlobalOrdinal> importList;
    Zoltan2::getImportList<GraphAdapter_t, GraphAdapter_t>(partitioningProblem.getSolution(), &graphAdapter, importList);
    // source map
    auto smap = origGraph->getRowMap();
    // target map
    auto tmap = rcp(new map_t(smap->getGlobalNumElements(), importList(), smap->getIndexBase(), smap->getComm()));
    // importer
    importer_ = Teuchos::rcp(new import_t(smap, tmap));
  }

  // ****************************************************************
  // Rebalance the matrix
  // ****************************************************************
  Teuchos::RCP<cm_t> newMatrix;
  origMatrix->importAndFillComplete(newMatrix, *importer_, importer_->getTargetMap(), importer_->getTargetMap());

  // ****************************************************************
  // Rebalance the lhs vector
  // ****************************************************************
  Teuchos::RCP<mv_t> newLHS = Teuchos::rcp(new mv_t(importer_->getTargetMap(), origLHS->getNumVectors(), false));
  newLHS->doImport(*origLHS, *importer_, Tpetra::INSERT);

  // ****************************************************************
  // Rebalance the rhs vector
  // ****************************************************************
  Teuchos::RCP<mv_t> newRHS = Teuchos::rcp(new mv_t(importer_->getTargetMap(), origRHS->getNumVectors(), false));
  newRHS->doImport(*origRHS, *importer_, Tpetra::INSERT);

  this->newObj_ = Teuchos::rcp<lp_t>( new lp_t(newMatrix, newLHS, newRHS) );

  return this->newObj_;
}


template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
Rebalance_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::fwd() {
  auto origProblem = this->origObj_;
  auto origRHS = origProblem->getRHS();
  TEUCHOS_ASSERT(this->isConstructed());
  auto rebalProblem = this->newObj_;
  auto rebalRHS = rebalProblem->getRHS();
  rebalRHS->doImport(*origRHS, *importer_, Tpetra::INSERT);
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
Rebalance_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::rvs() {
  auto origProblem = this->origObj_;
  auto origLHS = origProblem->getLHS();
  TEUCHOS_ASSERT(this->isConstructed());
  auto rebalProblem = this->newObj_;
  auto rebalLHS = rebalProblem->getLHS();
  origLHS->doExport(*rebalLHS, *importer_, Tpetra::INSERT);
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
