// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_TpetraReorderedMappingStrategy.hpp"

#include "Teko_BlockedReordering.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;

namespace Teko {
namespace TpetraHelpers {

TpetraReorderedMappingStrategy::TpetraReorderedMappingStrategy(
    const BlockReorderManager& brm, const Teuchos::RCP<const MappingStrategy>& map)
    : reorderManager_(brm), mapStrategy_(map) {
  rangeMap_  = mapStrategy_->rangeMap();
  domainMap_ = mapStrategy_->domainMap();
}

void TpetraReorderedMappingStrategy::copyTpetraIntoThyra(
    const Tpetra::MultiVector<ST, LO, GO, NT>& X,
    const Teuchos::Ptr<Thyra::MultiVectorBase<ST> >& thyra_X) const {
  using Teuchos::ptr_const_cast;
  using Teuchos::rcp_const_cast;

  // first flatten the vector: notice this just works on the block structure
  RCP<Thyra::ProductMultiVectorBase<ST> > prod_X =
      rcp_dynamic_cast<Thyra::ProductMultiVectorBase<ST> >(rcpFromRef(*thyra_X));
  RCP<Thyra::MultiVectorBase<ST> > flat_X = buildFlatMultiVector(reorderManager_, prod_X);

  // now use the underlying mapping strategy to copy the flat vector
  mapStrategy_->copyTpetraIntoThyra(X, flat_X.ptr());
}

void TpetraReorderedMappingStrategy::copyThyraIntoTpetra(
    const RCP<const Thyra::MultiVectorBase<ST> >& thyra_Y,
    Tpetra::MultiVector<ST, LO, GO, NT>& Y) const {
  // first flatten the vector: notice this just works on the block structure
  RCP<const Thyra::ProductMultiVectorBase<ST> > prod_Y =
      rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<ST> >(rcpFromRef(*thyra_Y));
  RCP<const Thyra::MultiVectorBase<ST> > flat_Y = buildFlatMultiVector(reorderManager_, prod_Y);

  // now use the underlying mapping strategy to copy the flat vector
  mapStrategy_->copyThyraIntoTpetra(flat_Y, Y);
}

}  // end namespace TpetraHelpers
}  // end namespace Teko
