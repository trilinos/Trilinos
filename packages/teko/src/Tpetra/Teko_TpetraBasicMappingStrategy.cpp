// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_TpetraBasicMappingStrategy.hpp"
#include "Teko_TpetraHelpers.hpp"

#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_DefaultSpmdMultiVector.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

namespace Teko {
namespace TpetraHelpers {

// Creates a strided mapping strategy. This class is useful
// for breaking up nodally ordered matrices (i.e. the unknowns
// in a FEM problem are ordered [u0,v0,p0,u1,v1,p1,...]). Current
// implimentation only supports a fixed number of variables
//
//    arguments:
//       vars - Number of different variables
//       map  - original Epetra_Map to be broken up
//       comm - Epetra_Comm object related to the map
//
BasicMappingStrategy::BasicMappingStrategy(const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& rMap,
                                           const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& dMap,
                                           const Teuchos::Comm<Thyra::Ordinal>& /* comm */) {
  rangeMap_  = rMap;
  domainMap_ = dMap;
}

// Virtual function defined in MappingStrategy.  This copies
// an Epetra_MultiVector into a Thyra::MultiVectorBase with
// blocking handled by the strides defined in the constructor.
//
//   arguments:
//      X       - source Epetra_MultiVector
//      thyra_X - destination Thyra::MultiVectorBase
//
void BasicMappingStrategy::copyTpetraIntoThyra(
    const Tpetra::MultiVector<ST, LO, GO, NT>& X,
    const Teuchos::Ptr<Thyra::MultiVectorBase<ST> >& thyra_X) const {
  // perform a simple copy
  // RCP<Thyra::DefaultSpmdMultiVector<ST> > vec
  //         = rcp_dynamic_cast<Thyra::DefaultSpmdMultiVector<ST> >(Teuchos::rcpFromRef(*thyra_X));
  RCP<Thyra::TpetraMultiVector<ST, LO, GO, NT> > vec =
      rcp_dynamic_cast<Thyra::TpetraMultiVector<ST, LO, GO, NT> >(Teuchos::rcpFromRef(*thyra_X));
  Teuchos::RCP<Tpetra::MultiVector<ST, LO, GO, NT> > ptrX =
      Teuchos::rcp_const_cast<Tpetra::MultiVector<ST, LO, GO, NT> >(Teuchos::rcpFromRef(X));
  fillDefaultSpmdMultiVector(vec, ptrX);
}

// Virtual function defined in MappingStrategy.  This copies
// an Epetra_MultiVector into a Thyra::MultiVectorBase with
// blocking handled by the strides defined in the constructor.
//
//   arguments:
//      thyra_Y - source Thyra::MultiVectorBase
//      Y       - destination Epetra_MultiVector
//
void BasicMappingStrategy::copyThyraIntoTpetra(
    const RCP<const Thyra::MultiVectorBase<ST> >& thyra_Y,
    Tpetra::MultiVector<ST, LO, GO, NT>& Y) const {
  RCP<const Tpetra::MultiVector<ST, LO, GO, NT> > tSrc =
      Thyra::TpetraOperatorVectorExtraction<ST, LO, GO, NT>::getConstTpetraMultiVector(thyra_Y);

  Y = *tSrc;
}

}  // namespace TpetraHelpers
}  // end namespace Teko
