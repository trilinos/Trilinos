// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_TpetraInverseFactoryOperator.hpp"

// Teko includes
#include "Teko_TpetraBasicMappingStrategy.hpp"

// Thyra includes
#include "Thyra_TpetraLinearOp.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;

namespace Teko {
namespace TpetraHelpers {

/** \brief Constructor that takes the InverseFactory that will
 *        build the operator.
 *
 * Constructor that takes the InverseFactory that will
 * build the operator.
 */
InverseFactoryOperator::InverseFactoryOperator(const Teuchos::RCP<const InverseFactory>& ifp)
    : inverseFactory_(ifp), firstBuildComplete_(false), setConstFwdOp_(true) {}

/** \brief Build the underlying data structure for the inverse operator.
 *
 * Build the underlying data structure for the inverse operator. This
 * permits the manipulation of the state object for an inverse operator.
 *
 * \param[in] clearOld If true any previously constructed
 *                     operator will be wiped out and
 *                     a new one created. If false, anoperator
 *                     will be created only if the current one is
 *                     empty (i.e. <code>initPreconditioner</code>
 *                     had not been called).
 */
void InverseFactoryOperator::initInverse(bool clearOld) {
  if (not clearOld) return;
  invOperator_ = Teuchos::null;
}

/** \brief Build this inverse operator from an Tpetra::Operator<ST,LO,GO,NT>
 * passed in to this object.
 *
 * Build this inverse opeerator from an Tpetra::Operator<ST,LO,GO,NT>
 * passed in to this object. If this Tpetra::Operator<ST,LO,GO,NT>
 * is an EpetraOperatorWrapper object then the block Thyra components
 * are extracted.
 *
 * \param[in] A The Epetra source operator.
 * \param[in] clear If true, than any previous state saved by the operator
 *                  is discarded.
 */
void InverseFactoryOperator::buildInverseOperator(
    const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> >& A, bool clear) {
  Teko_DEBUG_SCOPE("InverseFactoryOperator::buildInverseOperator", 10);

  // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
  RCP<const Thyra::LinearOpBase<ST> > thyraA = extractLinearOp(A);

  // set the mapping strategy
  SetMapStrategy(rcp(new InverseMappingStrategy(extractMappingStrategy(A))));

  initInverse(clear);

  // actually build the inverse operator
  invOperator_ = Teko::buildInverse(*inverseFactory_, thyraA);

  SetOperator(invOperator_, false);

  firstBuildComplete_ = true;

  if (setConstFwdOp_) fwdOp_ = A;

  setConstFwdOp_ = true;

  TEUCHOS_ASSERT(invOperator_ != Teuchos::null);
  TEUCHOS_ASSERT(getForwardOp() != Teuchos::null);
  TEUCHOS_ASSERT(getThyraOp() != Teuchos::null);
  TEUCHOS_ASSERT(getMapStrategy() != Teuchos::null);
  TEUCHOS_ASSERT(firstBuildComplete_ == true);
}

void InverseFactoryOperator::buildInverseOperator(
    const Teuchos::RCP<Tpetra::Operator<ST, LO, GO, NT> >& A, bool /* clear */) {
  setConstFwdOp_ = false;

  fwdOp_.initialize(A);

  buildInverseOperator(A.getConst());

  TEUCHOS_ASSERT(setConstFwdOp_ == true);
}

/** \brief Rebuild this inverse from an Tpetra::Operator<ST,LO,GO,NT>  passed
 * in this to object.
 *
 * Rebuild this inverse from an Tpetra::Operator<ST,LO,GO,NT>  passed
 * in this to object.  If <code>buildInverseOperator</code> has not been called
 * the inverse operator will be built instead. Otherwise efforts are taken
 * to only rebuild what is neccessary. Also, that this Tpetra::Operator<ST,LO,GO,NT>
 * may be an EpetraOperatorWrapper object, so the block Thyra components
 * can be extracted.
 *
 * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
 */
void InverseFactoryOperator::rebuildInverseOperator(
    const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> >& A) {
  Teko_DEBUG_SCOPE("InverseFactoryOperator::rebuildPreconditioner", 10);

  // if the inverse hasn't been built yet, rebuild from scratch
  if (not firstBuildComplete_) {
    buildInverseOperator(A, false);
    return;
  }

  RCP<const Thyra::LinearOpBase<ST> > thyraA = extractLinearOp(A);
  Teko::rebuildInverse(*inverseFactory_, thyraA, invOperator_);

  if (setConstFwdOp_) fwdOp_.initialize(A);

  SetOperator(invOperator_, false);

  setConstFwdOp_ = true;

  TEUCHOS_ASSERT(getForwardOp() != Teuchos::null);
  TEUCHOS_ASSERT(invOperator_ != Teuchos::null);
  TEUCHOS_ASSERT(getThyraOp() != Teuchos::null);
  TEUCHOS_ASSERT(firstBuildComplete_ == true);
}

void InverseFactoryOperator::rebuildInverseOperator(
    const Teuchos::RCP<Tpetra::Operator<ST, LO, GO, NT> >& A) {
  setConstFwdOp_ = false;

  fwdOp_.initialize(A);

  // build from constant epetra operator
  rebuildInverseOperator(A.getConst());

  TEUCHOS_ASSERT(setConstFwdOp_ == true);
}

Teuchos::RCP<const Thyra::LinearOpBase<ST> > InverseFactoryOperator::extractLinearOp(
    const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> >& A) const {
  // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
  const RCP<const TpetraOperatorWrapper>& eow = rcp_dynamic_cast<const TpetraOperatorWrapper>(A);

  // if it is an EpetraOperatorWrapper, then get the Thyra operator
  if (eow != Teuchos::null) return eow->getThyraOp();

  // otherwise wrap it up as a thyra operator
  return Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A->getDomainMap()), A);
}

Teuchos::RCP<const MappingStrategy> InverseFactoryOperator::extractMappingStrategy(
    const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> >& A) const {
  // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
  const RCP<const TpetraOperatorWrapper>& eow = rcp_dynamic_cast<const TpetraOperatorWrapper>(A);

  // if it is an EpetraOperatorWrapper, then get the Thyra operator
  if (eow != Teuchos::null) return eow->getMapStrategy();

  // otherwise wrap it up as a thyra operator
  RCP<const Tpetra::Map<LO, GO, NT> > range  = A->getRangeMap();
  RCP<const Tpetra::Map<LO, GO, NT> > domain = A->getDomainMap();
  return rcp(
      new BasicMappingStrategy(range, domain, *Thyra::convertTpetraToThyraComm(range->getComm())));
}

}  // end namespace TpetraHelpers
}  // end namespace Teko
