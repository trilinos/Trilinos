// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_TpetraBlockPreconditioner.hpp"
#include "Teko_Preconditioner.hpp"

// Thyra includes
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_TpetraLinearOp.hpp"

// Teuchos includes
#include "Teuchos_Time.hpp"

// Teko includes
#include "Teko_TpetraBasicMappingStrategy.hpp"
#include "Teko_Utilities.hpp"

namespace Teko {
namespace TpetraHelpers {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;

/** \brief Constructor that takes the BlockPreconditionerFactory that will
 *        build the preconditioner.
 *
 * Constructor that takes the BlockPreconditionerFactory that will
 * build the preconditioner.
 */
TpetraBlockPreconditioner::TpetraBlockPreconditioner(
    const Teuchos::RCP<const PreconditionerFactory>& bfp)
    : preconFactory_(bfp), firstBuildComplete_(false) {}

void TpetraBlockPreconditioner::initPreconditioner(bool clearOld) {
  if ((not clearOld) && preconObj_ != Teuchos::null) return;
  preconObj_ = preconFactory_->createPrec();
}

/** \brief Build this preconditioner from an Epetra_Operator
 * passed in to this object. It is assume that this Epetra_Operator
 *
 * Build this preconditioner from an Epetra_Operator
 * passed in to this object. It is assume that this Epetra_Operator
 * will be a EpetraOperatorWrapper object, so the block Thyra components
 * can be easily extracted.
 *
 * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
 *
 * \note This will clear any internal state stored by the state object
 */
void TpetraBlockPreconditioner::buildPreconditioner(
    const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> >& A, bool clear) {
  Teko_DEBUG_SCOPE("TBP::buildPreconditioner", 10);

  // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
  RCP<const Thyra::LinearOpBase<ST> > thyraA = extractLinearOp(A);

  // set the mapping strategy
  // SetMapStrategy(rcp(new InverseMappingStrategy(eow->getMapStrategy())));
  SetMapStrategy(rcp(new InverseMappingStrategy(extractMappingStrategy(A))));

  // build preconObj_
  initPreconditioner(clear);

  // actually build the preconditioner
  RCP<const Thyra::LinearOpSourceBase<ST> > lOpSrc = Thyra::defaultLinearOpSource(thyraA);
  preconFactory_->initializePrec(lOpSrc, &*preconObj_, Thyra::SUPPORT_SOLVE_UNSPECIFIED);

  // extract preconditioner operator
  RCP<const Thyra::LinearOpBase<ST> > preconditioner = preconObj_->getUnspecifiedPrecOp();

  SetOperator(preconditioner, false);

  firstBuildComplete_ = true;

  TEUCHOS_ASSERT(preconObj_ != Teuchos::null);
  TEUCHOS_ASSERT(getThyraOp() != Teuchos::null);
  TEUCHOS_ASSERT(getMapStrategy() != Teuchos::null);
  TEUCHOS_ASSERT(firstBuildComplete_ == true);
}

/** \brief Build this preconditioner from an Epetra_Operator
 * passed in to this object. It is assume that this Epetra_Operator
 *
 * Build this preconditioner from an Epetra_Operator
 * passed in to this object. It is assume that this Epetra_Operator
 * will be a EpetraOperatorWrapper object, so the block Thyra components
 * can be easily extracted.
 *
 * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
 * \param[in] src A vector that was used to build the source operator.
 *
 * \note This will clear any internal state stored by the state object
 */
void TpetraBlockPreconditioner::buildPreconditioner(
    const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> >& A,
    const Tpetra::MultiVector<ST, LO, GO, NT>& tpetra_mv, bool clear) {
  Teko_DEBUG_SCOPE("TBP::buildPreconditioner - with solution", 10);

  // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
  RCP<const Thyra::LinearOpBase<ST> > thyraA = extractLinearOp(A);

  // set the mapping strategy
  SetMapStrategy(rcp(new InverseMappingStrategy(extractMappingStrategy(A))));

  TEUCHOS_ASSERT(getMapStrategy() != Teuchos::null);

  // build the thyra version of the source multivector
  RCP<Thyra::MultiVectorBase<ST> > thyra_mv =
      Thyra::createMembers(thyraA->range(), tpetra_mv.getNumVectors());
  getMapStrategy()->copyTpetraIntoThyra(tpetra_mv, thyra_mv.ptr());

  // build preconObj_
  initPreconditioner(clear);

  // actually build the preconditioner
  preconFactory_->initializePrec(Thyra::defaultLinearOpSource(thyraA), thyra_mv, &*preconObj_,
                                 Thyra::SUPPORT_SOLVE_UNSPECIFIED);
  RCP<const Thyra::LinearOpBase<ST> > preconditioner = preconObj_->getUnspecifiedPrecOp();

  SetOperator(preconditioner, false);

  firstBuildComplete_ = true;

  TEUCHOS_ASSERT(preconObj_ != Teuchos::null);
  TEUCHOS_ASSERT(getThyraOp() != Teuchos::null);
  TEUCHOS_ASSERT(getMapStrategy() != Teuchos::null);
  TEUCHOS_ASSERT(firstBuildComplete_ == true);
}

/** \brief Rebuild this preconditioner from an Epetra_Operator passed
 * in this to object.
 *
 * Rebuild this preconditioner from an Epetra_Operator passed
 * in this to object.  If <code>buildPreconditioner</code> has not been called
 * the preconditioner will be built instead. Otherwise efforts are taken
 * to only rebuild what is neccessary. Also, it is assumed that this Epetra_Operator
 * will be an EpetraOperatorWrapper object, so the block Thyra components
 * can be easily extracted.
 *
 * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
 * \param[in] mv A vector that was used to build the source operator.
 */
void TpetraBlockPreconditioner::rebuildPreconditioner(
    const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> >& A) {
  Teko_DEBUG_SCOPE("TBP::rebuildPreconditioner", 10);

  // if the preconditioner hasn't been built yet, rebuild from scratch
  if (not firstBuildComplete_) {
    buildPreconditioner(A, false);
    return;
  }
  Teko_DEBUG_EXPR(Teuchos::Time timer(""));

  // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
  Teko_DEBUG_EXPR(timer.start(true));
  RCP<const Thyra::LinearOpBase<ST> > thyraA = extractLinearOp(A);
  Teko_DEBUG_EXPR(timer.stop());
  Teko_DEBUG_MSG("TBP::rebuild get thyraop time =  " << timer.totalElapsedTime(), 2);

  // reinitialize the preconditioner
  Teko_DEBUG_EXPR(timer.start(true));
  preconFactory_->initializePrec(Thyra::defaultLinearOpSource(thyraA), &*preconObj_,
                                 Thyra::SUPPORT_SOLVE_UNSPECIFIED);
  RCP<const Thyra::LinearOpBase<ST> > preconditioner = preconObj_->getUnspecifiedPrecOp();
  Teko_DEBUG_EXPR(timer.stop());
  Teko_DEBUG_MSG("TBP::rebuild initialize prec time =  " << timer.totalElapsedTime(), 2);

  Teko_DEBUG_EXPR(timer.start(true));
  SetOperator(preconditioner, false);
  Teko_DEBUG_EXPR(timer.stop());
  Teko_DEBUG_MSG("TBP::rebuild set operator time =  " << timer.totalElapsedTime(), 2);

  TEUCHOS_ASSERT(preconObj_ != Teuchos::null);
  TEUCHOS_ASSERT(getThyraOp() != Teuchos::null);
  TEUCHOS_ASSERT(firstBuildComplete_ == true);
}

/** \brief Rebuild this preconditioner from an Epetra_Operator passed
 * in this to object.
 *
 * Rebuild this preconditioner from an Epetra_Operator passed
 * in this to object.  If <code>buildPreconditioner</code> has not been called
 * the preconditioner will be built instead. Otherwise efforts are taken
 * to only rebuild what is neccessary. Also, it is assumed that this Epetra_Operator
 * will be an EpetraOperatorWrapper object, so the block Thyra components
 * can be easily extracted.
 *
 * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
 * \param[in] mv A vector that was used to build the source operator.
 */
void TpetraBlockPreconditioner::rebuildPreconditioner(
    const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> >& A,
    const Tpetra::MultiVector<ST, LO, GO, NT>& tpetra_mv) {
  Teko_DEBUG_SCOPE("TBP::rebuildPreconditioner - with solution", 10);

  // if the preconditioner hasn't been built yet, rebuild from scratch
  if (not firstBuildComplete_) {
    buildPreconditioner(A, tpetra_mv, false);
    return;
  }
  Teko_DEBUG_EXPR(Teuchos::Time timer(""));

  // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
  Teko_DEBUG_EXPR(timer.start(true));
  RCP<const Thyra::LinearOpBase<ST> > thyraA = extractLinearOp(A);
  Teko_DEBUG_EXPR(timer.stop());
  Teko_DEBUG_MSG("TBP::rebuild get thyraop time =  " << timer.totalElapsedTime(), 2);

  // build the thyra version of the source multivector
  Teko_DEBUG_EXPR(timer.start(true));
  RCP<Thyra::MultiVectorBase<ST> > thyra_mv =
      Thyra::createMembers(thyraA->range(), tpetra_mv.getNumVectors());
  getMapStrategy()->copyTpetraIntoThyra(tpetra_mv, thyra_mv.ptr());
  Teko_DEBUG_EXPR(timer.stop());
  Teko_DEBUG_MSG("TBP::rebuild vector copy time =  " << timer.totalElapsedTime(), 2);

  // reinitialize the preconditioner
  Teko_DEBUG_EXPR(timer.start(true));
  preconFactory_->initializePrec(Thyra::defaultLinearOpSource(thyraA), thyra_mv, &*preconObj_,
                                 Thyra::SUPPORT_SOLVE_UNSPECIFIED);
  RCP<const Thyra::LinearOpBase<ST> > preconditioner = preconObj_->getUnspecifiedPrecOp();
  Teko_DEBUG_EXPR(timer.stop());
  Teko_DEBUG_MSG("TBP::rebuild initialize prec time =  " << timer.totalElapsedTime(), 2);

  Teko_DEBUG_EXPR(timer.start(true));
  SetOperator(preconditioner, false);
  Teko_DEBUG_EXPR(timer.stop());
  Teko_DEBUG_MSG("TBP::rebuild set operator time =  " << timer.totalElapsedTime(), 2);

  TEUCHOS_ASSERT(preconObj_ != Teuchos::null);
  TEUCHOS_ASSERT(getThyraOp() != Teuchos::null);
  TEUCHOS_ASSERT(firstBuildComplete_ == true);
}

/** Try to get a <code>Teko::PreconditionerState</code> object. This method
 * attempts to cast its internal representation of a preconditioner
 * object to a <code>Teko::BlockPreconditioner</code> object.  If it suceeds a
 * state object is returned.  Otherwise, <code>Teuchos::null</code> is returned.
 *
 * \returns Get the state object associated with this preconditioner.
 *          If it doesn't exist for this type of preconditioner factory
 *          this method returns null.
 */
Teuchos::RCP<PreconditionerState> TpetraBlockPreconditioner::getPreconditionerState() {
  Teuchos::RCP<Preconditioner> bp = rcp_dynamic_cast<Preconditioner>(preconObj_);

  if (bp != Teuchos::null) return bp->getStateObject();

  return Teuchos::null;
}

/** Try to get a <code>Teko::PreconditionerState</code> object. This method
 * attempts to cast its internal representation of a preconditioner
 * object to a <code>Teko::Preconditioner</code> object.  If it suceeds a
 * state object is returned.  Otherwise, <code>Teuchos::null</code> is returned.
 *
 * \returns Get the state object associated with this preconditioner.
 *          If it doesn't exist for this type of preconditioner factory
 *          this method returns null.
 */
Teuchos::RCP<const PreconditionerState> TpetraBlockPreconditioner::getPreconditionerState() const {
  Teuchos::RCP<const Preconditioner> bp = rcp_dynamic_cast<const Preconditioner>(preconObj_);

  if (bp != Teuchos::null) return bp->getStateObject();

  return Teuchos::null;
}

Teuchos::RCP<const Thyra::LinearOpBase<ST> > TpetraBlockPreconditioner::extractLinearOp(
    const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> >& A) const {
  // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
  const RCP<const TpetraOperatorWrapper>& tow = rcp_dynamic_cast<const TpetraOperatorWrapper>(A);

  // if it is an EpetraOperatorWrapper, then get the Thyra operator
  if (tow != Teuchos::null) return tow->getThyraOp();

  // otherwise wrap it up as a thyra operator
  return Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A->getRangeMap()), A);
}

Teuchos::RCP<const MappingStrategy> TpetraBlockPreconditioner::extractMappingStrategy(
    const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> >& A) const {
  // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
  const RCP<const TpetraOperatorWrapper>& tow = rcp_dynamic_cast<const TpetraOperatorWrapper>(A);

  // if it is an EpetraOperatorWrapper, then get the Thyra operator
  if (tow != Teuchos::null) return tow->getMapStrategy();

  // otherwise wrap it up as a thyra operator
  RCP<const Tpetra::Map<LO, GO, NT> > range  = A->getRangeMap();
  RCP<const Tpetra::Map<LO, GO, NT> > domain = A->getDomainMap();
  return rcp(
      new BasicMappingStrategy(range, domain, *Thyra::convertTpetraToThyraComm(range->getComm())));
}

}  // namespace TpetraHelpers
}  // end namespace Teko
