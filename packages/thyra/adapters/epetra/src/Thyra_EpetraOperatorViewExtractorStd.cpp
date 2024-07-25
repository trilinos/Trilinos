// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_EpetraOperatorViewExtractorStd.hpp"
#include "Thyra_EpetraLinearOpBase.hpp"
#include "Thyra_ScaledAdjointLinearOpBase.hpp"

class Epetra_Operator;


namespace Thyra {


// Overridden from EpetraOperatorViewExtractorBase


bool EpetraOperatorViewExtractorStd::isCompatible( const LinearOpBase<double> &fwdOp ) const
{
  double wrappedScalar = 0.0;
  EOpTransp wrappedTransp = NOTRANS;
  const LinearOpBase<double> *wrappedFwdOp = NULL;
  ::Thyra::unwrap(fwdOp, &wrappedScalar, &wrappedTransp, &wrappedFwdOp);
  const EpetraLinearOpBase *eFwdOp =
    dynamic_cast<const EpetraLinearOpBase*>(wrappedFwdOp);
  if (!eFwdOp) {
    return false;
  }
  return true;
}


void EpetraOperatorViewExtractorStd::getNonconstEpetraOpView(
  const RCP<LinearOpBase<double> > &/* fwdOp */,
  const Ptr<RCP<Epetra_Operator> > &/* epetraOp */,
  const Ptr<EOpTransp> &/* epetraOpTransp */,
  const Ptr<EApplyEpetraOpAs> &/* epetraOpApplyAs */,
  const Ptr<EAdjointEpetraOp> &/* epetraOpAdjointSupport */,
    const Ptr<double> &/* epetraOpScalar */
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
  // ToDo: Implement once this is needed by just copying what is below and
  // removing the 'const' in the right places!
}


void EpetraOperatorViewExtractorStd::getEpetraOpView(
  const RCP<const LinearOpBase<double> > &fwdOp,
  const Ptr<RCP<const Epetra_Operator> > &epetraOp,
  const Ptr<EOpTransp> &epetraOpTransp,
  const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
  const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport,
    const Ptr<double> &epetraOpScalar
  ) const
{
  using Teuchos::outArg;
  double wrappedFwdOpScalar = 0.0;
  EOpTransp wrappedFwdOpTransp = NOTRANS;
  Teuchos::RCP<const LinearOpBase<double> > wrappedFwdOp; 
  unwrap(fwdOp,&wrappedFwdOpScalar, &wrappedFwdOpTransp, &wrappedFwdOp);
  Teuchos::RCP<const EpetraLinearOpBase> epetraFwdOp =
    Teuchos::rcp_dynamic_cast<const EpetraLinearOpBase>(wrappedFwdOp,true);
  EOpTransp epetra_epetraOpTransp;
  epetraFwdOp->getEpetraOpView(epetraOp, outArg(epetra_epetraOpTransp),
    epetraOpApplyAs, epetraOpAdjointSupport);
  *epetraOpTransp = trans_trans(real_trans(epetra_epetraOpTransp), wrappedFwdOpTransp);
  *epetraOpScalar = wrappedFwdOpScalar;
}


// ToDo: Refactor unwrap(...) to not use raw pointers!


} // namespace Thyra
