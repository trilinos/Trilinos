// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
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
  const EpetraLinearOpBase *eFwdOp = NULL;
  if( !(eFwdOp = dynamic_cast<const EpetraLinearOpBase*>(wrappedFwdOp)) )
    return false;
  return true;
}


void EpetraOperatorViewExtractorStd::getNonconstEpetraOpView(
  const RCP<LinearOpBase<double> > &fwdOp,
  const Ptr<RCP<Epetra_Operator> > &epetraOp,
  const Ptr<EOpTransp> &epetraOpTransp,
  const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
  const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport,
    const Ptr<double> &epetraOpScalar
  ) const
{
  TEST_FOR_EXCEPT(true);
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
