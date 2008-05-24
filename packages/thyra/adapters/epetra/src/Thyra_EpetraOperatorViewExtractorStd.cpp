// ***********************************************************************
// 
//               Thyra: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
  double                     wrappedScalar = 0.0;
  EOpTransp                    wrappedTransp = NOTRANS;
  const LinearOpBase<double> *wrappedFwdOp = NULL;
  ::Thyra::unwrap(fwdOp,&wrappedScalar,&wrappedTransp,&wrappedFwdOp);
  const EpetraLinearOpBase *eFwdOp = NULL;
  if( ! (eFwdOp = dynamic_cast<const EpetraLinearOpBase*>(wrappedFwdOp)) )
    return false;
  return true;
}

void EpetraOperatorViewExtractorStd::getEpetraOpView(
  const Teuchos::RCP<LinearOpBase<double> >   &fwdOp
  ,Teuchos::RCP<Epetra_Operator>              *epetraOp
  ,EOpTransp                                            *epetraOpTransp
  ,EApplyEpetraOpAs                                   *epetraOpApplyAs
  ,EAdjointEpetraOp                                   *epetraOpAdjointSupport
  ,double                                             *epetraOpScalar
  ) const
{
  TEST_FOR_EXCEPT(true);
  // ToDo: Implement once this is needed by just copying what is below and
  // removing the 'const' in the right places!
}

void EpetraOperatorViewExtractorStd::getEpetraOpView(
  const Teuchos::RCP<const LinearOpBase<double> >   &fwdOp
  ,Teuchos::RCP<const Epetra_Operator>              *epetraOp
  ,EOpTransp                                                  *epetraOpTransp
  ,EApplyEpetraOpAs                                         *epetraOpApplyAs
  ,EAdjointEpetraOp                                         *epetraOpAdjointSupport
  ,double                                                   *epetraOpScalar
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(epetraOp==NULL);
  TEST_FOR_EXCEPT(epetraOpTransp==NULL);
  TEST_FOR_EXCEPT(epetraOpApplyAs==NULL);
  TEST_FOR_EXCEPT(epetraOpAdjointSupport==NULL);
#endif // TEUCHOS_DEBUG
  double                                              wrappedFwdOpScalar = 0.0;
  EOpTransp                                             wrappedFwdOpTransp = NOTRANS;
  Teuchos::RCP<const LinearOpBase<double> >   wrappedFwdOp; 
  unwrap(fwdOp,&wrappedFwdOpScalar,&wrappedFwdOpTransp,&wrappedFwdOp);
  Teuchos::RCP<const EpetraLinearOpBase>
    epetraFwdOp = Teuchos::rcp_dynamic_cast<const EpetraLinearOpBase>(wrappedFwdOp,true);
  EOpTransp epetra_epetraOpTransp;
  epetraFwdOp->getEpetraOpView(epetraOp,&epetra_epetraOpTransp,epetraOpApplyAs,epetraOpAdjointSupport);
  *epetraOpTransp = trans_trans(real_trans(epetra_epetraOpTransp),wrappedFwdOpTransp);
  *epetraOpScalar = wrappedFwdOpScalar;
}

} // namespace Thyra
