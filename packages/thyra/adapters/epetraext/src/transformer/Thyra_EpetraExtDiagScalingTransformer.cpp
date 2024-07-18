// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_EpetraExtDiagScalingTransformer.hpp"
#include "Thyra_MultipliedLinearOpBase.hpp"
#include "Thyra_DiagonalLinearOpBase.hpp"
#include "Thyra_ScaledAdjointLinearOpBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"


namespace Thyra {


// Overridden from LinearOpTransformerBase


bool EpetraExtDiagScalingTransformer::isCompatible(
      const LinearOpBase<double> &op_in) const
{
   using Thyra::unwrap;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;
   using Teuchos::dyn_cast;

   const MultipliedLinearOpBase<double> &multi_op =
         dyn_cast<const MultipliedLinearOpBase<double> >(op_in);

   // this operation must only have two operands
   if(multi_op.numOps()!=2)
      return false;

   double scalar = 0.0;
   EOpTransp transp = NOTRANS;

   // get properties of first operator: Transpose, scaler multiply...etc
   RCP<const LinearOpBase<double> > A;
   unwrap( multi_op.getOp(0), &scalar, &transp, &A );
   if(transp!=NOTRANS) return false;

   // get properties of second operator: Transpose, scaler multiply...etc
   RCP<const LinearOpBase<double> > B;
   unwrap( multi_op.getOp(1), &scalar, &transp, &B );
   if(transp!=NOTRANS) return false;
 
   // see if exactly one operator is a diagonal linear op
   RCP<const DiagonalLinearOpBase<double> > dA 
      = rcp_dynamic_cast<const DiagonalLinearOpBase<double> >(A);
   RCP<const DiagonalLinearOpBase<double> > dB 
      = rcp_dynamic_cast<const DiagonalLinearOpBase<double> >(B);

   if(dA==Teuchos::null && dB!=Teuchos::null)
      return true;

   if(dA!=Teuchos::null && dB==Teuchos::null)
      return true;

   return false;
}


RCP<LinearOpBase<double> >
EpetraExtDiagScalingTransformer::createOutputOp() const
{
   return nonconstEpetraLinearOp();
}


void EpetraExtDiagScalingTransformer::transform(
   const LinearOpBase<double> &op_in,
   const Ptr<LinearOpBase<double> > &op_inout) const
{
   using Thyra::unwrap;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;
   using Teuchos::dyn_cast;

   //
   // A) Get the component Thyra objects for M = op(A) * op(B)
   //
 
   const MultipliedLinearOpBase<double> &multi_op =
         dyn_cast<const MultipliedLinearOpBase<double> >(op_in);

   TEUCHOS_ASSERT(multi_op.numOps()==2);
 
   // get properties of first operator: Transpose, scaler multiply...etc
   const RCP<const LinearOpBase<double> > op_A = multi_op.getOp(0);
   double A_scalar = 0.0;
   EOpTransp A_transp = NOTRANS;
   RCP<const LinearOpBase<double> > A;
   unwrap( op_A, &A_scalar, &A_transp, &A );
   TEUCHOS_ASSERT(A_transp==NOTRANS || A_transp==CONJTRANS); // sanity check
 
   // get properties of third operator: Transpose, scaler multiply...etc
   const RCP<const LinearOpBase<double> > op_B = multi_op.getOp(1);
   double B_scalar = 0.0;
   EOpTransp B_transp = NOTRANS;
   RCP<const LinearOpBase<double> > B;
   unwrap( op_B, &B_scalar, &B_transp, &B );
   TEUCHOS_ASSERT(B_transp==NOTRANS || B_transp==CONJTRANS); // sanity check

   //
   // B) Extract out the CrsMatrix and Diagonal objects
   //

   // see if exactly one operator is a diagonal linear op
   RCP<const DiagonalLinearOpBase<double> > dA 
      = rcp_dynamic_cast<const DiagonalLinearOpBase<double> >(A);
   RCP<const DiagonalLinearOpBase<double> > dB 
      = rcp_dynamic_cast<const DiagonalLinearOpBase<double> >(B);

   RCP<const Epetra_CrsMatrix> epetra_A;
   RCP<const Epetra_CrsMatrix> epetra_B;
   if(dA==Teuchos::null) {
      bool exactly_one_op_must_be_diagonal__dB_neq_null = dB!=Teuchos::null;
      TEUCHOS_ASSERT(exactly_one_op_must_be_diagonal__dB_neq_null);

      // convert second operator to an Epetra_CrsMatrix
      epetra_A = rcp_dynamic_cast<const Epetra_CrsMatrix>(get_Epetra_Operator(*A), true);
   }
   else if(dB==Teuchos::null) {
      bool exactly_one_op_must_be_diagonal__dA_neq_null = dA!=Teuchos::null;
      TEUCHOS_ASSERT(exactly_one_op_must_be_diagonal__dA_neq_null);

      // convert second operator to an Epetra_CrsMatrix
      epetra_B = rcp_dynamic_cast<const Epetra_CrsMatrix>(get_Epetra_Operator(*B), true);
   }
   else {
      bool exactly_one_op_must_be_diagonal=false;
      TEUCHOS_ASSERT(exactly_one_op_must_be_diagonal);
   }

   TEUCHOS_ASSERT( A_transp == NOTRANS ); // ToDo: Handle the transpose
   TEUCHOS_ASSERT( B_transp == NOTRANS ); // ToDo: Handle the transpose


   // get or allocate an object to use as the destination
   EpetraLinearOp & thyra_epetra_op_inout = dyn_cast<EpetraLinearOp>(*op_inout);
   RCP<Epetra_CrsMatrix>  epetra_op =
         rcp_dynamic_cast<Epetra_CrsMatrix>(thyra_epetra_op_inout.epetra_op());

   bool rightScale = dB!=Teuchos::null;

   if(rightScale) {
      RCP<const Epetra_Vector> v = get_Epetra_Vector(epetra_A->OperatorDomainMap(),dB->getDiag());

      // if needed allocate a new operator, otherwise use old one assuming
      // constant sparsity
      if(epetra_op==Teuchos::null)
         epetra_op = Teuchos::rcp(new Epetra_CrsMatrix(*epetra_A));
      else
         *epetra_op = *epetra_A;
      epetra_op->RightScale(*v);
   }
   else {
      RCP<const Epetra_Vector> v = get_Epetra_Vector(epetra_B->OperatorRangeMap(),dA->getDiag());

      // if needed allocate a new operator, otherwise use old one assuming
      // constant sparsity
      if(epetra_op==Teuchos::null)
         epetra_op = Teuchos::rcp(new Epetra_CrsMatrix(*epetra_B));
      else
         *epetra_op = *epetra_B;
      epetra_op->LeftScale(*v);
   }

   epetra_op->Scale(A_scalar*B_scalar);

   // set output operator to use newly create epetra_op
   thyra_epetra_op_inout.initialize(epetra_op);
}


} // namespace Thyra
