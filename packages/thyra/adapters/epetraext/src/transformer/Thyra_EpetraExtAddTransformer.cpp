// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_EpetraExtAddTransformer.hpp"
#include "Thyra_AddedLinearOpBase.hpp"
#include "Thyra_ScaledAdjointLinearOpBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DiagonalLinearOpBase.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_IdentityLinearOpBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_Assert.hpp"
#include "EpetraExt_ConfigDefs.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_MMHelpers.h"
#include "EpetraExt_Transpose_RowMatrix.h"


#include "EpetraExt_RowMatrixOut.h"


namespace Thyra {


// Overridden from LinearOpTransformerBase


bool EpetraExtAddTransformer::isCompatible(
      const LinearOpBase<double> &/* op_in */) const
{
   TEUCHOS_TEST_FOR_EXCEPT(true);
   TEUCHOS_UNREACHABLE_RETURN(false);
}


RCP<LinearOpBase<double> >
EpetraExtAddTransformer::createOutputOp() const
{
   return nonconstEpetraLinearOp();
}


void EpetraExtAddTransformer::transform(
   const LinearOpBase<double> &op_in,
   const Ptr<LinearOpBase<double> > &op_inout) const
{
   using Thyra::unwrap;
   using EpetraExt::MatrixMatrix;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;
   using Teuchos::dyn_cast;
 
   //
   // A) Get the component Thyra objects
   //
 
   const AddedLinearOpBase<double> & add_op =
         dyn_cast<const AddedLinearOpBase<double> >(op_in);

#ifdef TEUCHOS_DEBUG
   TEUCHOS_ASSERT_EQUALITY( add_op.numOps(), 2 );
#endif

 
   // get properties of first operator: Transpose, scaler multiply...etc
   const RCP<const LinearOpBase<double> > op_A = add_op.getOp(0);
   double A_scalar = 0.0;
   EOpTransp A_transp = NOTRANS;
   RCP<const LinearOpBase<double> > A;
   unwrap( op_A, &A_scalar, &A_transp, &A );
   TEUCHOS_ASSERT(A_transp==NOTRANS || A_transp==CONJTRANS); // sanity check
 
   // get properties of third operator: Transpose, scaler multiply...etc
   const RCP<const LinearOpBase<double> > op_B = add_op.getOp(1);
   double B_scalar = 0.0;
   EOpTransp B_transp = NOTRANS;
   RCP<const LinearOpBase<double> > B;
   unwrap( op_B, &B_scalar, &B_transp, &B );
   TEUCHOS_ASSERT(B_transp==NOTRANS || B_transp==CONJTRANS); // sanity check

   //
   // B) Extract out the Epetra_CrsMatrix objects and the vector
   //

   // first makre sure identity operators are represented as digaonal vectors
   if(rcp_dynamic_cast<const Thyra::IdentityLinearOpBase<double> >(A)!=Teuchos::null) {
     RCP<Thyra::VectorBase<double> > d = Thyra::createMember(A->domain(), "d");
     Thyra::V_S( d.ptr(), 1.0 ); // ToDo: Set ton != 1.0 and generalize
     A = Thyra::diagonal(d);
   }
   if(rcp_dynamic_cast<const Thyra::IdentityLinearOpBase<double> >(B)!=Teuchos::null) {
     RCP<Thyra::VectorBase<double> > d = Thyra::createMember(B->domain(), "d");
     Thyra::V_S( d.ptr(), 1.0 ); // ToDo: Set ton != 1.0 and generalize
     B = Thyra::diagonal(d);
   }
      

   // see if exactly one operator is a diagonal linear op
   RCP<const DiagonalLinearOpBase<double> > dA 
      = rcp_dynamic_cast<const DiagonalLinearOpBase<double> >(A);
   RCP<const DiagonalLinearOpBase<double> > dB 
      = rcp_dynamic_cast<const DiagonalLinearOpBase<double> >(B);

   // convert operators to Epetra_CrsMatrix
   RCP<const Epetra_CrsMatrix> epetra_A;
   RCP<const Epetra_CrsMatrix> epetra_B;
   if(dA==Teuchos::null)
      epetra_A = rcp_dynamic_cast<const Epetra_CrsMatrix>(get_Epetra_Operator(*A), true);
   if(dB==Teuchos::null)
      epetra_B = rcp_dynamic_cast<const Epetra_CrsMatrix>(get_Epetra_Operator(*B), true);

   //
   // C) Do the explicit addition
   //

   if(epetra_A!=Teuchos::null && epetra_B!=Teuchos::null) {
  
      // allocate space for final addition: 3 steps
      //   1. Get destination EpetraLinearOp
      //   2. Extract RCP to destination Epetra_CrsMatrix
      //   3. If neccessary, allocate new Epetra_CrsMatrix
      EpetraLinearOp &thyra_epetra_op_inout = dyn_cast<EpetraLinearOp>(*op_inout);
      RCP<Epetra_CrsMatrix>  epetra_op =
            rcp_dynamic_cast<Epetra_CrsMatrix>(thyra_epetra_op_inout.epetra_op());
      Epetra_CrsMatrix * epetra_op_raw = epetra_op.get();
   
      // perform addition
      const int add_epetra_B_err 
         = EpetraExt::MatrixMatrix::Add(*epetra_A,A_transp==CONJTRANS,A_scalar,*epetra_B,B_transp==CONJTRANS,B_scalar,epetra_op_raw);
      if(epetra_op==Teuchos::null)
         epetra_op = Teuchos::rcp(epetra_op_raw);
      
      TEUCHOS_ASSERT_EQUALITY( add_epetra_B_err, 0 );
     
      epetra_op->FillComplete(epetra_A->DomainMap(),epetra_A->RangeMap()); 
   
      // set output operator to use newly create epetra_op
      thyra_epetra_op_inout.initialize(epetra_op);
   }
   else if((dA!=Teuchos::null && epetra_B!=Teuchos::null) ||
           (dB!=Teuchos::null && epetra_A!=Teuchos::null)) { 

      // get unique addition values
      RCP<const Epetra_CrsMatrix> crsMat = (dA!=Teuchos::null) ? epetra_B : epetra_A;
      double matScalar = (dA!=Teuchos::null) ? B_scalar : A_scalar; 
      RCP<const DiagonalLinearOpBase<double> > diag = (dA!=Teuchos::null) ? dA : dB; 
      double diagScalar = (dA!=Teuchos::null) ? A_scalar : B_scalar; 

      TEUCHOS_ASSERT(crsMat!=Teuchos::null);
      TEUCHOS_ASSERT(diag!=Teuchos::null);

      // get or allocate an object to use as the destination
      EpetraLinearOp & thyra_epetra_op_inout = dyn_cast<EpetraLinearOp>(*op_inout);
      RCP<Epetra_CrsMatrix>  epetra_op =
            rcp_dynamic_cast<Epetra_CrsMatrix>(thyra_epetra_op_inout.epetra_op());

      if(epetra_op==Teuchos::null)
         epetra_op = Teuchos::rcp(new Epetra_CrsMatrix(*crsMat));
      else
         *epetra_op = *crsMat;
      
      // grab vector to add to diagonal
      RCP<const Epetra_Vector> v = get_Epetra_Vector(epetra_op->OperatorDomainMap(),diag->getDiag());

      if(matScalar!=1.0)
         epetra_op->Scale(matScalar);

      // grab digaonal from matrix, do summation, then replace the values
      RCP<Epetra_Vector> diagonal = rcp(new Epetra_Vector(epetra_op->OperatorDomainMap()));
      TEUCHOS_TEST_FOR_EXCEPTION(epetra_op->ExtractDiagonalCopy(*diagonal),std::runtime_error,
                                 "Thyra::EpetraExtractAddTransformer::transform ExtractDiagonalCopy failed!");;    
      diagonal->Update(diagScalar,*v,1.0); // no need to scale matrix, already scaled
      TEUCHOS_TEST_FOR_EXCEPTION(epetra_op->ReplaceDiagonalValues(*diagonal),std::runtime_error,
                                 "Thyra::EpetraExtractAddTransformer::transform ReplaceDiagonalValues failed!");;    

      // set output operator to use newly create epetra_op
      thyra_epetra_op_inout.initialize(epetra_op);
   }
   else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                                 "Your case of adding Epetra operators is not yet implemented! Contact the Thyra developers.");
   }
}


} // namespace Thyra
