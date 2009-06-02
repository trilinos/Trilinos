// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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


#include "Thyra_EpetraExtAddTransformer.hpp"
#include "Thyra_AddedLinearOpBase.hpp"
#include "Thyra_ScaledAdjointLinearOpBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_Assert.hpp"
#include "EpetraExt_ConfigDefs.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_MMHelpers.h"
#include "EpetraExt_Transpose_RowMatrix.h"


#include "EpetraExt_RowMatrixOut.h"


namespace Thyra {


int EE_Add(const Epetra_CrsMatrix & A,
                      bool transposeA,
                      double scalarA,
                      const Epetra_CrsMatrix & B,
                      bool transposeB,
                      double scalarB,
                      Epetra_CrsMatrix * & C)
{
  //
  //This method forms the matrix-matrix sum C = scalarA * op(A) + scalarB * op(B), where

  //A and B should already be Filled. C should be an empty pointer.

  if (!A.Filled() || !B.Filled() ) {
     std::cerr << "EpetraExt::MatrixMatrix::Add ERROR, input matrix A.Filled() or B.Filled() is false,"
               << "they are required to be true. (Result matrix C should be an empty pointer)" << std::endl;
     EPETRA_CHK_ERR(-1);
  }

  Epetra_CrsMatrix * Aprime = 0, * Bprime=0;
  EpetraExt::RowMatrix_Transpose * Atrans = 0,* Btrans = 0;

  //explicit tranpose A formed as necessary
  if( transposeA ) {
     Atrans = new EpetraExt::RowMatrix_Transpose();
     Aprime = &(dynamic_cast<Epetra_CrsMatrix&>(((*Atrans)(const_cast<Epetra_CrsMatrix&>(A)))));
  }
  else
     Aprime = const_cast<Epetra_CrsMatrix*>(&A);

  //explicit tranpose B formed as necessary
  if( transposeB ) {
     Btrans = new EpetraExt::RowMatrix_Transpose();
     Bprime = &(dynamic_cast<Epetra_CrsMatrix&>(((*Btrans)(const_cast<Epetra_CrsMatrix&>(B)))));
  }
  else
     Bprime = const_cast<Epetra_CrsMatrix*>(&B);

  // allocate the new matrix
  C = new Epetra_CrsMatrix(Copy,Aprime->RowMap(),0);

  // build arrays  for easy resuse
  int ierr = 0;
  Epetra_CrsMatrix * Mat[] = { Aprime,Bprime};
  double scalar[] = { scalarA, scalarB};

  // do a loop over each matrix to add
  for(int k=0;k<2;k++) {
     int MaxNumEntries = Mat[k]->MaxNumEntries();
     int NumEntries;
     int * Indices = new int[MaxNumEntries];
     double * Values = new double[MaxNumEntries];
   
     int NumMyRows = Mat[k]->NumMyRows();
     int Row, err;
     int ierr = 0;
   
     //Loop over B's rows and sum into
     for( int i = 0; i < NumMyRows; ++i ) {
        Row = Mat[k]->GRID(i);
        EPETRA_CHK_ERR( Mat[k]->ExtractGlobalRowCopy( Row, MaxNumEntries, NumEntries, Values, Indices));
   
        if( scalar[k] != 1.0 )
           for( int j = 0; j < NumEntries; ++j ) Values[j] *= scalar[k];
   
        err = C->InsertGlobalValues( Row, NumEntries, Values, Indices );
        assert( err == 0 || err == 1 || err == 3 );
        if (err < 0) ierr = err;
     }

     delete [] Indices;
     delete [] Values;
  }

  if( Atrans ) delete Atrans;
  if( Btrans ) delete Btrans;

  return(ierr);
}


// Overridden from LinearOpTransformerBase


bool EpetraExtAddTransformer::isCompatible(
      const LinearOpBase<double> &op_in) const
{
   TEST_FOR_EXCEPT(true);
   return false;
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
   // A) Get the component Thyra objects for M = op(B) * D * G
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
   
  // convert operators to Epetra_CrsMatrix
   const RCP<const Epetra_CrsMatrix> epetra_A =
         rcp_dynamic_cast<const Epetra_CrsMatrix>(get_Epetra_Operator(*A), true);
   const RCP<const Epetra_CrsMatrix> epetra_B =
         rcp_dynamic_cast<const Epetra_CrsMatrix>(get_Epetra_Operator(*B), true);
   
   const Epetra_Map op_inout_row_map 
         = (A_transp==CONJTRANS ? epetra_A->ColMap() : epetra_A->RowMap());
   const Epetra_Map op_inout_col_map 
         = (A_transp==CONJTRANS ? epetra_A->RowMap() : epetra_A->ColMap());
 
   //
   // C) Do the explicit multiplication
   //
  
/*
   // allocate space for final addition: 3 steps
   //   1. Get destination EpetraLinearOp
   //   2. Extract RCP to destination Epetra_CrsMatrix
   //   3. If neccessary, allocate new Epetra_CrsMatrix
   EpetraLinearOp &thyra_epetra_op_inout = dyn_cast<EpetraLinearOp>(*op_inout);
   RCP<Epetra_CrsMatrix>  epetra_op =
         rcp_dynamic_cast<Epetra_CrsMatrix>(thyra_epetra_op_inout.epetra_op());
   if(is_null(epetra_op)) {
      int maxEntries = std::max(epetra_A->GlobalMaxNumEntries(),epetra_B->GlobalMaxNumEntries());
      epetra_op = Teuchos::rcp(
            new Epetra_CrsMatrix(::Copy, op_inout_row_map, maxEntries));
   }
*/

/*
   // perform multiply: Its annoying I have to do this with two adds.
   // Problem is I can't find a copy of CrsMatrices that doesn't call FillComplete.
   // I want this sum operation to allow new entries into the sparse matrix
   // not to be resricted to the sparsity pattern of A or B
   // epetra_op = A_scalar*A
   const int add_epetra_A_err = EpetraExt::MatrixMatrix::Add(
     *epetra_A, A_transp==CONJTRANS, A_scalar, *epetra_op, 0.0 );
   TEUCHOS_ASSERT_EQUALITY( add_epetra_A_err, 0 );

   // epetra_op += B_Scalar*B
   const int add_epetra_B_err = EpetraExt::MatrixMatrix::Add(
     *epetra_B, A_transp==CONJTRANS, B_scalar, *epetra_op, 1.0 );
   TEUCHOS_ASSERT_EQUALITY( add_epetra_B_err, 0 );
*/
   EpetraLinearOp &thyra_epetra_op_inout = dyn_cast<EpetraLinearOp>(*op_inout);
    
   Epetra_CrsMatrix * ptrEpetra_op;
   const int add_epetra_B_err = EE_Add(*epetra_A,A_transp==CONJTRANS,A_scalar,*epetra_B,B_transp==CONJTRANS,B_scalar,ptrEpetra_op);
   TEUCHOS_ASSERT_EQUALITY( add_epetra_B_err, 0 );
   
   RCP<Epetra_CrsMatrix>  epetra_op = Teuchos::rcp(ptrEpetra_op);
   epetra_op->FillComplete();

   // set output operator to use newly create epetra_op
   thyra_epetra_op_inout.initialize(epetra_op);
}


} // namespace Thyra
