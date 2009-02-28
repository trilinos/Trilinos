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


#include "Thyra_EpetraExtDiagScaledMatProdTransformer.hpp"
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
#include "EpetraExt_MatrixMatrix.h"
#include "Teuchos_Assert.hpp"


namespace Thyra {


// Overridden from LinearOpTransformerBase


bool EpetraExtDiagScaledMatProdTransformer::isCompatible(
  const LinearOpBase<double> &op_in) const
{
  TEST_FOR_EXCEPT(true);
  return false;
}


RCP<LinearOpBase<double> >
EpetraExtDiagScaledMatProdTransformer::createOutputOp() const
{
  return nonconstEpetraLinearOp();
}


void EpetraExtDiagScaledMatProdTransformer::transform(
  const LinearOpBase<double> &op_in,
  const Ptr<LinearOpBase<double> > &op_inout
  ) const
{

  using Thyra::unwrap;
  using EpetraExt::MatrixMatrix;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::dyn_cast;

  //
  // A) Get the component Thyra objects for M = op(B) * D * G
  //

  const MultipliedLinearOpBase<double> &multi_op =
    dyn_cast<const MultipliedLinearOpBase<double> >(op_in);
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY( multi_op.numOps(), 3 );
#endif

  const RCP<const LinearOpBase<double> > op_B = multi_op.getOp(0);
  double B_scalar = 0.0;
  EOpTransp B_transp = NOTRANS;
  RCP<const LinearOpBase<double> > B;
  unwrap( op_B, &B_scalar, &B_transp, &B );

  const RCP<const LinearOpBase<double> > D = multi_op.getOp(1);
  const RCP<const VectorBase<double> > d =
    rcp_dynamic_cast<const DiagonalLinearOpBase<double> >(D, true)->getDiag();

  const RCP<const LinearOpBase<double> > G = multi_op.getOp(2);

  //
  // B) Extract out the Epetra_CrsMatrix objects and the vector
  //
  
  const RCP<const Epetra_CrsMatrix> epetra_B =
    rcp_dynamic_cast<const Epetra_CrsMatrix>(get_Epetra_Operator(*B), true);
  
  const RCP<const Epetra_CrsMatrix> epetra_G =
    rcp_dynamic_cast<const Epetra_CrsMatrix>(get_Epetra_Operator(*G), true);

  TEUCHOS_ASSERT( B_transp == NOTRANS ); // ToDo: Handle the transpose
  const Epetra_Map op_inout_row_map = epetra_B->RowMap();

  const Epetra_Map op_inout_col_map = epetra_G->ColMap();
  
  const RCP<const Epetra_Vector> epetra_d =
    get_Epetra_Vector(epetra_G->OperatorRangeMap(), d);

  //
  // C) Do the explicit multiplication
  //
 
  EpetraLinearOp &thyra_epetra_op_inout = dyn_cast<EpetraLinearOp>(*op_inout);
  RCP<Epetra_CrsMatrix>  epetra_op =
    rcp_dynamic_cast<Epetra_CrsMatrix>(thyra_epetra_op_inout.epetra_op());
  if (is_null(epetra_op)) {
    epetra_op = Teuchos::rcp(
      new Epetra_CrsMatrix(::Copy, op_inout_row_map, op_inout_col_map, 0)
      );
    // 2009/02/27: rabartl: Note: Above, the row map must be the right size
    // and distribution and the column map can not be arbitrary.
  }
   
  // ToDo: Implement the case where d != 1.0!

  TEUCHOS_ASSERT_INEQUALITY(
    MatrixMatrix::Multiply(
      *epetra_B, B_transp != NOTRANS,
      *epetra_G, false,
      *epetra_op
      ),
    >=, 0
    );

  thyra_epetra_op_inout.initialize(epetra_op);

}


} // namespace Thyra
