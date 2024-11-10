// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace Thyra {

template<>
Teuchos::RCP<Epetra_Operator>
get_Epetra_Operator( LinearOpBase<double> &op )
{
  EpetraLinearOp &thyra_epetra_op = Teuchos::dyn_cast<EpetraLinearOp>(op);
  return thyra_epetra_op.epetra_op();
}

template<>
Teuchos::RCP<const Epetra_Operator>
get_Epetra_Operator( const LinearOpBase<double> &op )
{
  const EpetraLinearOp &thyra_epetra_op = Teuchos::dyn_cast<const EpetraLinearOp>(op);
  return thyra_epetra_op.epetra_op();
}

} // namespace Thyra
