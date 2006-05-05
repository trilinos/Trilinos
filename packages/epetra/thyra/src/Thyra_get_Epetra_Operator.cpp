// @HEADER
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

#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Teuchos_dyn_cast.hpp"

Teuchos::RefCountPtr<Epetra_Operator>
Thyra::get_Epetra_Operator( LinearOpBase<double> &op )
{
  EpetraLinearOp &tsfcore_epetra_op = Teuchos::dyn_cast<EpetraLinearOp>(op);
  return tsfcore_epetra_op.epetra_op();
}

Teuchos::RefCountPtr<const Epetra_Operator>
Thyra::get_Epetra_Operator( const LinearOpBase<double> &op )
{
  const EpetraLinearOp &tsfcore_epetra_op = Teuchos::dyn_cast<const EpetraLinearOp>(op);
  return tsfcore_epetra_op.epetra_op();
}
