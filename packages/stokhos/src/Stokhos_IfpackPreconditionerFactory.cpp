// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_ConfigDefs.h"
#include "Stokhos_IfpackPreconditionerFactory.hpp"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_TestForException.hpp"
#ifdef HAVE_STOKHOS_IFPACK
#include "Ifpack.h"
#endif

Stokhos::IfpackPreconditionerFactory::
IfpackPreconditionerFactory(const Teuchos::RCP<Teuchos::ParameterList>& p) :
  precParams(p) 
{
}

Teuchos::RCP<Epetra_Operator> 
Stokhos::IfpackPreconditionerFactory::
compute(const Teuchos::RCP<Epetra_Operator>& op, bool compute_prec) {
#ifdef HAVE_STOKHOS_IFPACK
  Teuchos::RCP<Epetra_RowMatrix> mat = 
    Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(op, true);
  Ifpack Factory;
  std::string prec = precParams->get("Ifpack Preconditioner", "ILU");
  int overlap = precParams->get("Overlap", 0);
  Teuchos::RCP<Ifpack_Preconditioner> ifpackPrec = 
    Teuchos::rcp(Factory.Create(prec, mat.get(), overlap));
  ifpackPrec->SetParameters(*precParams);
  int err = ifpackPrec->Initialize();  
  if (compute_prec)
    err = ifpackPrec->Compute();
  return ifpackPrec;
#else
  TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Stokhos::IfpackPreconditionerFactory is available " <<
		     "only with configured with Ifpack support!");
  return Teuchos::null;
#endif // HAVE_STOKHOS_IFPACK
}

void
Stokhos::IfpackPreconditionerFactory::
recompute(const Teuchos::RCP<Epetra_Operator>& op,
	  const Teuchos::RCP<Epetra_Operator>& prec_op) {
#ifdef HAVE_STOKHOS_IFPACK
  // Copy matrix represented by "op" into underlying matrix in Ifpack
  Teuchos::RCP<Epetra_CrsMatrix> mat = 
    Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(op, true);
  Teuchos::RCP<Ifpack_Preconditioner> ifpackPrec = 
    Teuchos::rcp_dynamic_cast<Ifpack_Preconditioner>(prec_op);
  const Epetra_RowMatrix& prec_mat = ifpackPrec->Matrix();
  const Epetra_CrsMatrix& prec_crs_mat = 
    dynamic_cast<const Epetra_CrsMatrix&>(prec_mat);
  Epetra_CrsMatrix& non_const_prec_crs_mat = 
    const_cast<Epetra_CrsMatrix&>(prec_crs_mat);
  non_const_prec_crs_mat = *mat;
  
  // Compute preconditioenr
  int err = ifpackPrec->Compute();
#else
  TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Stokhos::IfpackPreconditionerFactory is available " <<
		     "only with configured with Ifpack support!");
#endif // HAVE_STOKHOS_IFPACK
}
