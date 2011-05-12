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
#include "Stokhos_MLPreconditionerFactory.hpp"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_TestForException.hpp"
#ifdef HAVE_STOKHOS_ML
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#endif

Stokhos::MLPreconditionerFactory::
MLPreconditionerFactory(const Teuchos::RCP<Teuchos::ParameterList>& p) :
  precParams(p) 
{
#ifdef HAVE_STOKHOS_ML
  // Set default parameters
  if (precParams->isType<std::string>("default values")) {
    Teuchos::ParameterList ml_defaults;
    ML_Epetra::SetDefaults(precParams->get<std::string>("default values"), 
			   ml_defaults);
    precParams->setParametersNotAlreadySet(ml_defaults);
  }
#endif
}

Teuchos::RCP<Epetra_Operator> 
Stokhos::MLPreconditionerFactory::
compute(const Teuchos::RCP<Epetra_Operator>& op, bool compute_prec) {
#ifdef HAVE_STOKHOS_ML
  Teuchos::RCP<Epetra_RowMatrix> mat = 
    Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(op, true);
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> ml_prec =
    Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*mat, *precParams,
							 compute_prec));
  //ml_prec->PrintUnused(0);
  return ml_prec;
#else
  TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Stokhos::MLPreconditionerFactory is available " <<
		     "only with configured with ML support!");
  return Teuchos::null;
#endif // HAVE_STOKHOS_ML
}

void
Stokhos::MLPreconditionerFactory::
recompute(const Teuchos::RCP<Epetra_Operator>& op,
	  const Teuchos::RCP<Epetra_Operator>& prec_op) {
#ifdef HAVE_STOKHOS_ML
  // Copy matrix represented by "op" into underlying matrix in ML
  Teuchos::RCP<Epetra_CrsMatrix> mat = 
    Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(op, true);
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> ml_prec =
    Teuchos::rcp_dynamic_cast<ML_Epetra::MultiLevelPreconditioner>(prec_op);
  const Epetra_RowMatrix& prec_mat = ml_prec->RowMatrix();
  const Epetra_CrsMatrix& prec_crs_mat = 
    dynamic_cast<const Epetra_CrsMatrix&>(prec_mat);
  Epetra_CrsMatrix& non_const_prec_crs_mat = 
    const_cast<Epetra_CrsMatrix&>(prec_crs_mat);
  non_const_prec_crs_mat = *mat;

  // Compute preconditioner
  ml_prec->ComputePreconditioner();
#else
  TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Stokhos::MLPreconditionerFactory is available " <<
		     "only with configured with ML support!");
#endif // HAVE_STOKHOS_ML
}
