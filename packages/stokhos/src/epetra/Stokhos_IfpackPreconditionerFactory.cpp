// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_ConfigDefs.h"
#include "Stokhos_IfpackPreconditionerFactory.hpp"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_Assert.hpp"
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
  ifpackPrec->Initialize();  
  if (compute_prec)
    ifpackPrec->Compute();
  return ifpackPrec;
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
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
  ifpackPrec->Compute();
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Stokhos::IfpackPreconditionerFactory is available " <<
		     "only with configured with Ifpack support!");
#endif // HAVE_STOKHOS_IFPACK
}
