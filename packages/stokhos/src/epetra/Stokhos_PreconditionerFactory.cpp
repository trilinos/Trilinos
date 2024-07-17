// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_PreconditionerFactory.hpp"
#include "Stokhos_IfpackPreconditionerFactory.hpp"
#include "Stokhos_MLPreconditionerFactory.hpp"
#include "Teuchos_Assert.hpp"

Stokhos::PreconditionerFactory::
PreconditionerFactory(const std::string& prec_name,
		      const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  if (prec_name == "Ifpack")
    factory = 
      Teuchos::rcp(new Stokhos::IfpackPreconditionerFactory(params));
  else if (prec_name == "ML")
    factory = 
      Teuchos::rcp(new Stokhos::MLPreconditionerFactory(params));
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Unknown preconditioner type " << prec_name
		       << ".  Valid choices are \"Ifpack\" and \"ML\".");
}

Teuchos::RCP<Epetra_Operator> 
Stokhos::PreconditionerFactory::
compute(const Teuchos::RCP<Epetra_Operator>& mat, bool compute_prec)
{
  return factory->compute(mat, compute_prec);
}

void
Stokhos::PreconditionerFactory::
recompute(const Teuchos::RCP<Epetra_Operator>& mat,
	  const Teuchos::RCP<Epetra_Operator>& prec)
{
  factory->recompute(mat, prec);
}
