// $Id$ 
// $Source$ 
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

#include "Stokhos_PreconditionerFactory.hpp"
#include "Stokhos_IfpackPreconditionerFactory.hpp"
#include "Stokhos_MLPreconditionerFactory.hpp"
#include "Teuchos_TestForException.hpp"

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
    TEST_FOR_EXCEPTION(true, std::logic_error,
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
