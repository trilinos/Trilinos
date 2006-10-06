// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER
#include "Teuchos_TestForException.hpp"
#include "Sacado_ConfigDefs.h"
#include "FEApp_DiscretizationFactory.hpp"
#include "FEApp_CZeroDiscretization.hpp"

FEApp::DiscretizationFactory::DiscretizationFactory(
	    const Teuchos::RefCountPtr<Teuchos::ParameterList>& discParams_) :
  discParams(discParams_)
{
}

Teuchos::RefCountPtr<FEApp::AbstractDiscretization>
FEApp::DiscretizationFactory::create(
		  const std::vector<double>& coords,
		  unsigned int num_equations,
	          const Teuchos::RefCountPtr<const Epetra_Comm>& epetra_comm)
{
  Teuchos::RefCountPtr<FEApp::AbstractDiscretization> strategy;

  std::string& method = discParams->get("Method", "C Zero");
  if (method == "C Zero") {
    strategy = Teuchos::rcp(new FEApp::CZeroDiscretization(coords, 
							   num_equations, 
							   epetra_comm,
							   discParams));
  }
  else {
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		       std::endl << 
		       "Error!  Unknown discretization method " << method << 
		       "!" << std::endl << "Supplied parameter list is " << 
		       std::endl << *discParams);
  }

  return strategy;
}
