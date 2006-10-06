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
#include "FEApp_QuadratureFactory.hpp"
#include "FEApp_GaussianQuadrature2.hpp"

FEApp::QuadratureFactory::QuadratureFactory(
	    const Teuchos::RefCountPtr<Teuchos::ParameterList>& quadParams_) :
  quadParams(quadParams_)
{
}

Teuchos::RefCountPtr<FEApp::AbstractQuadrature>
FEApp::QuadratureFactory::create()
{
  Teuchos::RefCountPtr<FEApp::AbstractQuadrature> strategy;

  std::string& method = quadParams->get("Method", "Gaussian");
  if (method == "Gaussian") {
    int num_points = quadParams->get("Num Points", 2);

    if (num_points == 2) {
      strategy = Teuchos::rcp(new FEApp::GaussianQuadrature2);
    }
    else {
      TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
			 std::endl << 
			 "Error!  Number of quadrature points = " << 
			 num_points << 
			 " is not supported for Gaussian quadrature!" << 
			 std::endl << "Supplied parameter list is " << 
			 std::endl << *quadParams);
    }
  }
  else {
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		       std::endl << 
		       "Error!  Unknown quadrature method " << method << 
		       "!" << std::endl << "Supplied parameter list is " << 
		       std::endl << *quadParams);
  }

  return strategy;
}
