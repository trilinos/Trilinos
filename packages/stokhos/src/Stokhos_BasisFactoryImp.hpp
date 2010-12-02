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

#include <sstream>
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"

#include "Stokhos_LegendreBasis.hpp"
#include "Stokhos_ClenshawCurtisLegendreBasis.hpp"
#include "Stokhos_HermiteBasis.hpp"
#include "Stokhos_RysBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"

template <typename ordinal_type, typename value_type>
Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >
Stokhos::BasisFactory<ordinal_type, value_type>::
create(Teuchos::ParameterList& sgParams)
{
  Teuchos::ParameterList& basisParams = sgParams.sublist("Basis");

  // Check if basis is already there
  Teuchos::RCP< const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > basis
    = basisParams.template get< Teuchos::RCP< const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > >("Stochastic Galerkin Basis", Teuchos::null);
  if (basis != Teuchos::null) 
    return basis;

  ordinal_type dimension = basisParams.get("Dimension", 1);
  Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(dimension);
  for (ordinal_type i=0; i<dimension; i++) {
    std::ostringstream ss;
    ss << "Basis " << i;
    Teuchos::ParameterList& bp = basisParams.sublist(ss.str());
    std::string type = bp.get("Type","Legendre");
    ordinal_type order = bp.get("Order", 3);
    bool normalize = bp.get("Normalize", false);
    if (type == "Legendre")
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(order, normalize));
    else if (type == "Clenshaw-Curtis") {
      bool isotropic = bp.get("Isotropic", false);
      bases[i] = Teuchos::rcp(new Stokhos::ClenshawCurtisLegendreBasis<ordinal_type,value_type>(order, normalize, isotropic));
    }
    else if (type == "Hermite")
      bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<ordinal_type,value_type>(order, normalize));
    else if (type == "Rys") {
      value_type cut = bp.get("Weight Cut", 1.0);
      bases[i] = Teuchos::rcp(new Stokhos::RysBasis<ordinal_type,value_type>(order, cut, normalize));
    }
    else
      TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
			 std::endl << 
			 "Error!  Stokhos::BasisFactory::create():  " <<
			 "Invalid basis type  " << type << std::endl);
    
    
  }
  basis = 
    Teuchos::rcp(new Stokhos::CompletePolynomialBasis<ordinal_type,value_type>(bases));
  basisParams.set("Stochastic Galerkin Basis", basis);
  
  return basis;
}
