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
#include "Stokhos_JacobiBasis.hpp"
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
  bool isotropic = basisParams.get("Isotropic", false);

  Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(dimension);
  for (ordinal_type i=0; i<dimension; i++) {
    if (isotropic)
      bases[i] = create1DBasis(basisParams);
    else {
      std::ostringstream ss;
      ss << "Basis " << i;
      Teuchos::ParameterList& bp = basisParams.sublist(ss.str());
      bases[i] = create1DBasis(bp);
    }
  }

  value_type drop = basisParams.get("Cijk Drop Tolerance", 1e-15);
  bool use_old = basisParams.get("Use Old Cijk Algorithm", false);
  basis = 
    Teuchos::rcp(new Stokhos::CompletePolynomialBasis<ordinal_type,value_type>(bases, drop, use_old));
  basisParams.set("Stochastic Galerkin Basis", basis);
  
  return basis;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> >
Stokhos::BasisFactory<ordinal_type, value_type>::
create1DBasis(Teuchos::ParameterList& bp)
{
  Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > basis;

  std::string type = bp.get("Type","Legendre");
  ordinal_type order = bp.get("Order", 3);
  bool normalize = bp.get("Normalize", false);
  if (type == "Legendre")
    basis = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(order, normalize));
  else if (type == "Clenshaw-Curtis") {
    bool isotropic = bp.get("Isotropic", false);
    basis = Teuchos::rcp(new Stokhos::ClenshawCurtisLegendreBasis<ordinal_type,value_type>(order, normalize, isotropic));
  }
  else if (type == "Hermite")
    basis = Teuchos::rcp(new Stokhos::HermiteBasis<ordinal_type,value_type>(order, normalize));
  else if (type == "Jacobi") {
    value_type alpha = bp.get<value_type>("Jacobi Alpha");
    value_type beta = bp.get<value_type>("Jacobi Beta");
    basis = Teuchos::rcp(new Stokhos::JacobiBasis<ordinal_type,value_type>(order, alpha, beta, normalize));
  }
  else if (type == "Rys") {
    value_type cut = bp.get("Weight Cut", 1.0);
    basis = Teuchos::rcp(new Stokhos::RysBasis<ordinal_type,value_type>(order, cut, normalize));
  }
  else
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		       std::endl << 
		       "Error!  Stokhos::BasisFactory::create1DBasis():  " <<
		       "Invalid basis type  " << type << std::endl);

  if (bp.isType<int>("Sparse Grid Rule"))
    basis->setSparseGridRule(bp.get<int>("Sparse Grid Rule"));
  if (bp.isType<int>("Sparse Grid Growth Rule"))
    basis->setSparseGridRule(bp.get<int>("Sparse Grid Growth Rule"));

  return basis;
}
