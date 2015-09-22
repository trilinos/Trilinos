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

#include <sstream>
#include "Teuchos_Assert.hpp"
#include "Teuchos_Array.hpp"

#include "Stokhos_LegendreBasis.hpp"
#include "Stokhos_ClenshawCurtisLegendreBasis.hpp"
#include "Stokhos_GaussPattersonLegendreBasis.hpp"
#include "Stokhos_HermiteBasis.hpp"
#include "Stokhos_JacobiBasis.hpp"
#include "Stokhos_RysBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_TensorProductBasis.hpp"
#include "Stokhos_TotalOrderBasis.hpp"
#include "Stokhos_SmolyakBasis.hpp"

template <typename ordinal_type, typename value_type>
Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >
Stokhos::BasisFactory<ordinal_type, value_type>::
create(Teuchos::ParameterList& sgParams)
{
  Teuchos::ParameterList& basisParams = sgParams.sublist("Basis");

  // Check if basis is already there
  Teuchos::RCP< const OrthogPolyBasis<ordinal_type,value_type> > basis
    = basisParams.template get< Teuchos::RCP< const OrthogPolyBasis<ordinal_type,value_type> > >("Stochastic Galerkin Basis", Teuchos::null);
  if (basis != Teuchos::null) 
    return basis;

  ordinal_type dimension = basisParams.get("Dimension", 1);
  bool isotropic = basisParams.get("Isotropic", false);

  Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(dimension);
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

  std::string type = basisParams.get("Multivariate Type", "Complete");
  value_type drop = basisParams.get("Cijk Drop Tolerance", 1e-12);
  std::string ordering = basisParams.get("Coefficient Ordering", "Total");
  
  if (type == "Complete") {
    bool use_old = basisParams.get("Use Old Cijk Algorithm", false);
    basis = 
      Teuchos::rcp(new CompletePolynomialBasis<ordinal_type,value_type>(
		     bases, drop, use_old));
  }

  else if (type == "Tensor Product") {
    if (ordering == "Total")
       basis = 
	 Teuchos::rcp(new TensorProductBasis<ordinal_type,value_type,TotalOrderLess<MultiIndex<ordinal_type> > >(bases, drop));
    else if (ordering == "Lexicographical")
      basis = 
	 Teuchos::rcp(new TensorProductBasis<ordinal_type,value_type,LexographicLess<MultiIndex<ordinal_type> > >(bases, drop));
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
	true, Teuchos::Exceptions::InvalidParameter,
	std::endl << "Invalid coefficient ordering  " << ordering << std::endl);
  }

  else if (type == "Total Order") {
    if (ordering == "Total")
       basis = 
	 Teuchos::rcp(new TotalOrderBasis<ordinal_type,value_type,TotalOrderLess<MultiIndex<ordinal_type> > >(bases, drop));
    else if (ordering == "Lexicographical")
      basis = 
	 Teuchos::rcp(new TotalOrderBasis<ordinal_type,value_type,LexographicLess<MultiIndex<ordinal_type> > >(bases, drop));
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
	true, Teuchos::Exceptions::InvalidParameter,
	std::endl << "Invalid coefficient ordering  " << ordering << std::endl);
  }

  else if (type == "Smolyak") {
    ordinal_type order = basisParams.template get<ordinal_type>("Order");
    TotalOrderIndexSet<ordinal_type> index_set(dimension, order);
    if (ordering == "Total")
       basis = 
	 Teuchos::rcp(new SmolyakBasis<ordinal_type,value_type,TotalOrderLess<MultiIndex<ordinal_type> > >(bases, index_set, drop));
    else if (ordering == "Lexicographical")
      basis = 
	Teuchos::rcp(new SmolyakBasis<ordinal_type,value_type,LexographicLess<MultiIndex<ordinal_type> > >(bases, index_set, drop));
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
	true, Teuchos::Exceptions::InvalidParameter,
	std::endl << "Invalid coefficient ordering  " << ordering << std::endl);
  }

  else
    TEUCHOS_TEST_FOR_EXCEPTION(
	true, Teuchos::Exceptions::InvalidParameter,
	std::endl << "Invalid multivariate basis type " << type << std::endl);


  basisParams.set("Stochastic Galerkin Basis", basis);
  
  return basis;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> >
Stokhos::BasisFactory<ordinal_type, value_type>::
create1DBasis(Teuchos::ParameterList& bp)
{
  Teuchos::RCP<OneDOrthogPolyBasis<ordinal_type,value_type> > basis;

  std::string type = bp.get("Type","Legendre");
  ordinal_type order = bp.get("Order", 3);
  bool normalize = bp.get("Normalize", false);
  bool isotropic = bp.get("Isotropic", false);

  std::string growth_string = bp.get("Growth Policy", "Slow");
  GrowthPolicy growth;
  if (growth_string == "Slow")
    growth = SLOW_GROWTH;
  else if (growth_string == "Moderate")
    growth = MODERATE_GROWTH;
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
	true, Teuchos::Exceptions::InvalidParameter,
	std::endl << "Invalid growth policy  " << growth_string << std::endl);

  if (type == "Legendre")
    basis = Teuchos::rcp(new LegendreBasis<ordinal_type,value_type>(order, normalize, growth));
  else if (type == "Clenshaw-Curtis") {
    basis = Teuchos::rcp(new ClenshawCurtisLegendreBasis<ordinal_type,value_type>(order, normalize, isotropic));
  }
  else if (type == "Gauss-Patterson") {
    basis = Teuchos::rcp(new GaussPattersonLegendreBasis<ordinal_type,value_type>(order, normalize, isotropic));
  }
  else if (type == "Hermite")
    basis = Teuchos::rcp(new HermiteBasis<ordinal_type,value_type>(order, normalize, growth));
  else if (type == "Jacobi") {
    value_type alpha = bp.get<value_type>("Jacobi Alpha");
    value_type beta = bp.get<value_type>("Jacobi Beta");
    basis = Teuchos::rcp(new JacobiBasis<ordinal_type,value_type>(order, alpha, beta, normalize, growth));
  }
  else if (type == "Rys") {
    value_type cut = bp.get("Weight Cut", 1.0);
    basis = Teuchos::rcp(new RysBasis<ordinal_type,value_type>(order, cut, normalize, growth));
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, Teuchos::Exceptions::InvalidParameter,
      std::endl << "Invalid basis type  " << type << std::endl);

  return basis;
}
