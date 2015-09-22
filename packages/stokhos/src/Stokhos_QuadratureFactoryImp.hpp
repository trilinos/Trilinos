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

#include "Teuchos_Assert.hpp"

#include "Stokhos_BasisFactory.hpp"
#include "Stokhos_TensorProductQuadrature.hpp"
#include "Stokhos_SparseGridQuadrature.hpp"
#include "Stokhos_ProductBasis.hpp"

template <typename ordinal_type, typename value_type>
Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >
Stokhos::QuadratureFactory<ordinal_type, value_type>::
create(Teuchos::ParameterList& sgParams)
{
  // Check if quadrature is already there
  Teuchos::ParameterList& quadParams = sgParams.sublist("Quadrature");
  Teuchos::RCP<const Stokhos::Quadrature<ordinal_type,value_type> > quad = 
    quadParams.template get< Teuchos::RCP<const Stokhos::Quadrature<ordinal_type,value_type> > >("Stochastic Galerkin Quadrature", Teuchos::null);
  if (quad != Teuchos::null)
    return quad;

  // Get basis
  Teuchos::ParameterList& basisParams = sgParams.sublist("Basis");
  Teuchos::RCP< const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > basis;
  if (basisParams.template isType< Teuchos::RCP< const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > >("Stochastic Galerkin Basis"))
    basis = basisParams.template get< Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > >("Stochastic Galerkin Basis");
  else
    basis = Stokhos::BasisFactory<ordinal_type,value_type>::create(sgParams);
  Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > product_basis = Teuchos::rcp_dynamic_cast<const Stokhos::ProductBasis<ordinal_type,value_type> >(basis, true);

  // Create quadrature
  std::string quad_type = quadParams.get("Type", "Tensor Product");
  if (quad_type == "Tensor Product") {
    if (quadParams.isType<ordinal_type>("Quadrature Order")) {
      ordinal_type order = quadParams.get<ordinal_type>("Quadrature Order");
      quad = 
	Teuchos::rcp(new Stokhos::TensorProductQuadrature<ordinal_type,value_type>(product_basis, order));
    }
    else {
      quad = 
	Teuchos::rcp(new Stokhos::TensorProductQuadrature<ordinal_type,value_type>(product_basis));
    }
  }
  else if (quad_type == "Sparse Grid") {
#ifdef HAVE_STOKHOS_DAKOTA
    ordinal_type level = quadParams.get("Sparse Grid Level", 0);
    value_type dup_tol = quadParams.get("Duplicate Tolerance", 1e-12);
    ordinal_type growth = quadParams.get<ordinal_type>(
      "Growth Rule", Pecos::SLOW_RESTRICTED_GROWTH);
    quad = 
      Teuchos::rcp(new Stokhos::SparseGridQuadrature<ordinal_type,value_type>(
		     product_basis, level, dup_tol, growth));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		       std::endl << 
		       "Error!  Stokhos::QuadratureFactory::create():  " <<
		       "Sparse Grid Quadrature requires Dakota!" << std::endl);
#endif
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		       std::endl << 
		       "Error!  Stokhos::QuadratureFactory::create():  " <<
		       "Invalid quadrature type  " << quad_type << std::endl);

  quadParams.set("Stochastic Galerkin Quadrature", quad);
  return quad;
}
