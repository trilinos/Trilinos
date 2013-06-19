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

#include "Teuchos_TestForException.hpp"
#include "Stokhos_MonomialProjGramSchmidtPCEBasis.hpp"
#include "Stokhos_MonomialProjGramSchmidtPCEBasis2.hpp"
#include "Stokhos_MonomialGramSchmidtPCEBasis.hpp"
#include "Stokhos_ProductLanczosPCEBasis.hpp"
#include "Stokhos_ProductLanczosGramSchmidtPCEBasis.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::ReducedBasisFactory<ordinal_type, value_type>::
ReducedBasisFactory(
  const Teuchos::ParameterList& params_) :
  params(params_),
  reduction_method(params.get("Reduced Basis Method", 
			      "Monomial Proj Gram-Schmidt"))
{
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::ReducedPCEBasis<ordinal_type, value_type> >
Stokhos::ReducedBasisFactory<ordinal_type, value_type>::
createReducedBasis(
  ordinal_type p,
  const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
  const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
  const Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk) const
{
  // Compute reduced basis
  Teuchos::RCP<Stokhos::ReducedPCEBasis<ordinal_type, value_type> > red_basis;

  if (reduction_method == "Monomial Proj Gram-Schmidt")
    red_basis = Teuchos::rcp(new Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type,value_type>(p, pce, quad, params));

  else if (reduction_method == "Monomial Proj Gram-Schmidt2")
    red_basis = Teuchos::rcp(new Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type,value_type>(p, pce, quad, params));

  else if (reduction_method == "Monomial Gram-Schmidt")
    red_basis = Teuchos::rcp(new Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type,value_type>(p, pce, quad, params));

  else if (reduction_method == "Product Lanczos")
    red_basis = Teuchos::rcp(new Stokhos::ProductLanczosPCEBasis<ordinal_type,value_type>(p, pce, quad, Cijk, params));

  else if (reduction_method == "Product Lanczos Gram-Schmidt")
    red_basis = Teuchos::rcp(new Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type,value_type>(p, pce, quad, Cijk, params));

  else
    TEUCHOS_TEST_FOR_EXCEPTION(
	true, std::logic_error, 
	"Invalid reduced basis method " << reduction_method);

  return red_basis;
}
