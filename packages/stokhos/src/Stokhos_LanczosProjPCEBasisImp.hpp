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

#include "Teuchos_TestForException.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_TimeMonitor.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
LanczosProjPCEBasis(
  ordinal_type p,
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& pce,
  const Stokhos::Sparse3Tensor<ordinal_type, value_type>& Cijk,
  bool normalize) :
  RecurrenceBasis<ordinal_type, value_type>("Lanczos-proj PCE", p, normalize),
  pce_sz(pce.basis()->size()),
  Cijk_matrix(pce_sz,pce_sz),
  weights(Teuchos::Copy, 
	  const_cast<value_type*>(pce.basis()->norm_squared().getRawPtr()), 
	  pce_sz),
  u0(pce_sz),
  lanczos_vecs(pce_sz, p+1)
{
  u0[0] = value_type(1);
  
  // Compute matrix
  typedef Stokhos::Sparse3Tensor<ordinal_type, value_type> Cijk_type;
  for (typename Cijk_type::k_iterator k_it = Cijk.k_begin();
       k_it != Cijk.k_end(); ++k_it) {
    ordinal_type k = index(k_it);
    for (typename Cijk_type::kj_iterator j_it = Cijk.j_begin(k_it); 
	 j_it != Cijk.j_end(k_it); ++j_it) {
      ordinal_type j = index(j_it);
      value_type val = 0;
      for (typename Cijk_type::kji_iterator i_it = Cijk.i_begin(j_it);
	   i_it != Cijk.i_end(j_it); ++i_it) {
	ordinal_type i = index(i_it);
	value_type c = value(i_it);
	val += pce[i]*c;
      }
      Cijk_matrix(k,j) = val / weights[k];
    }
  }
  
  // Setup of rest of recurrence basis
  this->setup();
}

template <typename ordinal_type, typename value_type>
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
~LanczosProjPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > 
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
cloneWithOrder(ordinal_type p) const
{
   return 
     Teuchos::rcp(new Stokhos::LanczosProjPCEBasis<ordinal_type,value_type>(
		    p,*this));
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
transformCoeffsFromLanczos(const value_type *in, value_type *out) const
{
  Teuchos::BLAS<ordinal_type, value_type> blas;
  blas.GEMV(Teuchos::NO_TRANS, pce_sz, this->p+1, 
	    value_type(1.0), lanczos_vecs.values(), pce_sz, 
	    in, ordinal_type(1), value_type(0.0), out, ordinal_type(1));
}

template <typename ordinal_type, typename value_type>
bool
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
computeRecurrenceCoefficients(ordinal_type n,
			      Teuchos::Array<value_type>& alpha,
			      Teuchos::Array<value_type>& beta,
			      Teuchos::Array<value_type>& delta,
			      Teuchos::Array<value_type>& gamma) const
{
  Teuchos::Array<value_type> nrm(n);
  vectorspace_type vs(weights);
  operator_type A(Cijk_matrix);

  // Create space to store lanczos vectors -- use lanczos_vecs if 
  // we are requesting p+1 vectors
  Teuchos::RCP<matrix_type> lv;
  if (n == this->p+1)
    lv = Teuchos::rcp(&lanczos_vecs, false);
  else
    lv = Teuchos::rcp(new matrix_type(pce_sz,n));
  lanczos_type::compute(n, vs, A, u0, *lv, alpha, beta, nrm);
  for (ordinal_type i=0; i<n; i++) {
    delta[i] = value_type(1.0);
    gamma[i] = value_type(1.0);
  }

  return false;
}

template <typename ordinal_type, typename value_type> 
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
LanczosProjPCEBasis(ordinal_type p, const LanczosProjPCEBasis& basis) :
  RecurrenceBasis<ordinal_type, value_type>("Lanczos-proj PCE", p, false),
  pce_sz(basis.pce_sz),
  Cijk_matrix(basis.Cijk_matrix),
  weights(basis.weights),
  u0(basis.u0),
  lanczos_vecs(basis.lanczos_vecs)
{
  this->setup();
}
