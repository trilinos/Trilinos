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
  RecurrenceBasis<ordinal_type, value_type>("Lanczos-proj PCE", p, normalize)
{
  ordinal_type sz = pce.basis()->size();
  u0.resize(sz);
  //u0.putScalar(value_type(1));
  u0[0] = value_type(1);
  weights = pce.basis()->norm_squared();
  Cijk_matrix.shape(sz,sz);

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
  
  // Compute coefficients via Lanczos
  lanczos_type::compute(sz, p+1, Cijk_matrix, weights, u0, lanczos_vecs,
			this->alpha, this->beta, this->norms);
  for (ordinal_type i=0; i<=p; i++)
    this->delta[i] = value_type(1.0);
  
  // Setup rest of recurrence basis
  //this->setup();
  this->gamma[0] = value_type(1);
  for (ordinal_type k=1; k<=p; k++) {
    this->gamma[k] = value_type(1);
  } 

  //If you want normalized polynomials, set gamma and reset the norms to 1.
  if( normalize ) {
    for (ordinal_type k=0; k<=p; k++) {
      this->gamma[k] = value_type(1)/std::sqrt(this->norms[k]);
      this->norms[k] = value_type(1);
    }
  }
}

template <typename ordinal_type, typename value_type>
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
~LanczosProjPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
getQuadPoints(ordinal_type quad_order,
	      Teuchos::Array<value_type>& quad_points,
	      Teuchos::Array<value_type>& quad_weights,
	      Teuchos::Array< Teuchos::Array<value_type> >& quad_values) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::LanczosProjPCEBasis -- compute Gauss points");
#endif

  // Call base class
  ordinal_type num_points = 
    static_cast<ordinal_type>(std::ceil((quad_order+1)/2.0));

  // We can't reliably generate quadrature points of order > 2*p
  //std::cout << "quad_order = " << quad_order << ", 2*p = " << 2*this->p << std::endl;
  // if (quad_order > 2*this->p)
  //   quad_order = 2*this->p;
  Stokhos::RecurrenceBasis<ordinal_type,value_type>::getQuadPoints(quad_order, 
								   quad_points, 
								   quad_weights,
								   quad_values);

  // Fill in the rest of the points with zero weight
  if (quad_weights.size() < num_points) {
    ordinal_type old_size = quad_weights.size();
    quad_weights.resize(num_points);
    quad_points.resize(num_points);
    quad_values.resize(num_points);
    for (ordinal_type i=old_size; i<num_points; i++) {
      quad_weights[i] = value_type(0);
      quad_points[i] = quad_points[0];
      quad_values[i].resize(this->p+1);
      evaluateBases(quad_points[i], quad_values[i]);
    }
  }
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
computeRecurrenceCoefficients(ordinal_type n,
			      Teuchos::Array<value_type>& alpha,
			      Teuchos::Array<value_type>& beta,
			      Teuchos::Array<value_type>& delta) const
{
  ordinal_type sz = weights.size();
  Teuchos::Array<value_type> nrm(n);
  Teuchos::Array<vector_type> lv;
  lanczos_type::compute(sz, n, Cijk_matrix, weights, u0, lv,
			alpha, beta, nrm);
  for (ordinal_type i=0; i<n; i++) {
    delta[i] = value_type(1.0);
  }
}

template <typename ordinal_type, typename value_type> 
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
LanczosProjPCEBasis(ordinal_type p, const LanczosProjPCEBasis& basis) :
  RecurrenceBasis<ordinal_type, value_type>("Lanczos-proj PCE", p, false),
  Cijk_matrix(basis.Cijk_matrix),
  weights(basis.weights),
  u0(basis.u0)
{
  // Compute coefficients via Lanczos
  ordinal_type sz = weights.size();
  lanczos_type::compute(sz, p+1, Cijk_matrix, weights, u0, lanczos_vecs,
			this->alpha, this->beta, this->norms);
  for (ordinal_type i=0; i<=p; i++)
    this->delta[i] = value_type(1.0);
  
  // Setup rest of recurrence basis
  //this->setup();
  this->gamma[0] = value_type(1);
  for (ordinal_type k=1; k<=p; k++) {
    this->gamma[k] = value_type(1);
  } 

  //If you want normalized polynomials, set gamma and reset the norms to 1.
  if( basis.normalize ) {
    for (ordinal_type k=0; k<=p; k++) {
      this->gamma[k] = value_type(1)/std::sqrt(this->norms[k]);
      this->norms[k] = value_type(1);
    }
  }
}
