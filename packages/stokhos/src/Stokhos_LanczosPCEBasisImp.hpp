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
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
LanczosPCEBasis(
   ordinal_type p,
   const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& pce,
   const Stokhos::Quadrature<ordinal_type, value_type>& quad,
   bool normalize) :
  RecurrenceBasis<ordinal_type, value_type>("Lanczos PCE", p, normalize),
  nqp(quad.size()),
  pce_weights(Teuchos::Copy, 
	      const_cast<value_type*>(quad.getQuadWeights().getRawPtr()), 
	      nqp),
  pce_vals(nqp),
  u0(nqp),
  lanczos_vecs(nqp, p+1),
  fromStieltjesMat(pce.size(), p+1)
{
  // Evaluate PCE at quad points
  const Teuchos::Array< Teuchos::Array<value_type> >& quad_points =
    quad.getQuadPoints();
  const Teuchos::Array< Teuchos::Array<value_type> >& basis_values =
    quad.getBasisAtQuadPoints();
  for (ordinal_type i=0; i<nqp; i++) {
    pce_vals[i] = pce.evaluate(quad_points[i], basis_values[i]);
    u0[i] = value_type(1);
  }

  // Setup rest of basis
  this->setup();

  // Compute transformation matrix back to original basis
  ordinal_type sz = pce.size();
  fromStieltjesMat.putScalar(0.0);
  for (ordinal_type i=0; i<sz; i++) {
    for (ordinal_type j=0; j<=p; j++) {
      for (ordinal_type k=0; k<nqp; k++)
	fromStieltjesMat(i,j) += 
	  pce_weights[k]*lanczos_vecs(k,j)*basis_values[k][i];
      fromStieltjesMat(i,j) /= pce.basis()->norm_squared(i);
    }
  }
}

template <typename ordinal_type, typename value_type>
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
~LanczosPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
getQuadPoints(ordinal_type quad_order,
	      Teuchos::Array<value_type>& quad_points,
	      Teuchos::Array<value_type>& quad_weights,
	      Teuchos::Array< Teuchos::Array<value_type> >& quad_values) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::LanczosPCEBasis -- compute Gauss points");
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
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
cloneWithOrder(ordinal_type p) const
{
   return Teuchos::rcp(new Stokhos::LanczosPCEBasis<ordinal_type,value_type>(
			 p,*this));
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
transformCoeffsFromLanczos(const value_type *in, value_type *out) const
{
  Teuchos::BLAS<ordinal_type, value_type> blas;
  ordinal_type sz = fromStieltjesMat.numRows();
  blas.GEMV(Teuchos::NO_TRANS, sz, this->p+1, 
	    value_type(1.0), fromStieltjesMat.values(), sz, 
	    in, ordinal_type(1), value_type(0.0), out, ordinal_type(1));
}

template <typename ordinal_type, typename value_type>
bool
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
computeRecurrenceCoefficients(ordinal_type n,
			      Teuchos::Array<value_type>& alpha,
			      Teuchos::Array<value_type>& beta,
			      Teuchos::Array<value_type>& delta,
			      Teuchos::Array<value_type>& gamma) const
{
  Teuchos::Array<value_type> nrm(n);
  vectorspace_type vs(pce_weights);
  operator_type A(pce_vals);
  
  // Create space to store lanczos vectors -- use lanczos_vecs if 
  // we are requesting p+1 vectors
  Teuchos::RCP<matrix_type> lv;
  if (n == this->p+1)
    lv = Teuchos::rcp(&lanczos_vecs, false);
  else
    lv = Teuchos::rcp(new matrix_type(nqp,n));

  lanczos_type::compute(n, vs, A, u0, *lv, alpha, beta, nrm);
  for (ordinal_type i=0; i<n; i++) {
    delta[i] = value_type(1.0);
    gamma[i] = value_type(1.0);
  }

  return false;
}

template <typename ordinal_type, typename value_type> 
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
LanczosPCEBasis(ordinal_type p, const LanczosPCEBasis& basis) :
  RecurrenceBasis<ordinal_type, value_type>("Lanczos PCE", p, false),
  nqp(basis.nqp),
  pce_weights(basis.pce_weights),
  pce_vals(basis.pce_vals),
  u0(basis.u0),
  lanczos_vecs(basis.lanczos_vecs)
{
  this->setup();
}
