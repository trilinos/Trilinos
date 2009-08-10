// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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

#include "Teuchos_BLAS.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
GramSchmidtBasis(
   const Teuchos::RCP<const OrthogPolyBasis<ordinal_type,value_type> >& basis_,
   const std::vector< std::vector<value_type> >& points,
   const std::vector<value_type>& weights_,
   const value_type& sparse_tol_) :
  name("Gram Schmidt Basis"),
  basis(basis_),
  weights(weights_),
  sparse_tol(sparse_tol_),
  p(basis->order()),
  d(basis->dimension()),
  sz(basis->size()),
  norms(sz),
  gs_mat(sz,sz),
  basis_pts(sz)
{
  // Get quadrature data
  ordinal_type nqp = weights.size();
  std::vector< std::vector<value_type> > values(nqp);
  for (ordinal_type k=0; k<nqp; k++)
    values[k] = basis->evaluateBases(points[k]);

  // Compute all inner products
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> inner_product(sz,sz);
  inner_product.putScalar(0.0);
  for (ordinal_type i=0; i<sz; i++) {
    for (ordinal_type j=0; j<=i; j++) {
      value_type t = 0.0;
      for (ordinal_type k=0; k<nqp; k++)
	t += weights[k]*values[k][i]*values[k][j];
      inner_product(i,j) = t;
    }
  }

  // Classical Gram-Schmidt algorithm:
  // u_i = v_i - \sum_{j<i} (v_i,u_j)/(u_j,u_j) u_j
  // u_j = \sum_{k<=i} a_{jk} v_k
  // => u_i = v_i - \sum_{j<i}\sum_{k<=j} (v_i,u_j)/(u_j,u_j)*a_{jk}*v_k
  for (ordinal_type i=0; i<sz; i++) {

    // a_{ii} = 1.0
    gs_mat(i,i) = 1.0;

    for (ordinal_type j=0; j<i; j++) {

      // compute t = (v_i,u_j)/(u_j,u_j)
      value_type t = 0.0;
      for (ordinal_type k=0; k<=j; k++)
	t += gs_mat(j,k)*inner_product(i,k);
      t /= norms[j];
      
      // substract contribution to a_{ik}:  t*a_{jk}
      for (ordinal_type k=0; k<=j; k++)
	gs_mat(i,k) -= t*gs_mat(j,k);
    }

    // compute (u_i,u_i) = \sum_{j,k<=i} a_{ij}*a_{ik}*(v_j,v_k)
    value_type nrm = 0.0;
    for (ordinal_type j=0; j<=i; j++) {
      for (ordinal_type k=0; k<=j; k++)
	nrm += gs_mat(i,j)*gs_mat(i,k)*inner_product(j,k);
      for (ordinal_type k=j+1; k<=i; k++)
      	nrm += gs_mat(i,j)*gs_mat(i,k)*inner_product(k,j);
    }
    norms[i] = nrm;

  }

  basis_values.resize(nqp);
  for (ordinal_type k=0; k<nqp; k++) {
    basis_values[k].resize(sz);
    for (ordinal_type i=0; i<sz; i++) {
      value_type t = 0.0;
      for (ordinal_type j=0; j<=i; j++)
	t += gs_mat(i,j)*values[k][j];
      basis_values[k][i] = t;
    }
  }
}

template <typename ordinal_type, typename value_type>
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
~GramSchmidtBasis()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
dimension() const
{
  return d;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
size() const
{
  return sz;
}

template <typename ordinal_type, typename value_type>
const std::vector<value_type>&
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
getTripleProductTensor() const
{ 
  // Compute Cijk = < \Psi_i \Psi_j \Psi_k >
  if (Cijk == Teuchos::null) {
    Cijk = Teuchos::rcp(new Sparse3Tensor<ordinal_type, value_type>(sz));
    // Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> > cijk 
    //   = basis->getTripleProductTensor();
    // for (ordinal_type j=0; j<sz; j++) {
    //   for (ordinal_type i=0; i<sz; i++) {
    // 	for (ordinal_type k=0; k<sz; k++) {
    // 	  value_type t = 0.0;
    // 	  for (ordinal_type n=0; n<=k; n++) {
    // 	    //ordinal_type num_m = Cijk->num_j(n);
    // 	    const Teuchos::Array<ordinal_type>& m_indices = cijk->Jindices(n);
    // 	    ordinal_type num_m = m_indices.size();
    // 	    for (ordinal_type midx=0; midx<num_m; midx++) {
    // 	      ordinal_type m = m_indices[midx];
    // 	      //ordinal_type m = cijk->j_index(n,midx);
    // 	      const Teuchos::Array<ordinal_type>& l_indices = 
    // 		cijk->Iindices(n,midx);
    // 	      ordinal_type num_l = l_indices.size();
    // 	      const Teuchos::Array<value_type>& vals = cijk->values(n,midx);
    // 	      for (ordinal_type lidx=0; lidx<num_l; lidx++) {
    // 		ordinal_type l = l_indices[lidx];
    // 		t += gs_mat(i,l)*gs_mat(j,m)*gs_mat(k,n)*vals[lidx];
    // 	      }
    // 	    }
    // 	  }
    // 	  //if (t > sparse_tol)
    // 	    Cijk->add_term(i,j,k,t);
    // 	}
    //   }
    // }

    ordinal_type nqp = weights.size();
    for (ordinal_type j=0; j<sz; j++) {
      for (ordinal_type i=0; i<sz; i++) {
    	for (ordinal_type k=0; k<sz; k++) {
	  value_type t = 0.0;
	  for (ordinal_type l=0; l<nqp; l++)
	    t += weights[l]*basis_values[l][i]*basis_values[l][j]*basis_values[l][k];
	  if (std::abs(t) > sparse_tol)
	    Cijk->add_term(i,j,k,t);
	}
      }
    }
  }

  return Cijk;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< const Stokhos::Dense3Tensor<ordinal_type, value_type> >
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
getDerivTripleProductTensor() const
{
  TEST_FOR_EXCEPTION(true, std::logic_error,
  		     "Stokhos::GramSchmidtBasis::getDerivTripleProductTensor():"
  		     << "  Method not implemented!");
  return Teuchos::null;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< const Teuchos::SerialDenseMatrix<ordinal_type, value_type> >
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
getDerivDoubleProductTensor() const
{
  TEST_FOR_EXCEPTION(true, std::logic_error,
  		     "Stokhos::GramSchmidtBasis::getDerivDoubleProductTensor():"
  		     << "  Method not implemented!");
  return Teuchos::null;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
projectProduct(ordinal_type i, ordinal_type j, std::vector<value_type>& coeffs) const
{
  TEST_FOR_EXCEPTION(true, std::logic_error,
  		     "Stokhos::GramSchmidtBasis::projectProduct():"
  		     << "  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
projectDerivative(ordinal_type i, std::vector<value_type>& coeffs) const
{
  TEST_FOR_EXCEPTION(true, std::logic_error,
  		     "Stokhos::GramSchmidtBasis::projectDerivative():"
  		     << "  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
evaluateZero(ordinal_type i) const
{
  value_type z = 0.0;
  for (ordinal_type j=0; j<sz; j++)
    z += gs_mat(i,j)*basis->evaluateZero(j);

  return z;
}

template <typename ordinal_type, typename value_type>
const std::vector<value_type>&
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
evaluateBases(const std::vector<value_type>& point) const
{
  const std::vector<value_type>& basis_vals = basis->evaluateBases(point);
  for (ordinal_type i=0; i<sz; i++) {
    value_type t = 0.0;
    for (ordinal_type j=0; j<sz; j++)
      t += gs_mat(i,j)*basis_vals[j];
    basis_pts[i] = t;
  }

  return basis_pts;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
print(std::ostream& os) const
{
  os << "Gram-Schmidt basis of order " << p << ", dimension " << d 
     << ", and size " << sz << ".  Matrix coefficients:\n";
  os << gs_mat << std::endl;
  os << "Basis vector norms (squared):\n\t";
  for (ordinal_type i=0; i<sz; i++)
    os << norms[i] << " ";
  os << "\n";
  os << "Underlying basis:\n";
  os << *basis;
}

template <typename ordinal_type, typename value_type>
std::vector<ordinal_type>
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
getTerm(ordinal_type i) const
{
  TEST_FOR_EXCEPTION(true, std::logic_error,
  		     "Stokhos::GramSchmidtBasis::getTerm():"
  		     << "  Method not implemented!");
  return basis->getTerm(i);
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
getIndex(const std::vector<ordinal_type>& term) const
{
  TEST_FOR_EXCEPTION(true, std::logic_error,
  		     "Stokhos::GramSchmidtBasis::getIndex():"
  		     << "  Method not implemented!");
  return basis->getIndex(term);
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
const std::vector< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> > >&
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
getCoordinateBases() const
{
  return basis->getCoordinateBases();
}

template <typename ordinal_type, typename value_type>
void
Stokhos::GramSchmidtBasis<ordinal_type, value_type>::
transformCoeffs(const value_type *in, value_type *out) const
{
  Teuchos::BLAS<ordinal_type, value_type> blas;
  for (ordinal_type i=0; i<sz; i++)
    out[i] = in[i];
  blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::TRANS,
	    Teuchos::UNIT_DIAG, sz, 1, 1.0, gs_mat.values(), sz, out, sz);
}
