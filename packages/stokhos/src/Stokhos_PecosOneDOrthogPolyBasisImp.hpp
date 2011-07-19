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

template <typename ordinal_type, typename value_type>
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type, value_type>::
PecosOneDOrthogPolyBasis(
  const Teuchos::RCP<Pecos::OrthogonalPolynomial>& pecosPoly_,
  const std::string& name_, ordinal_type p_) :
  pecosPoly(pecosPoly_),
  name(name_),
  p(p_),
  sparse_grid_rule(Pecos::GOLUB_WELSCH),
  sparse_grid_growth_rule(Pecos::DEFAULT_GROWTH),
  norms(p+1, value_type(0.0))
{
  for (ordinal_type i=0; i<=p; i++)
    norms[i] = pecosPoly->norm_squared(i);
}

template <typename ordinal_type, typename value_type>
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type, value_type>::
PecosOneDOrthogPolyBasis(
  ordinal_type p_, const PecosOneDOrthogPolyBasis& basis) :
  pecosPoly(basis.pecosPoly),
  name(basis.name),
  p(p_),
  sparse_grid_rule(basis.sparse_grid_rule),
  sparse_grid_growth_rule(basis.sparse_grid_growth_rule),
  norms(p+1, value_type(0.0))
{
  for (ordinal_type i=0; i<=p; i++)
    norms[i] = pecosPoly->norm_squared(i);
}

template <typename ordinal_type, typename value_type>
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type, value_type>::
~PecosOneDOrthogPolyBasis()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type, value_type>::
size() const
{
  return p+1;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> >
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type, value_type>::
computeTripleProductTensor() const
{
  // Compute Cijk = < \Psi_i \Psi_j \Psi_k >
  ordinal_type sz = size();
  Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> > Cijk = 
    Teuchos::rcp(new Dense3Tensor<ordinal_type, value_type>(sz));
  Teuchos::Array<value_type> points, weights;
  Teuchos::Array< Teuchos::Array<value_type> > values;
  getQuadPoints(3*p, points, weights, values);

  for (ordinal_type i=0; i<sz; i++) {
    for (ordinal_type j=0; j<sz; j++) {
      for (ordinal_type k=0; k<sz; k++) {
	value_type triple_product = 0;
	for (ordinal_type l=0; l<static_cast<ordinal_type>(points.size());
	     l++){
	  triple_product += 
	    weights[l]*(values[l][i])*(values[l][j])*(values[l][k]);
	}
	(*Cijk)(i,j,k) = triple_product;
      }
    }
  }

  return Cijk;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type, value_type>::
computeSparseTripleProductTensor(ordinal_type order) const
{
  // Compute Cijk = < \Psi_i \Psi_j \Psi_k >
  value_type sparse_tol = 1.0e-15;
  ordinal_type sz = size();
  Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk = 
    Teuchos::rcp(new Sparse3Tensor<ordinal_type, value_type>());
  Teuchos::Array<value_type> points, weights;
  Teuchos::Array< Teuchos::Array<value_type> > values;
  getQuadPoints(3*p, points, weights, values);

  for (ordinal_type i=0; i<sz; i++) {
    for (ordinal_type j=0; j<sz; j++) {
      for (ordinal_type k=0; k<order; k++) {
	value_type triple_product = 0;
	for (ordinal_type l=0; l<static_cast<ordinal_type>(points.size());
	     l++){
	  triple_product += 
	    weights[l]*(values[l][i])*(values[l][j])*(values[l][k]);
	}
	if (std::abs(triple_product/norms[i]) > sparse_tol)
	  Cijk->add_term(i,j,k,triple_product);
      }
    }
  }
  Cijk->fillComplete();

  return Cijk;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> >
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type, value_type>::
computeDerivDoubleProductTensor() const
{
  // Compute Bij = < \Psi_i' \Psi_j >
  Teuchos::Array<value_type> points, weights;
  Teuchos::Array< Teuchos::Array<value_type> > values, derivs;
  getQuadPoints(2*p, points, weights, values);
  ordinal_type nqp = weights.size();
  derivs.resize(nqp);
  ordinal_type sz = size();
  for (ordinal_type i=0; i<nqp; i++) {
    derivs[i].resize(sz);
    evaluateBasesAndDerivatives(points[i], values[i], derivs[i]);
  }
  Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > Bij = 
    Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type>(sz,sz));
  for (ordinal_type i=0; i<sz; i++) {
    for (ordinal_type j=0; j<sz; j++) {
      value_type b = value_type(0.0);
      for (int qp=0; qp<nqp; qp++)
	b += weights[qp]*derivs[qp][i]*values[qp][j];
      (*Bij)(i,j) = b;
    }
  }

  return Bij;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type,value_type>::
evaluateBases(const value_type& x, Teuchos::Array<value_type>& basis_pts) const
{
  for (ordinal_type i=0; i<=p; i++)
    basis_pts[i] = pecosPoly->type1_value(x, i);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type,value_type>::
evaluateBasesAndDerivatives(const value_type& x, 
			    Teuchos::Array<value_type>& vals,
			    Teuchos::Array<value_type>& derivs) const
{
  for (ordinal_type i=0; i<=p; i++) {
    vals[i] = pecosPoly->type1_value(x, i);
    derivs[i] = pecosPoly->type1_gradient(x, i);
  }
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type,value_type>::
evaluate(const value_type& x, ordinal_type k) const
{
  return pecosPoly->type1_value(x, k);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type, value_type>::
print(std::ostream& os) const
{
  os << "Pecos " << name << " basis of order " << p << "." << std::endl;
  os << "Basis polynomial norms (squared):\n\t";
  for (ordinal_type i=0; i<=p; i++)
    os << norms[i] << " ";
  os << std::endl;
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type,value_type>::
getQuadPoints(ordinal_type quad_order,
	      Teuchos::Array<value_type>& quad_points,
	      Teuchos::Array<value_type>& quad_weights,
	      Teuchos::Array< Teuchos::Array<value_type> >& quad_values) const
{
  ordinal_type num_points = 
    static_cast<ordinal_type>(std::ceil((quad_order+1)/2.0));
  const Pecos::RealArray& gp = pecosPoly->collocation_points(num_points);
  const Pecos::RealArray& gw = pecosPoly->type1_collocation_weights(num_points); 
  quad_points.resize(num_points);
  quad_weights.resize(num_points);
  for (ordinal_type i=0; i<num_points; i++) {
    quad_points[i] = gp[i];
    quad_weights[i] = gw[i];
  }
  
  // Evalute basis at gauss points
  quad_values.resize(num_points);
  for (ordinal_type i=0; i<num_points; i++) {
    quad_values[i].resize(p+1);
    evaluateBases(quad_points[i], quad_values[i]);
  }
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > 
Stokhos::PecosOneDOrthogPolyBasis<ordinal_type,value_type>::
cloneWithOrder(ordinal_type p) const
{
  return 
    Teuchos::rcp(new Stokhos::PecosOneDOrthogPolyBasis<ordinal_type,value_type>(p,*this));
}
