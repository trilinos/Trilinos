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

template <typename ordinal_type, typename value_type>
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
CompletePolynomialBasis(
	const Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, value_type> > >& bases_,
	const value_type& sparse_tol_,
	const Teuchos::RCP< Teuchos::Array<value_type> >& deriv_coeffs_) :
  p(0),
  d(0),
  sz(0),
  bases(bases_),
  sparse_tol(sparse_tol_),
  deriv_coeffs(deriv_coeffs_),
  norms(),
  terms()
{
  // Compute total dimension -- we are assuming each basis is 1-D
  d = bases.size();

  // Compute total order
  for (ordinal_type i=0; i<static_cast<ordinal_type>(bases.size()); i++)
    if (bases[i]->order() > p)
      p = bases[i]->order();

  // Compute basis terms
  compute_terms();
    
  // Compute norms
  norms.resize(sz);
  value_type nrm;
  for (ordinal_type k=0; k<sz; k++) {
    nrm = value_type(1.0);
    for (ordinal_type i=0; i<static_cast<ordinal_type>(bases.size()); i++)
      nrm = nrm * bases[i]->norm_squared(terms[k][i]);
    norms[k] = nrm;
  }

  // Create name
  name = "Complete polynomial basis (";
  for (ordinal_type i=0; i<static_cast<ordinal_type>(bases.size())-1; i++)
    name += bases[i]->getName() + ", ";
  name += bases[bases.size()-1]->getName() + ")";

  // Allocate array for basis evaluation
  basis_eval_tmp.resize(bases.size());
  for (ordinal_type j=0; j<static_cast<ordinal_type>(bases.size()); j++)
    basis_eval_tmp[j].resize(bases[j]->order()+1);

  // Set up deriv_coeffs
  if (deriv_coeffs == Teuchos::null) {
    deriv_coeffs = Teuchos::rcp(new Teuchos::Array<value_type>(bases.size()));
    for (ordinal_type j=0; j<static_cast<ordinal_type>(bases.size()); j++)
      (*deriv_coeffs)[j] = value_type(1.0);
  }
}

template <typename ordinal_type, typename value_type>
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
~CompletePolynomialBasis()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
dimension() const
{
  return d;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
size() const
{
  return sz;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
getTripleProductTensor() const
{/*
  ordinal_type sz = size();
  
  // Compute Cijk = < \Psi_i \Psi_j \Psi_k >
  if (Cijk == Teuchos::null) {
    Teuchos::Array<value_type> points, weights;
    Teuchos::Array< Teuchos::Array<value_type> > values;
    getQuadPoints(ceil((3*p+1)/2), points, weights, values);
    Cijk = Teuchos::rcp(new Sparse3Tensor<ordinal_type, value_type>(sz));
    
    
    for (ordinal_type i=0; i<sz; i++) {
      for (ordinal_type j=0; j<sz; j++) {
	for (ordinal_type k=0; k<sz; k++) {
          value_type triple_product = 0;
	  for (ordinal_type l=0; l<points.size();l++){
             triple_product = triple_product + weights[l]*(values[l][i])*(values[l][j])*(values[l][k]);
             //if(k == 0) std::cout<< "values[0]["<<l<<"] = "<<values[l][k] <<"\n";
          }
          if(std::abs(triple_product)>sparse_tol){
            Cijk->add_term(i,j,k,triple_product);
          //if(i == 0 && j == 0 && k == 0) std::cout<< "C000 = " << (*Cijk)(i,j,k) << "\n";
          }
	}
      }
    }
  }
  
  return Cijk;

  */
  
  // Compute Cijk = < \Psi_i \Psi_j \Psi_k >
  if (Cijk == Teuchos::null) {
    Cijk = Teuchos::rcp(new Sparse3Tensor<ordinal_type, value_type>(sz));
    Teuchos::Array<value_type> a(sz);
    for (ordinal_type j=0; j<sz; j++) {
      for (ordinal_type i=0; i<sz; i++) {
	projectProduct(i, j, a);
	for (ordinal_type k=0; k<sz; k++) {
          //std::cout << "k = " << k << " " << a[k] << "\n";
	  if (std::abs(a[k]) > sparse_tol) {
	    Cijk->add_term(i,j,k,a[k]);
	  }
	}
      }
    }
  }

  return Cijk;
  
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
getLowOrderTripleProductTensor(ordinal_type order) const
{
  // Compute Cijk = < \Psi_i \Psi_j \Psi_k >
  if (Cijk == Teuchos::null) {
    Cijk = Teuchos::rcp(new Sparse3Tensor<ordinal_type, value_type>(sz));
    Teuchos::Array<value_type> a(order);
    for (ordinal_type i=0; i<sz; i++) {
      for (ordinal_type j=0; j<sz; j++) {
	projectProduct(i, j, a);
	for (ordinal_type k=0; k<order; k++) {
	  if (std::abs(a[k]) > sparse_tol) {
	    Cijk->add_term(i,j,k,a[k]);
	  }
	}
      }
    }
  }

  return Cijk;
  
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< const Stokhos::Dense3Tensor<ordinal_type, value_type> >
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
getDerivTripleProductTensor() const
{
  // Compute Dijk = < \Psi_i \Psi_j \Psi_k' >
  if (Dijk == Teuchos::null) {

    // Compute Cijk & Bij
    getTripleProductTensor();
    getDerivDoubleProductTensor();

    Dijk = Teuchos::rcp(new Dense3Tensor<ordinal_type, value_type>(sz));
    for (ordinal_type i=0; i<sz; i++)
      for (ordinal_type j=0; j<sz; j++)
	for (ordinal_type k=0; k<sz; k++)
	  (*Dijk)(i,j,k) = value_type(0.0);

    ordinal_type i,j;
    value_type c;
    for (ordinal_type k=0; k<sz; k++) {
      for (ordinal_type m=0; m<sz; m++) {
	ordinal_type n = Cijk->num_values(m);
	for (ordinal_type l=0; l<n; l++) {
	  Cijk->value(m,l,i,j,c);
	  (*Dijk)(i,j,k) += (*Bij)(m,k)*c/norms[m];
	}
      }
    }
  }

  return Dijk;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< const Teuchos::SerialDenseMatrix<ordinal_type, value_type> >
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
getDerivDoubleProductTensor() const
{
  // Compute Bij = < \Psi_i \Psi_j' >
  if (Bij == Teuchos::null) {
    Bij = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type>(sz,sz));
    for (ordinal_type i=0; i<sz; i++) {
      Teuchos::Array<value_type> b(sz);
      projectDerivative(i, b);
      for (ordinal_type j=0; j<sz; j++)
	(*Bij)(i,j) = b[j]*norms[j];
    }
  }

  return Bij;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
projectProduct(ordinal_type i, ordinal_type j, Teuchos::Array<value_type>& coeffs) const
{
  // Create products
  if (Cijk_1d.size() == 0) {
    Cijk_1d.resize(bases.size());
    for (ordinal_type i=0; i<static_cast<ordinal_type>(bases.size()); i++)
      Cijk_1d[i] = bases[i]->getTripleProductTensor();
  }

  // Commenting out this modification for now.  There seem to be cases where
  // truncating at sparse_tol gives incorrect answers

  // double val;
  // for (ordinal_type k=0; k<coeffs.size(); k++) {
  //   value_type c = value_type(1.0);
  //   for (ordinal_type l=0; l<static_cast<ordinal_type>(bases.size()); l++) {
  //     val = (*Cijk_1d[l])(terms[i][l],terms[j][l],terms[k][l]);
  //     if(val < sparse_tol){
  //       c = 0;
  //       break;
  //     }else c *= val;      
  //   }
  //   coeffs[k] = c;
  // }

  double val;
  for (ordinal_type k=0; k<static_cast<ordinal_type>(coeffs.size()); k++) {
    value_type c = value_type(1.0);
    for (ordinal_type l=0; l<static_cast<ordinal_type>(bases.size()); l++) {
      val = (*Cijk_1d[l])(terms[i][l],terms[j][l],terms[k][l]);
      c *= val;      
    }
    coeffs[k] = c;
  }
}


template <typename ordinal_type, typename value_type>
void
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
projectDerivative(ordinal_type i, Teuchos::Array<value_type>& coeffs) const
{
  // Create products
  if (Bij_1d.size() == 0) {
    Bij_1d.resize(bases.size());
    for (ordinal_type i=0; i<static_cast<ordinal_type>(bases.size()); i++)
      Bij_1d[i] = bases[i]->getDerivDoubleProductTensor();
  }

  // Initialize
  for (ordinal_type j=0; j<static_cast<ordinal_type>(coeffs.size()); j++)
    coeffs[j] = value_type(0.0);

  for (ordinal_type k=0; k<sz; k++) {
    value_type t = value_type(1.0);
    for (ordinal_type j=0; j<d; j++) {
      bool is_zero = false;
      for (ordinal_type l=0; l<d; l++) {
	if (l != j && terms[i][l] != terms[k][l])
	  is_zero = true;
	if (l != j)
	  t *= bases[l]->norm_squared(terms[k][l]);
      }
      if (!is_zero)
        coeffs[k] += 
          t*(*deriv_coeffs)[j]*(*Bij_1d[j])(terms[k][j],terms[i][j]);
    }
    coeffs[k] /= norms[k];
  }
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
evaluateZero(ordinal_type i) const
{
  // z = psi_{i_1}(0) * ... * psi_{i_d}(0) where i_1,...,i_d are the basis
  // terms for coefficient i
  value_type z = value_type(1.0);
  for (ordinal_type j=0; j<static_cast<ordinal_type>(bases.size()); j++)
    z = z * bases[j]->evaluateZero(terms[i][j]);

  return z;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
evaluateBases(const Teuchos::Array<value_type>& point,
	      Teuchos::Array<value_type>& basis_vals) const
{
  for (ordinal_type j=0; j<static_cast<ordinal_type>(bases.size()); j++)
    bases[j]->evaluateBases(point[j], basis_eval_tmp[j]);

  // Only evaluate basis upto number of terms included in basis_pts
  for (ordinal_type i=0; i<sz; i++) {
    value_type t = value_type(1.0);
    for (ordinal_type j=0; j<static_cast<ordinal_type>(bases.size()); j++)
      t *= basis_eval_tmp[j][terms[i][j]];
    basis_vals[i] = t;
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
print(std::ostream& os) const
{
  os << "Complete basis of order " << p << ", dimension " << d 
     << ", and size " << sz << ".  Component bases:\n";
  for (ordinal_type i=0; i<static_cast<ordinal_type>(bases.size()); i++)
    os << *bases[i];
  os << "Basis vector norms (squared):\n\t";
  for (ordinal_type i=0; i<static_cast<ordinal_type>(norms.size()); i++)
    os << norms[i] << " ";
  os << "\n";
}

template <typename ordinal_type, typename value_type>
Teuchos::Array<ordinal_type>
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
getTerm(ordinal_type i) const
{
  return terms[i];
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
getIndex(const Teuchos::Array<ordinal_type>& term) const
{
  return compute_index(term);
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> > >&
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
getCoordinateBases() const
{
  return bases;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
compute_num_terms(ordinal_type dim, ordinal_type ord) const
{
  ordinal_type num = 1;
  
  // Use the formula (p+d)!/(p!d!) = (d+p)...(d+1)/p!
  if (dim >= ord) {
    for (ordinal_type i=1; i<=ord; i++) {
      num *= dim+i;
      num /= i;
    }
  }

  // Use the formula (p+d)!/(p!d!) = (p+d)...(p+1)/d!
  else {
    for (ordinal_type i=1; i<=dim; i++) {
      num *= ord+i;
      num /= i;
    }
  }

  return num;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
compute_terms()
{
  // The approach here for ordering the terms is inductive on the total
  // order p.  We get the terms of total order p from the terms of total
  // order p-1 by incrementing the orders of the first dimension by 1.
  // We then increment the orders of the second dimension by 1 for all of the
  // terms whose first dimension order is 0.  We then repeat for the third
  // dimension whose first and second dimension orders are 0, and so on.
  // How this is done is most easily illustrated by an example of dimension 3:
  //
  // Order  terms   cnt  Order  terms   cnt
  //   0    0 0 0          4    4 0 0  15 5 1
  //                            3 1 0
  //   1    1 0 0  3 2 1        3 0 1
  //        0 1 0               2 2 0
  //        0 0 1               2 1 1
  //                            2 0 2
  //   2    2 0 0  6 3 1        1 3 0
  //        1 1 0               1 2 1
  //        1 0 1               1 1 2
  //        0 2 0               1 0 3
  //        0 1 1               0 4 0
  //        0 0 2               0 3 1
  //                            0 2 2
  //   3    3 0 0  10 4 1       0 1 3
  //        2 1 0               0 0 4
  //        2 0 1
  //        1 2 0
  //        1 1 1
  //        1 0 2
  //        0 3 0
  //        0 2 1
  //        0 1 2
  //        0 0 3

  // First compute total size from (d+p)!/(d!p!)
  sz = compute_num_terms(d,p);

  // Allocate storage and initialize
  terms.resize(sz);
  for (ordinal_type i=0; i<sz; i++) {
    terms[i].resize(d);
    for (ordinal_type j=0; j<d; j++)
      terms[i][j] = 0;
  }

  if (p == 0)
    return;

  // The array "cnt" stores the number of terms we need to increment for each
  // dimension.  
  Teuchos::Array<ordinal_type> cnt(d);

  // Set order 1 terms
  for (ordinal_type j=0; j<d; j++) {
    terms[j+1][j] = 1;
    cnt[j] = d-j;
  }

  // Stores index of previous order block
  ordinal_type prev = 1;

  // Stores index of the term we are working on
  ordinal_type cur = d+1;

  // Loop over orders
  for (ordinal_type k=2; k<=p; k++) {

    // Loop over dimensions
    for (ordinal_type j=0; j<d; j++) {

      // Increment orders of cnt[j] terms for dimension j
      for (ordinal_type i=0; i<cnt[j]; i++) {
        terms[cur] = terms[prev+i];
        ++terms[cur][j];
        ++cur;
      }

      // Move forward the index of the previous order block.  If we aren't
      // at the last dimension, the amount we move forward is
      // cnt[j]-cnt[j+1].  If we are at the last dimension, we just increment
      // by 1
      if (j < d-1)
        prev += cnt[j]-cnt[j+1];
      else
        ++prev;
    }

    // Compute the number of terms we must increment for the new order
    // For each dimension j, this is just number plus the sum of the number 
    // of terms for the remaining d-j dimensions
    for (ordinal_type j=0; j<d; j++)
      for (ordinal_type i=j+1; i<d; i++)
        cnt[j] += cnt[i];

  }
  
}

template <typename ordinal_type, typename value_type>
ordinal_type 
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
compute_index(const Teuchos::Array<ordinal_type>& term) const
{
  // The approach here for computing the index is essentially recursive
  // on the number of dimensions.  Given the basis orders for each dimenion
  // in "term", we add the orders to get the total order, and then compute
  // the number of terms in an order-1 expansion.  That tells us which
  // order block "term" lies in the global "terms" array.  We then compute
  // the total order in the last d-1 dimensions and the number of terms in
  // an order-1, d-1 expansion, adding this to the previous offset.  We
  // repeat this until there are no dimensions left, which provides the index.
  //
  // For efficiency, we actually work from the last dimension to the first
  // to reduce the number of operations to compute the total order.

  ordinal_type index = 0;
  int dim = term.size();
  ordinal_type ord = 0;
  for (int i=dim-1; i>=0; i--) {

    // compute order
    ord += term[i];

    // compute number of terms for order-1
    if (ord > 0)
      index += compute_num_terms(dim-i, ord-1);
  }

  return index;
}
