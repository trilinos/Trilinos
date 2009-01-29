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

#include "Teuchos_TestForException.hpp"

template <typename T>
Stokhos::CompletePolynomialBasis<T>::
CompletePolynomialBasis(
	const std::vector< Teuchos::RCP<const OrthogPolyBasis<T> > >& bases_,
	const std::vector<T>& deriv_coeffs_) :
  p(0),
  d(0),
  sz(0),
  bases(bases_),
  deriv_coeffs(deriv_coeffs_),
  norms(),
  terms()
{
  // Compute total dimension -- we are assuming each basis is 1-D
  d = bases.size();

  // Compute total order
  for (unsigned int i=0; i<bases.size(); i++)
    if (bases[i]->order() > p)
      p = bases[i]->order();

  // Compute basis terms
  compute_terms();

  // Create triple products
  Cijk.resize(bases.size());
  for (unsigned int i=0; i<bases.size(); i++)
    Cijk[i] = Teuchos::rcp(new TripleProduct< OrthogPolyBasis<T> >(bases[i]));
  
  // Compute norms
  norms.resize(sz);
  T nrm;
  for (unsigned int k=0; k<sz; k++) {
    nrm = T(1.0);
    for (unsigned int i=0; i<bases.size(); i++)
      nrm = nrm * Cijk[i]->norm_squared(terms[k][i]);
    norms[k] = nrm;
  }

  // Create name
  name = "Complete polynomial basis (";
  for (unsigned int i=0; i<bases.size()-1; i++)
    name += bases[i]->getName() + ", ";
  name += bases[bases.size()-1]->getName() + ")";
}

template <typename T>
Stokhos::CompletePolynomialBasis<T>::
~CompletePolynomialBasis()
{
}

template <typename T>
unsigned int
Stokhos::CompletePolynomialBasis<T>::
order() const
{
  return p;
}

template <typename T>
unsigned int
Stokhos::CompletePolynomialBasis<T>::
dimension() const
{
  return d;
}

template <typename T>
unsigned int
Stokhos::CompletePolynomialBasis<T>::
size() const
{
  return sz;
}

template <typename T>
const std::vector<T>&
Stokhos::CompletePolynomialBasis<T>::
norm_squared() const
{
  return norms;
}

template <typename T>
void
Stokhos::CompletePolynomialBasis<T>::
projectPoly(const Polynomial<T>& poly, std::vector<T>& coeffs) const
{
  TEST_FOR_EXCEPTION(true, std::logic_error, "Stokhos::CompletePolynomialBasis<T>::projectPoly is not implemented! (for good reason, it is very expensive!)");
}

template <typename T>
void
Stokhos::CompletePolynomialBasis<T>::
projectProduct(unsigned int i, unsigned int j, std::vector<T>& coeffs) const
{
  for (unsigned int k=0; k<sz; k++) {
    T c = T(1.0);
    for (unsigned int l=0; l<bases.size(); l++) {
      c *= Cijk[l]->triple_value(terms[i][l],terms[j][l],terms[k][l]);
      //Cijk[l]->norm_squared(terms[k][l]);
    }
    coeffs[k] = c;
  }
}

// template <typename T>
// void
// Stokhos::CompletePolynomialBasis<T>::
// projectDerivative(unsigned int i, std::vector<T>& coeffs) const
// {
//   // Initialize
//   for (unsigned int j=0; j<coeffs.size(); j++)
//     coeffs[j] = T(0.0);

//   // Project derivative of each basis polynomial for term i
//   std::vector<T> bases_coeffs(p+1);
//   std::vector<unsigned int> term(d);
//   unsigned int index;
//   for (unsigned int j=0; j<d; j++) {
//     bases[j]->projectDerivative(terms[i][j],bases_coeffs);

//     term = terms[i];
//     for (unsigned int k=0; k<terms[i][j]; k++) {
//       term[j] = k;
//       index = compute_index(term);
//       coeffs[index] += deriv_coeffs[j]*bases_coeffs[k];
//     }
//   }
// }

template <typename T>
void
Stokhos::CompletePolynomialBasis<T>::
projectDerivative(unsigned int i, std::vector<T>& coeffs) const
{
  // Initialize
  for (unsigned int j=0; j<coeffs.size(); j++)
    coeffs[j] = T(0.0);

  for (unsigned int k=0; k<sz; k++) {
    T t = T(1.0);
    for (unsigned int j=0; j<d; j++) {
      bool is_zero = false;
      for (unsigned int l=0; l<d; l++) {
	if (l != j && terms[i][l] != terms[k][l])
	  is_zero = true;
	if (l != j)
	  t *= Cijk[l]->norm_squared(terms[k][l]);
      }
      if (!is_zero)
	coeffs[k] += 
	  t*deriv_coeffs[j]*Cijk[j]->double_deriv(terms[k][j],terms[i][j]);
    }
    coeffs[k] /= norms[k];
  }
}

template <typename T>
Stokhos::Polynomial<T>
Stokhos::CompletePolynomialBasis<T>::
toStandardBasis(const T coeffs[], unsigned int n) const
{
  TEST_FOR_EXCEPTION(true, std::logic_error, "Stokhos::CompletePolynomialBasis<T>::toStandardBasis is not implemented! (for good reason, it is very expensive!)");

  return Polynomial<T>(sz);
}

template <typename T>
T
Stokhos::CompletePolynomialBasis<T>::
evaluateZero(unsigned int i) const
{
  // z = psi_{i_1}(0) * ... * psi_{i_d}(0) where i_1,...,i_d are the basis
  // terms for coefficient i
  T z = T(1.0);
  for (unsigned int j=0; j<bases.size(); j++)
    z = z * bases[j]->evaluateZero(terms[i][j]);

  return z;
}

template <typename T>
void
Stokhos::CompletePolynomialBasis<T>::
evaluateBases(const std::vector<T>& point, std::vector<T>& basis_pts) const
{
  std::vector< std::vector<T> > tmp(bases.size());
  for (unsigned int j=0; j<bases.size(); j++) {
    std::vector<T> pt(1);
    pt[0] = point[j];
    tmp[j].resize(bases[j]->order()+1);
    bases[j]->evaluateBases(pt, tmp[j]);
  }

  for (unsigned int i=0; i<sz; i++) {
    T t = T(1.0);
    for (unsigned int j=0; j<bases.size(); j++)
      t *= tmp[j][terms[i][j]];
    basis_pts[i] = t;
  }
}

template <typename T>
void
Stokhos::CompletePolynomialBasis<T>::
print(std::ostream& os) const
{
  os << "Complete basis of order " << p << ", dimension " << d 
     << ", and size " << sz << ".  Component bases:\n";
  for (unsigned int i=0; i<bases.size(); i++)
    os << *bases[i];
  os << "Basis vector norms (squared):\n\t";
  for (unsigned int i=0; i<norms.size(); i++)
    os << norms[i] << " ";
  os << "\n";
}

template <typename T>
std::vector<unsigned int>
Stokhos::CompletePolynomialBasis<T>::
getTerm(unsigned int i) const
{
  return terms[i];
}

template <typename T>
unsigned int
Stokhos::CompletePolynomialBasis<T>::
getIndex(const std::vector<unsigned int>& term) const
{
  return compute_index(term);
}

template <typename T>
const std::string&
Stokhos::CompletePolynomialBasis<T>::
getName() const
{
  return name;
}

template <typename T>
unsigned int
Stokhos::CompletePolynomialBasis<T>::
compute_num_terms(unsigned int dim, unsigned int ord) const
{
  unsigned int num = 1;
  
  // Use the formula (p+d)!/(p!d!) = (d+p)...(d+1)/p!
  if (dim >= ord) {
    for (unsigned int i=1; i<=ord; i++) {
      num *= dim+i;
      num /= i;
    }
  }

  // Use the formula (p+d)!/(p!d!) = (p+d)...(p+1)/d!
  else {
    for (unsigned int i=1; i<=dim; i++) {
      num *= ord+i;
      num /= i;
    }
  }

  return num;
}

template <typename T>
void
Stokhos::CompletePolynomialBasis<T>::
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
  for (unsigned int i=0; i<sz; i++) {
    terms[i].resize(d);
    for (unsigned int j=0; j<d; j++)
      terms[i][j] = 0;
  }

  if (p == 0)
    return;

  // The array "cnt" stores the number of terms we need to increment for each
  // dimension.  
  std::vector<unsigned int> cnt(d);

  // Set order 1 terms
  for (unsigned int j=0; j<d; j++) {
    terms[j+1][j] = 1;
    cnt[j] = d-j;
  }

  // Stores index of previous order block
  unsigned int prev = 1;

  // Stores index of the term we are working on
  unsigned int cur = d+1;

  // Loop over orders
  for (unsigned int k=2; k<=p; k++) {

    // Loop over dimensions
    for (unsigned int j=0; j<d; j++) {

      // Increment orders of cnt[j] terms for dimension j
      for (unsigned int i=0; i<cnt[j]; i++) {
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
    for (unsigned int j=0; j<d; j++)
      for (unsigned int i=j+1; i<d; i++)
	cnt[j] += cnt[i];

  }
  
}

template <typename T>
unsigned int 
Stokhos::CompletePolynomialBasis<T>::
compute_index(const std::vector<unsigned int>& term) const
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

  unsigned int index = 0;
  int dim = term.size();
  unsigned int ord = 0;
  for (int i=dim-1; i>=0; i--) {

    // compute order
    ord += term[i];

    // compute number of terms for order-1
    if (ord > 0)
      index += compute_num_terms(dim-i, ord-1);
  }

  return index;
}
