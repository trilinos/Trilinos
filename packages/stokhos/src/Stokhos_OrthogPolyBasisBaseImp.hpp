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

template <typename T>
Stokhos::OrthogPolyBasisBase<T>::
OrthogPolyBasisBase(const std::string& name_, unsigned int p_) :
  name(name_),
  p(p_),
  basis(p+1, Polynomial<T>(p)),
  double_basis(2*p+1, Polynomial<T>(2*p)),
  norms(p+1, T(0.0))
{
}

template <typename T>
Stokhos::OrthogPolyBasisBase<T>::
~OrthogPolyBasisBase()
{
}

template <typename T>
unsigned int
Stokhos::OrthogPolyBasisBase<T>::
order() const
{
  return p;
}

template <typename T>
unsigned int
Stokhos::OrthogPolyBasisBase<T>::
dimension() const
{
  return 1;
}

template <typename T>
unsigned int
Stokhos::OrthogPolyBasisBase<T>::
size() const
{
  return p+1;
}

template <typename T>
const std::vector<T>&
Stokhos::OrthogPolyBasisBase<T>::
norm_squared() const
{
  return norms;
}

template <typename T>
void
Stokhos::OrthogPolyBasisBase<T>::
projectProduct(unsigned int i, unsigned int j, std::vector<T>& coeffs) const
{
  // Important note:  To get the correct value for <\Psi_i \Psi_j \Psi_k>,
  // we have to expand \Psi_i*\Psi_j in the full order_i+order_j basis, not
  // just the order p basis.  There for the basis needs to be of size 2*p
  Polynomial<T> pij(2*p);
  std::vector<T> a(2*(p+1));

  // Multiply ith and jth basis polynomial
  pij.multiply(1.0, basis[i], basis[j], 0.0);

  // Project onto basis
  projectPoly(pij, a);
  for (unsigned int k=0; k<=p; k++)
    coeffs[k] = a[k]*norms[k];
}

template <typename T>
Stokhos::Polynomial<T>
Stokhos::OrthogPolyBasisBase<T>::
toStandardBasis(const T coeffs[], unsigned int n) const
{
  unsigned int px = p;
  if (n < p+1)
    px = n-1;
  Polynomial<T> x(px);

  for (unsigned int i=0; i<=px; i++)
    x.add(coeffs[i], basis[i], 1.0);

  return x;
}

template <typename T>
T
Stokhos::OrthogPolyBasisBase<T>::
evaluateZero(unsigned int i) const
{
  return basis[i].coeff(0);
}

template <typename T>
void
Stokhos::OrthogPolyBasisBase<T>::
print(std::ostream& os) const
{
  os << name << " basis of order " << p << ".  Basis vectors:\n";
  for (unsigned int i=0; i<basis.size(); i++)
    os << "\t" << basis[i];
  os << "Basis vector norms (squared):\n\t";
  for (unsigned int i=0; i<norms.size(); i++)
    os << norms[i] << " ";
  os << "\n";
}

template <typename T>
std::vector<unsigned int>
Stokhos::OrthogPolyBasisBase<T>::
getTerm(unsigned int i) const
{
  std::vector<unsigned int> t(1);
  t[0] = i;
  return t;
}

template <typename T>
unsigned int
Stokhos::OrthogPolyBasisBase<T>::
getIndex(const std::vector<unsigned int>& term) const
{
  return term[0];
}

template <typename T>
const std::string&
Stokhos::OrthogPolyBasisBase<T>::
getName() const
{
  return name;
}
