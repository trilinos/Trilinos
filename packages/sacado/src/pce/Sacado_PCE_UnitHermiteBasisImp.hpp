// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

template <typename T>
Sacado::PCE::UnitHermiteBasis<T>::
UnitHermiteBasis(unsigned int degree) :
  d(degree),
  basis(d+1, StandardPoly<T>(d)),
  norms(d+1, T(0.0))
{
  // Fill in basis coefficients
  basis[0].coeff(0) = T(1.0);
  if (d >= 1)
    basis[1].coeff(1) = T(2.0);
  for (unsigned int k=2; k<=d; k++) {
    basis[k].coeff(0) = -T(2.0)*T(k-1)*basis[k-2].coeff(0);
    for (unsigned int i=1; i<=k; i++)
      basis[k].coeff(i) = 
	T(2.0)*(basis[k-1].coeff(i-1) - T(k-1)*basis[k-2].coeff(i));
  }

  // Compute norms
  norms[0] = 2.0*std::sqrt(std::atan(1.0));  // = sqrt(pi)
  for (unsigned int k=1; k<=d; k++)
    norms[k] = T(2.0)*T(k)*norms[k-1];

  // Rescale to unit norm
  for (unsigned int k=0; k<=d; k++) {
    for (unsigned int i=0; i<=d; i++)
      basis[k].coeff(i) /= std::sqrt(norms[k]);
    norms[k] = 1.0;
  }
}

template <typename T>
Sacado::PCE::UnitHermiteBasis<T>::
UnitHermiteBasis(const Sacado::PCE::UnitHermiteBasis<T>& b) :
  d(b.d),
  basis(b.basis),
  norms(b.norms)
{
}

template <typename T>
Sacado::PCE::UnitHermiteBasis<T>& 
Sacado::PCE::UnitHermiteBasis<T>::
operator=(const Sacado::PCE::UnitHermiteBasis<T>& b)
{
  if (this != &b) {
    d = b.d;
    basis = b.basis;
    norms = b.norms;
  }
  return *this;
}

template <typename T>
Sacado::PCE::UnitHermiteBasis<T>::
~UnitHermiteBasis()
{
}

template <typename T>
unsigned int
Sacado::PCE::UnitHermiteBasis<T>::
size() const
{
  return d+1;
}

template <typename T>
const std::vector<T>&
Sacado::PCE::UnitHermiteBasis<T>::
norm_squared() const
{
  return norms;
}

template <typename T>
T
Sacado::PCE::UnitHermiteBasis<T>::
derivCoeff(unsigned int i) const
{
  return std::sqrt(2.0*i);
}

template <typename T>
void
Sacado::PCE::UnitHermiteBasis<T>::
project(const Sacado::PCE::StandardPoly<T>& p, std::vector<T>& coeffs) const
{
  // Initialize
  for (unsigned int i=0; i<=d; i++)
    coeffs[i] = T(0.0);

  unsigned int dp = p.degree();
  if (dp > d)
    dp = d;

  // Handle degree 0 case
  if (dp == 0) {
    coeffs[0] = p.coeff(0);
    return;
  }

  // Temporary array
  std::vector<T> h(dp+1,T(0.));

  T pi = 4.0*std::atan(1.0);
  T piq = std::sqrt(std::sqrt(pi));
  T sq2 = std::sqrt(2.0);

  coeffs[0] = piq*p.coeff(dp-1);
  coeffs[1] = piq/sq2*p.coeff(dp);
  unsigned int dc = 1;

  for (int k=dp-2; k>=0; --k) {

    // Multiply by t
    h[0] = coeffs[1]/sq2 + piq*p.coeff(k);
    for (unsigned int i=1; i<=dc-1; i++) {
      h[i] = (std::sqrt(T(i))*coeffs[i-1] + std::sqrt(T(i+1))*coeffs[i+1])/sq2;
    }
    h[dc] = std::sqrt(T(dc))/sq2*coeffs[dc-1];
    h[dc+1] = std::sqrt(T(dc+1))/sq2*coeffs[dc];
    
    // Copy into coeffs
    for (unsigned int i=0; i<=dc+1; i++)
      coeffs[i] = h[i];

    dc = dc+1;
  }
}

template <typename T>
Sacado::PCE::StandardPoly<T>
Sacado::PCE::UnitHermiteBasis<T>::
toStandardBasis(const T coeffs[], unsigned int n) const
{
  unsigned int dp = d;
  if (n < d+1)
    dp = n-1;

  StandardPoly<T> p(dp);

  for (unsigned int i=0; i<=dp; i++)
    p.add(coeffs[i], basis[i], 1.0);

  return p;
}

template <typename T>
const Sacado::PCE::StandardPoly<T>&
Sacado::PCE::UnitHermiteBasis<T>::
getBasisPoly(unsigned int i) const
{
  return basis[i];
}

template <typename T>
void
Sacado::PCE::UnitHermiteBasis<T>::
print(std::ostream& os) const
{
  os << "Hermite basis of degree " << d << ".  Basis vectors:\n";
  for (unsigned int i=0; i<basis.size(); i++)
    os << "\t" << basis[i];
  os << "Basis vector norms (squared):\n\t";
  for (unsigned int i=0; i<norms.size(); i++)
    os << norms[i] << " ";
  os << "\n";
}
