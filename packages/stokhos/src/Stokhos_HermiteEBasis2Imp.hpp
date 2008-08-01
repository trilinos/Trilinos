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
Stokhos::HermiteEBasis2<T>::
HermiteEBasis2(unsigned int p) :
  OrthogPolyBasisBase<T>("HermiteE2",p)
{
  // Fill in basis coefficients
  this->basis[0].coeff(0) = T(1.0);
  if (this->p >= 1)
    this->basis[1].coeff(1) = T(1.0);
  for (unsigned int k=2; k<=this->p; k++) {
    this->basis[k].coeff(0) = -T(k-1)*(this->basis[k-2].coeff(0));
    for (unsigned int i=1; i<=k; i++)
      this->basis[k].coeff(i) = 
	this->basis[k-1].coeff(i-1) - T(k-1)*(this->basis[k-2].coeff(i));
  }

  // Compute norms
  this->norms[0] = T(1.0);
  for (unsigned int k=1; k<=this->p; k++)
    this->norms[k] = T(k)*(this->norms[k-1]);
}

template <typename T>
Stokhos::HermiteEBasis2<T>::
~HermiteEBasis2()
{
}

template <typename T>
void
Stokhos::HermiteEBasis2<T>::
projectPoly(const Stokhos::Polynomial<T>& x, std::vector<T>& coeffs) const
{
  // Initialize
  for (unsigned int i=0; i<=this->p; i++)
    coeffs[i] = T(0.0);

  unsigned int px = x.degree();
  if (px > this->p && px != this->p*2)
    px = this->p;

  // Handle degree 0 case
  if (px == 0) {
    coeffs[0] = x.coeff(0);
    return;
  }

  // Temporary array
  std::vector<T> h(px+1,T(0.));

  coeffs[0] = x.coeff(px-1);
  coeffs[1] = x.coeff(px);
  unsigned int pc = 1;

  for (int k=px-2; k>=0; --k) {

    // Multiply by t
    h[0] = coeffs[1] + x.coeff(k);
    for (unsigned int i=1; i<=pc-1; i++) {
      h[i] = coeffs[i-1] + T(i+1)*coeffs[i+1];
    }
    h[pc] = coeffs[pc-1];
    h[pc+1] = coeffs[pc];
    
    // Copy into coeffs
    for (unsigned int i=0; i<=pc+1; i++)
      coeffs[i] = h[i];

    pc = pc+1;
  }
}

template <typename T>
void
Stokhos::HermiteEBasis2<T>::
projectDerivative(unsigned int i, std::vector<T>& coeffs) const
{
  for (unsigned int j=0; j<coeffs.size(); j++)
    coeffs[j] = T(0.0);
  if (i > 0)
    coeffs[i-1] = T(i);
}

template <typename T>
void
Stokhos::HermiteEBasis2<T>::
evaluateBases(const std::vector<T>& point, std::vector<T>& basis_pts) const
{
  const T& x = point[0];

  // Evaluate basis polynomials He(x) using 3 term recurrence
  // He_0(x) = 1
  // He_1(x) = x
  // He_i(x) = x*He_{i-1}(x) - (i-1)*He_{i-2}(x), i=2,3,...
  basis_pts[0] = T(1.0);
  if (this->p >= 1)
    basis_pts[1] = x;
  for (unsigned int i=2; i<=this->p; i++)
    basis_pts[i] = x*basis_pts[i-1] - T(i-1)*basis_pts[i-2];
}
