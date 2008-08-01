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
Stokhos::LegendreBasis<T>::
LegendreBasis(unsigned int p) :
  OrthogPolyBasisBase<T>("Legendre",p),
  deriv_coeffs(p+1)
{
  // Fill in basis coefficients
  T c1, c2;
  this->basis[0].coeff(0) = T(1.0);
  if (this->p >= 1)
    this->basis[1].coeff(1) = T(1.0);
  for (unsigned int k=2; k<=this->p; k++) {
    c1 = T( (2.0*k-1.0) / k );
    c2 = T( (k-1.0) / k );
    this->basis[k].coeff(0) = -c2*(this->basis[k-2].coeff(0));
    for (unsigned int i=1; i<=k; i++)
      this->basis[k].coeff(i) = 
	c1*(this->basis[k-1].coeff(i-1)) - c2*(this->basis[k-2].coeff(i));
  }

  // Compute norms
  for (unsigned int k=0; k<=this->p; k++)
    this->norms[k] = T(1.0/(2.0*k+1.0));

  // Compute derivative coefficients
  deriv_coeffs[0].resize(p+1,T(0.0));
  if (p > 1) {
    deriv_coeffs[1].resize(p+1,T(0.0));
    deriv_coeffs[1][0] = T(1.0);
  }
  if (p > 2) {
    deriv_coeffs[2].resize(p+1,T(0.0));
    deriv_coeffs[2][1] = T(3.0);
  }
  for (unsigned int k=3; k<=p; k++) {
    deriv_coeffs[k].resize(p+1,T(0.0));
    deriv_coeffs[k][0] = T(1.0/3.0)*deriv_coeffs[k-1][1];
    for (unsigned int i=1; i<=k-3; i++)
      deriv_coeffs[k][i] = T(i/(2.0*i - 1.0))*deriv_coeffs[k-1][i-1] + 
	T((i+1.0)/(2.0*i + 3.0))*deriv_coeffs[k-1][i+1];
    deriv_coeffs[k][k-2] = T((k-2.0)/(2.0*k-5.0))*deriv_coeffs[k-1][k-3];
    deriv_coeffs[k][k-1] = T((k-1.0)/(2.0*k-3.0))*deriv_coeffs[k-1][k-2] + k;
  }
}

template <typename T>
Stokhos::LegendreBasis<T>::
~LegendreBasis()
{
}

template <typename T>
void
Stokhos::LegendreBasis<T>::
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
    h[0] = T(1.0/3.0)*coeffs[1] + x.coeff(k);
    for (unsigned int i=1; i<=pc-1; i++) {
      h[i] = T(i / (2.0*i - 1.0))*coeffs[i-1] + 
	T((i+1.0)/(2.0*i + 3.0))*coeffs[i+1];
    }
    h[pc] = T( pc / (2.0*pc - 1.0) )*coeffs[pc-1];
    h[pc+1] = T( (pc+1.0)/(2.0*pc + 1.0) )*coeffs[pc];
    
    // Copy into coeffs
    for (unsigned int i=0; i<=pc+1; i++)
      coeffs[i] = h[i];

    pc = pc+1;
  }
}

template <typename T>
void
Stokhos::LegendreBasis<T>::
projectDerivative(unsigned int i, std::vector<T>& coeffs) const
{
  coeffs = deriv_coeffs[i];
}

template <typename T>
void
Stokhos::LegendreBasis<T>::
evaluateBases(const std::vector<T>& point, std::vector<T>& basis_pts) const
{
  const T& x = point[0];

  // Evaluate basis polynomials P(x) using 3 term recurrence
  // P_0(x) = 1
  // P_1(x) = x
  // P_i(x) = (2*i-1)/i*x*P_{i-1}(x) - (i-1)/i*P_{i-2}(x), i=2,3,...
  basis_pts[0] = T(1.0);
  if (this->p >= 1)
    basis_pts[1] = x;
  for (unsigned int i=2; i<=this->p; i++)
    basis_pts[i] = T(2*i-1)/T(i)*x*basis_pts[i-1] - T(i-1)/T(i)*basis_pts[i-2];
}
