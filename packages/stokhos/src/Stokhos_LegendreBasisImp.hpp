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

#include "Stokhos_ConfigDefs.h"
#ifdef HAVE_STOKHOS_DAKOTA
#include "sandia_rules.H"
#elif HAVE_STOKHOS_FORUQTK
#include "Stokhos_gaussq.h"
#else
#include "Teuchos_TestForException.hpp"
#endif

template <typename ordinal_type, typename value_type>
Stokhos::LegendreBasis<ordinal_type, value_type>::
LegendreBasis(ordinal_type p) :
  OneDOrthogPolyBasisBase<ordinal_type, value_type>("Legendre",p),
  deriv_coeffs(p+1)
{
  // Fill in basis coefficients
  value_type c1, c2;
  this->basis[0].coeff(0) = value_type(1.0);
  if (this->p >= 1)
    this->basis[1].coeff(1) = value_type(1.0);
  for (ordinal_type k=2; k<=this->p; k++) {
    c1 = value_type( (2.0*k-1.0) / k );
    c2 = value_type( (k-1.0) / k );
    this->basis[k].coeff(0) = -c2*(this->basis[k-2].coeff(0));
    for (ordinal_type i=1; i<=k; i++)
      this->basis[k].coeff(i) = 
        c1*(this->basis[k-1].coeff(i-1)) - c2*(this->basis[k-2].coeff(i));
  }

  // Compute norms
  for (ordinal_type k=0; k<=this->p; k++)
    this->norms[k] = value_type(1.0/(2.0*k+1.0));

  // Compute derivative coefficients
  deriv_coeffs[0].resize(p+1,value_type(0.0));
  if (p >= 1) {
    deriv_coeffs[1].resize(p+1,value_type(0.0));
    deriv_coeffs[1][0] = value_type(1.0);
  }
  if (p >= 2) {
    deriv_coeffs[2].resize(p+1,value_type(0.0));
    deriv_coeffs[2][1] = value_type(3.0);
  }
  for (ordinal_type k=3; k<=p; k++) {
    deriv_coeffs[k].resize(p+1,value_type(0.0));
    deriv_coeffs[k][0] = value_type(1.0/3.0)*deriv_coeffs[k-1][1];
    for (ordinal_type i=1; i<=k-3; i++)
      deriv_coeffs[k][i] = value_type(i/(2.0*i - 1.0))*deriv_coeffs[k-1][i-1] + 
        value_type((i+1.0)/(2.0*i + 3.0))*deriv_coeffs[k-1][i+1];
    deriv_coeffs[k][k-2] = value_type((k-2.0)/(2.0*k-5.0))*deriv_coeffs[k-1][k-3];
    deriv_coeffs[k][k-1] = value_type((k-1.0)/(2.0*k-3.0))*deriv_coeffs[k-1][k-2] + k;
  }
}

template <typename ordinal_type, typename value_type>
Stokhos::LegendreBasis<ordinal_type, value_type>::
~LegendreBasis()
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LegendreBasis<ordinal_type, value_type>::
projectPoly(const Stokhos::Polynomial<value_type>& x, std::vector<value_type>& coeffs) const
{
  // Initialize
  for (ordinal_type i=0; i<=this->p; i++)
    coeffs[i] = value_type(0.0);

  ordinal_type px = x.degree();
  if (px > this->p && px != this->p*2)
    px = this->p;

  // Handle degree 0 case
  if (px == 0) {
    coeffs[0] = x.coeff(0);
    return;
  }

  // Temporary array
  std::vector<value_type> h(px+1,value_type(0.));

  coeffs[0] = x.coeff(px-1);
  coeffs[1] = x.coeff(px);
  ordinal_type pc = 1;

  for (int k=px-2; k>=0; --k) {

    // Multiply by t
    h[0] = value_type(1.0/3.0)*coeffs[1] + x.coeff(k);
    for (ordinal_type i=1; i<=pc-1; i++) {
      h[i] = value_type(i / (2.0*i - 1.0))*coeffs[i-1] + 
        value_type((i+1.0)/(2.0*i + 3.0))*coeffs[i+1];
    }
    h[pc] = value_type( pc / (2.0*pc - 1.0) )*coeffs[pc-1];
    h[pc+1] = value_type( (pc+1.0)/(2.0*pc + 1.0) )*coeffs[pc];
    
    // Copy into coeffs
    for (ordinal_type i=0; i<=pc+1; i++)
      coeffs[i] = h[i];

    pc = pc+1;
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LegendreBasis<ordinal_type, value_type>::
projectDerivative(ordinal_type i, std::vector<value_type>& coeffs) const
{
  coeffs = deriv_coeffs[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LegendreBasis<ordinal_type, value_type>::
evaluateBases(const value_type& x, std::vector<value_type>& basis_pts) const
{
  // Evaluate basis polynomials P(x) using 3 term recurrence
  // P_0(x) = 1
  // P_1(x) = x
  // P_i(x) = (2*i-1)/i*x*P_{i-1}(x) - (i-1)/i*P_{i-2}(x), i=2,3,...
  basis_pts[0] = value_type(1.0);
  if (this->p >= 1)
    basis_pts[1] = x;
  for (ordinal_type i=2; i<=this->p; i++)
    basis_pts[i] = value_type(2*i-1)/value_type(i)*x*basis_pts[i-1] - value_type(i-1)/value_type(i)*basis_pts[i-2];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LegendreBasis<ordinal_type, value_type>::
getQuadPoints(ordinal_type quad_order,
	      std::vector<value_type>& quad_points,
	      std::vector<value_type>& quad_weights,
	      std::vector< std::vector<value_type> >& quad_values) const
{
  // Compute gauss points, weights
  ordinal_type n = static_cast<ordinal_type>(std::ceil((quad_order+1)/2.0));
  std::vector<double> x(n), w(n);
  double alpha = 0.0;
  double beta = 0.0;
  
#ifdef HAVE_STOKHOS_DAKOTA
  webbur::legendre_compute(n, alpha, beta, &x[0], &w[0]);
#elif HAVE_STOKHOS_FORUQTK
  int kind = 1;
  int kpts = 0;
  double endpts[2] = {0.0, 0.0};
  std::vector<double> b(n);
  int ni = n;
  GAUSSQ_F77(&kind, &ni, &alpha, &beta, &kpts, endpts, &b[0], &x[0], &w[0]);
#else
  TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Stokhos::LegendreBasis::getQuadPoints():  "
		     << " Must have Dakota or ForUQTK enabled for quadrature!");
#endif

  quad_points.resize(n);
  quad_weights.resize(n);
  for (ordinal_type i=0; i<n; i++) {
    quad_points[i] = x[i];
    quad_weights[i] = w[i];
  }

  // Evalute basis at gauss points
  quad_values.resize(n);
  for (ordinal_type i=0; i<n; i++) {
    quad_values[i].resize(this->p+1);
    evaluateBases(quad_points[i], quad_values[i]);
  }

  // Scale weights by density at gauss points
  value_type weight_factor = getQuadWeightFactor();
  value_type point_factor = getQuadPointFactor();
  for (ordinal_type i=0; i<n; i++) {
    quad_weights[i] *= weight_factor;
    quad_points[i] *= point_factor;
  }

//   std::cout << "Legendre quadrature points, weights, values = " << std::endl;
//   for (int i=0; i<n; i++) {
//     std::cout << "\t" << quad_points[i] 
//               << "\t" << quad_weights[i];
//     for (ordinal_type j=0; j<p+1; j++)
//       std::cout << "\t" << quad_values[i][j];
//     cout << std::endl;
//   }
}
