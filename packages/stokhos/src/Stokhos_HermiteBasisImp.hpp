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
#else
#ifdef HAVE_STOKHOS_FORUQTK
#include "Stokhos_gaussq.h"
#else
#include "Teuchos_TestForException.hpp"
#endif
#endif

template <typename ordinal_type, typename value_type>
Stokhos::HermiteBasis<ordinal_type,value_type>::
HermiteBasis(ordinal_type p) :
  OneDOrthogPolyBasisBase<ordinal_type,value_type>("Hermite",p)
{
  // Fill in basis coefficients
  this->basis[0].coeff(0) = value_type(1.0);
  if (this->p >= 1)
    this->basis[1].coeff(1) = value_type(1.0);
  for (ordinal_type k=2; k<=this->p; k++) {
    this->basis[k].coeff(0) = -value_type(k-1)*(this->basis[k-2].coeff(0));
    for (ordinal_type i=1; i<=k; i++)
      this->basis[k].coeff(i) = 
	this->basis[k-1].coeff(i-1) - value_type(k-1)*(this->basis[k-2].coeff(i));
  }

  // Compute norms
  this->norms[0] = value_type(1.0);
  for (ordinal_type k=1; k<=this->p; k++)
    this->norms[k] = value_type(k)*(this->norms[k-1]);
}

template <typename ordinal_type, typename value_type>
Stokhos::HermiteBasis<ordinal_type,value_type>::
~HermiteBasis()
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::HermiteBasis<ordinal_type,value_type>::
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
    h[0] = coeffs[1] + x.coeff(k);
    for (ordinal_type i=1; i<=pc-1; i++) {
      h[i] = coeffs[i-1] + value_type(i+1)*coeffs[i+1];
    }
    h[pc] = coeffs[pc-1];
    h[pc+1] = coeffs[pc];
    
    // Copy into coeffs
    for (ordinal_type i=0; i<=pc+1; i++)
      coeffs[i] = h[i];

    pc = pc+1;
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::HermiteBasis<ordinal_type,value_type>::
projectDerivative(ordinal_type i, std::vector<value_type>& coeffs) const
{
  for (ordinal_type j=0; j<static_cast<ordinal_type>(coeffs.size()); j++)
    coeffs[j] = value_type(0.0);
  if (i > 0)
    coeffs[i-1] = value_type(i);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::HermiteBasis<ordinal_type,value_type>::
evaluateBases(const value_type& x, std::vector<value_type>& basis_pts) const
{
  // Evaluate basis polynomials He(x) using 3 term recurrence
  // He_0(x) = 1
  // He_1(x) = x
  // He_i(x) = x*He_{i-1}(x) - (i-1)*He_{i-2}(x), i=2,3,...
  basis_pts[0] = value_type(1.0);
  if (this->p >= 1)
    basis_pts[1] = x;
  for (ordinal_type i=2; i<=this->p; i++)
    basis_pts[i] = x*basis_pts[i-1] - value_type(i-1)*basis_pts[i-2];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::HermiteBasis<ordinal_type,value_type>::
getQuadPoints(ordinal_type quad_order,
	      std::vector<value_type>& quad_points,
	      std::vector<value_type>& quad_weights,
	      std::vector< std::vector<value_type> >& quad_values) const
{
  // Compute gauss points, weights
  ordinal_type n = static_cast<ordinal_type>(std::ceil((quad_order+1)/2.0));
  std::vector<double> x(n), w(n);
  
#ifdef HAVE_STOKHOS_DAKOTA
  webbur::hermite_compute(n, &x[0], &w[0]);
#else
#ifdef HAVE_STOKHOS_FORUQTK
  int kind = 4;
  int kpts = 0;
  double endpts[2] = {0.0, 0.0};
  std::vector<double> b(n);
  int ni = n;
  double alpha = 0.0;
  double beta = 0.0;
  GAUSSQ_F77(&kind, &ni, &alpha, &beta, &kpts, endpts, &b[0], &x[0], &w[0]);
#else
  TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Stokhos::HermiteBasis::getQuadPoints():  "
		     << " Must have Dakota or ForUQTK enabled for quadrature!");
#endif
#endif

  quad_points.resize(n);
  quad_weights.resize(n);
  for (ordinal_type i=0; i<n; i++) {
    quad_points[i] = x[i];
    quad_weights[i] = w[i];
  }

  // gaussq gives points and weights for w(x) = exp(-x*x).  We need points
  // and weights for w(x) = exp(-x*x/2)/sqrt(2*pi)
  value_type weight_factor = getQuadWeightFactor();
  value_type point_factor = getQuadPointFactor();
  for (ordinal_type i=0; i<n; i++) {
    quad_weights[i] *= weight_factor;
    quad_points[i] *= point_factor;
  }

  // Evalute basis at gauss points
  quad_values.resize(n);
  for (ordinal_type i=0; i<n; i++) {
    quad_values[i].resize(this->p+1);
    evaluateBases(quad_points[i], quad_values[i]);
  }

  // std::cout << "Hermite quadrature points, weights, values = " << std::endl;
  // for (ordinal_type i=0; i<n; i++) {
  //   std::cout << "\t" << quad_points[i] 
  //             << "\t" << quad_weights[i];
  //   for (ordinal_type j=0; j<p+1; j++)
  //     std::cout << "\t" << quad_values[i][j];
  //   std::cout << std::endl;
  // }
}
