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

#ifdef HAVE_STOKHOS_DAKOTA
#include "sandia_rules.hpp"
#endif
#include "Teuchos_TestForException.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::GaussPattersonLegendreBasis<ordinal_type, value_type>::
GaussPattersonLegendreBasis(ordinal_type p, bool normalize, bool isotropic_) :
  LegendreBasis<ordinal_type, value_type>(p, normalize),
  isotropic(isotropic_)
{
#ifdef HAVE_STOKHOS_DAKOTA
  this->setSparseGridGrowthRule(webbur::level_to_order_exp_gp);
#endif
}

template <typename ordinal_type, typename value_type>
Stokhos::GaussPattersonLegendreBasis<ordinal_type, value_type>::
GaussPattersonLegendreBasis(ordinal_type p, 
			    const GaussPattersonLegendreBasis& basis) :
  LegendreBasis<ordinal_type, value_type>(p, basis),
  isotropic(basis.isotropic)
{
}

template <typename ordinal_type, typename value_type>
Stokhos::GaussPattersonLegendreBasis<ordinal_type, value_type>::
~GaussPattersonLegendreBasis()
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::GaussPattersonLegendreBasis<ordinal_type,value_type>::
getQuadPoints(ordinal_type quad_order,
	      Teuchos::Array<value_type>& quad_points,
	      Teuchos::Array<value_type>& quad_weights,
	      Teuchos::Array< Teuchos::Array<value_type> >& quad_values) const
{
#ifdef HAVE_STOKHOS_DAKOTA
  // Gauss-Patterson points have the following structure 
  // (cf. http://people.sc.fsu.edu/~jburkardt/f_src/patterson_rule/patterson_rule.html):
  //  Level l  Num points n   Precision p
  //  -----------------------------------
  //    0        1                1
  //    1        3                5
  //    2        7                11
  //    3        15               23
  //    4        31               47
  //    5        63               95
  //    6        127              191
  //    7        255              383
  // Thus for l > 0, n = 2^{l+1}-1 and p = 3*2^l-1.  So for a given quadrature
  // order p, we find the smallest l s.t. 3*s^l-1 >= p and then compute the
  // number of points n from the above.  In this case, l = ceil(log2((p+1)/3))
  ordinal_type num_points;
  if (quad_order <= ordinal_type(1))
    num_points = 1;
  else {
    ordinal_type l = std::ceil(std::log((quad_order+1.0)/3.0)/std::log(2.0));
    num_points = (1 << (l+1)) - 1; // std::pow(2,l+1)-1;
  }
  
  quad_points.resize(num_points);
  quad_weights.resize(num_points);
  quad_values.resize(num_points);

  webbur::patterson_lookup(num_points, &quad_points[0], &quad_weights[0]);

  for (ordinal_type i=0; i<num_points; i++) {
    quad_weights[i] *= 0.5;  // scale to unit measure
    quad_values[i].resize(this->p+1);
    this->evaluateBases(quad_points[i], quad_values[i]);
  }

#else
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error, "Clenshaw-Curtis requires TriKota to be enabled!");
#endif
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::GaussPattersonLegendreBasis<ordinal_type,value_type>::
quadDegreeOfExactness(ordinal_type n) const
{
  // Based on the above structure, we find the largest l s.t. 2^{l+1}-1 <= n,
  // which is floor(log2(n+1)-1) and compute p = 3*2^l-1
  if (n == ordinal_type(1))
    return 1;
  ordinal_type l = std::floor(std::log(n+1.0)/std::log(2.0)-1.0);
  return (3 << l) - 1; // 3*std::pow(2,l)-1;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > 
Stokhos::GaussPattersonLegendreBasis<ordinal_type,value_type>::
cloneWithOrder(ordinal_type p) const
{
  return 
    Teuchos::rcp(new Stokhos::GaussPattersonLegendreBasis<ordinal_type,value_type>(p,*this));
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::GaussPattersonLegendreBasis<ordinal_type,value_type>::
coefficientGrowth(ordinal_type n) const
{
  // Gauss-Patterson rules have precision 3*2^l-1, which is odd.
  // Since discrete orthogonality requires integrating polynomials of
  // order 2*p, setting p = 3*2^{l-1}-1 will yield the largest p such that
  // 2*p <= 3*2^l-1
  if (n == 0) 
    return 0;
  return (3 << (n-1)) - 1; // 3*std::pow(2,n-1) - 1;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::GaussPattersonLegendreBasis<ordinal_type,value_type>::
pointGrowth(ordinal_type n) const
{
  return n;
}
