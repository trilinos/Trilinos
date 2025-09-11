// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef HAVE_STOKHOS_DAKOTA
#include "sandia_rules.hpp"
#endif
#include "Teuchos_TestForException.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::ClenshawCurtisLegendreBasis<ordinal_type, value_type>::
ClenshawCurtisLegendreBasis(ordinal_type ap, bool anormalize, bool isotropic_) :
  LegendreBasis<ordinal_type, value_type>(ap, anormalize),
  isotropic(isotropic_)
{
#ifdef HAVE_STOKHOS_DAKOTA
  this->setSparseGridGrowthRule(webbur::level_to_order_exp_cc);
#endif
}

template <typename ordinal_type, typename value_type>
Stokhos::ClenshawCurtisLegendreBasis<ordinal_type, value_type>::
ClenshawCurtisLegendreBasis(ordinal_type ap, 
			    const ClenshawCurtisLegendreBasis& basis) :
  LegendreBasis<ordinal_type, value_type>(ap, basis),
  isotropic(basis.isotropic)
{
}

template <typename ordinal_type, typename value_type>
Stokhos::ClenshawCurtisLegendreBasis<ordinal_type, value_type>::
~ClenshawCurtisLegendreBasis()
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ClenshawCurtisLegendreBasis<ordinal_type,value_type>::
getQuadPoints(ordinal_type quad_order,
	      Teuchos::Array<value_type>& quad_points,
	      Teuchos::Array<value_type>& quad_weights,
	      Teuchos::Array< Teuchos::Array<value_type> >& quad_values) const
{
#ifdef HAVE_STOKHOS_DAKOTA
  ordinal_type num_points;
  if (quad_order % ordinal_type(2) == ordinal_type(1))
    num_points = quad_order;
  else
    num_points = quad_order+1;
  quad_points.resize(num_points);
  quad_weights.resize(num_points);
  quad_values.resize(num_points);

  webbur::clenshaw_curtis_compute(
    num_points, &quad_points[0], &quad_weights[0]);

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
Stokhos::ClenshawCurtisLegendreBasis<ordinal_type,value_type>::
quadDegreeOfExactness(ordinal_type n) const
{
  if (n % ordinal_type(2) == ordinal_type(1))
    return n;
  return n-ordinal_type(1);
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > 
Stokhos::ClenshawCurtisLegendreBasis<ordinal_type,value_type>::
cloneWithOrder(ordinal_type ap) const
{
  return 
    Teuchos::rcp(new Stokhos::ClenshawCurtisLegendreBasis<ordinal_type,value_type>(ap,*this));
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::ClenshawCurtisLegendreBasis<ordinal_type,value_type>::
coefficientGrowth(ordinal_type n) const
{
  if (n == ordinal_type(0)) 
    return 0;
  return (1 << (n-1)); // std::pow(2,n-1);
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::ClenshawCurtisLegendreBasis<ordinal_type,value_type>::
pointGrowth(ordinal_type n) const
{
  return n;
}
