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

template <typename ordinal_type, typename value_type>
Stokhos::ClenshawCurtisLegendreBasis<ordinal_type, value_type>::
ClenshawCurtisLegendreBasis(ordinal_type p, bool normalize, bool isotropic_) :
  RecurrenceBasis<ordinal_type, value_type>("Clenshaw-Curtis Legendre", 
					    p, normalize),
  isotropic(isotropic_)
{
  // Compute coefficients in 3-term recurrsion
  computeRecurrenceCoefficients(p+1, this->alpha, this->beta, this->delta);

  // Setup rest of recurrence basis
  this->setup();

#ifdef HAVE_STOKHOS_DAKOTA
  this->setSparseGridRule(Pecos::CLENSHAW_CURTIS);
  if (isotropic)
    this->setSparseGridGrowthRule(Pecos::FULL_EXPONENTIAL);
  else
    this->setSparseGridGrowthRule(Pecos::MODERATE_EXPONENTIAL);
#endif
}

template <typename ordinal_type, typename value_type>
Stokhos::ClenshawCurtisLegendreBasis<ordinal_type, value_type>::
ClenshawCurtisLegendreBasis(ordinal_type p, 
			    const ClenshawCurtisLegendreBasis& basis) :
  RecurrenceBasis<ordinal_type, value_type>(p, basis),
  isotropic(basis.isotropic)
{
  // Compute coefficients in 3-term recurrsion
  computeRecurrenceCoefficients(p+1, this->alpha, this->beta, this->delta);

  // Setup rest of recurrence basis
  this->setup();
}

template <typename ordinal_type, typename value_type>
Stokhos::ClenshawCurtisLegendreBasis<ordinal_type, value_type>::
~ClenshawCurtisLegendreBasis()
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ClenshawCurtisLegendreBasis<ordinal_type, value_type>::
computeRecurrenceCoefficients(ordinal_type n,
			      Teuchos::Array<value_type>& alpha,
			      Teuchos::Array<value_type>& beta,
			      Teuchos::Array<value_type>& delta) const
{
  // Legendre 3 term recurrence:
  // P_0(x) = 1
  // P_1(x) = x
  // P_i(x) = (2*i-1)/i*x*P_{i-1}(x) - (i-1)/i*P_{i-2}(x), i=2,3,...
  alpha[0] = 0.0;
  beta[0] = 1.0;
  delta[0] = 1.0;
  for (ordinal_type i=1; i<n; i++) {
    alpha[i] = 0.0;
    //beta[i] = value_type(i*i) / value_type((2*i-1)*(2*i+1));
    beta[i] = value_type(i) / value_type(i+1);
    delta[i] = value_type(2*i+1) / value_type(i+1);
  }
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > 
Stokhos::ClenshawCurtisLegendreBasis<ordinal_type,value_type>::
cloneWithOrder(ordinal_type p) const
{
  return 
    Teuchos::rcp(new Stokhos::ClenshawCurtisLegendreBasis<ordinal_type,value_type>(p,*this));
}
