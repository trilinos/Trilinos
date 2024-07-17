// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

template <typename ordinal_type, typename value_type>
Stokhos::HermiteBasis<ordinal_type,value_type>::
HermiteBasis(ordinal_type p, bool normalize, Stokhos::GrowthPolicy growth) :
  RecurrenceBasis<ordinal_type,value_type>("Hermite", p, normalize, growth)
{
  this->setup();

#ifdef HAVE_STOKHOS_DAKOTA
  this->setSparseGridGrowthRule(webbur::level_to_order_linear_wn);
#endif
}

template <typename ordinal_type, typename value_type>
Stokhos::HermiteBasis<ordinal_type, value_type>::
HermiteBasis(ordinal_type p, const HermiteBasis& basis) :
  RecurrenceBasis<ordinal_type, value_type>(p, basis)
{
  // Compute coefficients in 3-term recurrsion
  computeRecurrenceCoefficients(p+1, this->alpha, this->beta, this->delta,
				this->gamma);

  // Setup rest of recurrence basis
  this->setup();
}

template <typename ordinal_type, typename value_type>
Stokhos::HermiteBasis<ordinal_type,value_type>::
~HermiteBasis()
{
}

template <typename ordinal_type, typename value_type>
bool
Stokhos::HermiteBasis<ordinal_type,value_type>::
computeRecurrenceCoefficients(ordinal_type n,
			      Teuchos::Array<value_type>& alpha,
			      Teuchos::Array<value_type>& beta,
			      Teuchos::Array<value_type>& delta,
			      Teuchos::Array<value_type>& gamma) const
{
  // Hermite 3 term recurrence:
  // He_0(x) = 1
  // He_1(x) = x
  // He_i(x) = x*He_{i-1}(x) - (i-1)*He_{i-2}(x), i=2,3,...
  alpha[0] = 0.0;
  beta[0] = 1.0;
  delta[0] = 1.0;
  gamma[0] = 1.0;
  for (ordinal_type i=1; i<n; i++) {
    alpha[i] = 0.0;
    beta[i] = value_type(i);
    delta[i] = 1.0;
    gamma[i] = 1.0;
  }

  return false;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > 
Stokhos::HermiteBasis<ordinal_type,value_type>::
cloneWithOrder(ordinal_type p) const
{
  return 
    Teuchos::rcp(new Stokhos::HermiteBasis<ordinal_type,value_type>(p,*this));
}
