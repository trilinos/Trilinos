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
HermiteBasis(ordinal_type ap, bool anormalize, Stokhos::GrowthPolicy agrowth) :
  RecurrenceBasis<ordinal_type,value_type>("Hermite", ap, anormalize, agrowth)
{
  this->setup();

#ifdef HAVE_STOKHOS_DAKOTA
  this->setSparseGridGrowthRule(webbur::level_to_order_linear_wn);
#endif
}

template <typename ordinal_type, typename value_type>
Stokhos::HermiteBasis<ordinal_type, value_type>::
HermiteBasis(ordinal_type ap, const HermiteBasis& basis) :
  RecurrenceBasis<ordinal_type, value_type>(ap, basis)
{
  // Compute coefficients in 3-term recurrsion
  computeRecurrenceCoefficients(ap+1, this->alpha, this->beta, this->delta,
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
			      Teuchos::Array<value_type>& aalpha,
			      Teuchos::Array<value_type>& abeta,
			      Teuchos::Array<value_type>& adelta,
			      Teuchos::Array<value_type>& agamma) const
{
  // Hermite 3 term recurrence:
  // He_0(x) = 1
  // He_1(x) = x
  // He_i(x) = x*He_{i-1}(x) - (i-1)*He_{i-2}(x), i=2,3,...
  aalpha[0] = 0.0;
  abeta[0] = 1.0;
  adelta[0] = 1.0;
  agamma[0] = 1.0;
  for (ordinal_type i=1; i<n; i++) {
    aalpha[i] = 0.0;
    abeta[i] = value_type(i);
    adelta[i] = 1.0;
    agamma[i] = 1.0;
  }

  return false;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > 
Stokhos::HermiteBasis<ordinal_type,value_type>::
cloneWithOrder(ordinal_type ap) const
{
  return 
    Teuchos::rcp(new Stokhos::HermiteBasis<ordinal_type,value_type>(ap,*this));
}
