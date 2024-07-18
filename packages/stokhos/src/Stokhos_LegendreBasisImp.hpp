// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

template <typename ordinal_type, typename value_type>
Stokhos::LegendreBasis<ordinal_type, value_type>::
LegendreBasis(ordinal_type p, bool normalize, Stokhos::GrowthPolicy growth) :
  RecurrenceBasis<ordinal_type, value_type>("Legendre", p, normalize, growth)
{
  this->setup();

#ifdef HAVE_STOKHOS_DAKOTA
  this->setSparseGridGrowthRule(webbur::level_to_order_linear_wn);
#endif
}

template <typename ordinal_type, typename value_type>
Stokhos::LegendreBasis<ordinal_type, value_type>::
LegendreBasis(ordinal_type p, const LegendreBasis& basis) :
  RecurrenceBasis<ordinal_type, value_type>(p, basis)
{
  // Compute coefficients in 3-term recurrsion
  computeRecurrenceCoefficients(p+1, this->alpha, this->beta, this->delta,
                                this->gamma);

  // Setup rest of recurrence basis
  this->setup();
}

template <typename ordinal_type, typename value_type>
Stokhos::LegendreBasis<ordinal_type, value_type>::
~LegendreBasis()
{
}

template <typename ordinal_type, typename value_type>
bool
Stokhos::LegendreBasis<ordinal_type, value_type>::
computeRecurrenceCoefficients(ordinal_type n,
                              Teuchos::Array<value_type>& theAlpha,
                              Teuchos::Array<value_type>& theBeta,
                              Teuchos::Array<value_type>& theDelta,
                              Teuchos::Array<value_type>& theGamma) const
{
  // Legendre 3 term recurrence:
  // P_0(x) = 1
  // P_1(x) = x
  // P_i(x) = (2*i-1)/i*x*P_{i-1}(x) - (i-1)/i*P_{i-2}(x), i=2,3,...
  theAlpha[0] = 0.0;
  theBeta[0] = 1.0;
  theDelta[0] = 1.0;
  theGamma[0] = 1.0;
  for (ordinal_type i=1; i<n; i++) {
    theAlpha[i] = 0.0;
    //theBeta[i] = value_type(i*i) / value_type((2*i-1)*(2*i+1));
    theBeta[i] = value_type(i) / value_type(i+1);
    theDelta[i] = value_type(2*i+1) / value_type(i+1);
    theGamma[i] = 1.0;
  }

  return false;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> >
Stokhos::LegendreBasis<ordinal_type,value_type>::
cloneWithOrder(ordinal_type pp) const
{
  return
    Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(pp,*this));
}
