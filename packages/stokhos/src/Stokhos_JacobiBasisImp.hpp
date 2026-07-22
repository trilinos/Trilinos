// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

template <typename ordinal_type, typename value_type>
Stokhos::JacobiBasis<ordinal_type, value_type>::
JacobiBasis(ordinal_type ap,  
  value_type alphaIndex, 
  value_type betaIndex, bool anormalize, Stokhos::GrowthPolicy agrowth) :
  RecurrenceBasis<ordinal_type, value_type>("Jacobi", ap, anormalize, agrowth),
  alphaIndex_(alphaIndex),
  betaIndex_(betaIndex)
{
  this->setup();

#ifdef HAVE_STOKHOS_DAKOTA
  this->setSparseGridGrowthRule(webbur::level_to_order_linear_wn);
#endif
}

template <typename ordinal_type, typename value_type>
Stokhos::JacobiBasis<ordinal_type, value_type>::
JacobiBasis(ordinal_type ap, const JacobiBasis& basis) :
  RecurrenceBasis<ordinal_type, value_type>(ap, basis),
  alphaIndex_(basis.alphaIndex_),
  betaIndex_(basis.betaIndex_)
{
  // Compute coefficients in 3-term recurrsion
  computeRecurrenceCoefficients(ap+1, this->alpha, this->beta, this->delta,
				this->gamma);

  // Setup rest of recurrence basis
  this->setup();
}

template <typename ordinal_type, typename value_type>
Stokhos::JacobiBasis<ordinal_type, value_type>::
~JacobiBasis()
{
}

template <typename ordinal_type, typename value_type>
bool
Stokhos::JacobiBasis<ordinal_type, value_type>::
computeRecurrenceCoefficients(ordinal_type n,
			      Teuchos::Array<value_type>& aalpha,
			      Teuchos::Array<value_type>& abeta,
			      Teuchos::Array<value_type>& adelta,
			      Teuchos::Array<value_type>& agamma) const
{
  value_type a = alphaIndex_;
  value_type b = betaIndex_;

  if (a==0.0 && b==0.0)
  {
    aalpha[0] = 0.0;
    abeta[0] = 1.0;
    adelta[0] = 1.0;
    agamma[0] = 1.0;
  }
  else
  {
    aalpha[0] = getB(0)/getA(0);
    abeta[0] = 1.0;
    adelta[0] = getC(0)/getA(0);
    agamma[0] = 1.0;
  }
  for (ordinal_type i=1; i<n; i++) 
  {
    aalpha[i] = getB(i)/getA(i);
    abeta[i] = getD(i)/getA(i);
    adelta[i] = getC(i)/getA(i);
    agamma[i] = 1.0;
  }

  return false;
}


template <typename ordinal_type, typename value_type>
value_type Stokhos::JacobiBasis<ordinal_type, value_type>::getA(int n) const
{
  return 2*(n+1)*(n+alphaIndex_+betaIndex_+1)*(2*n+alphaIndex_+betaIndex_);
}

template <typename ordinal_type, typename value_type>
value_type Stokhos::JacobiBasis<ordinal_type, value_type>::getB(int n) const
{
  value_type a = alphaIndex_;
  value_type b = betaIndex_;
  return -(2*n+a+b+1)*(a*a-b*b);
}

template <typename ordinal_type, typename value_type>
value_type Stokhos::JacobiBasis<ordinal_type, value_type>::getC(int n) const
{
  value_type a = alphaIndex_;
  value_type b = betaIndex_;
  return poch3(2*n+a+b);
}

template <typename ordinal_type, typename value_type>
value_type Stokhos::JacobiBasis<ordinal_type, value_type>::getD(int n) const
{
  value_type a = alphaIndex_;
  value_type b = betaIndex_;
  return 2*(n+a)*(n+b)*(2*n + a + b + 2);
}

template <typename ordinal_type, typename value_type>
value_type Stokhos::JacobiBasis<ordinal_type, value_type>::poch3(value_type m) const
{
  return (m+2)*(m+1)*m;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > 
Stokhos::JacobiBasis<ordinal_type,value_type>::
cloneWithOrder(ordinal_type ap) const
{
  return 
    Teuchos::rcp(new Stokhos::JacobiBasis<ordinal_type,value_type>(ap,*this));
}
