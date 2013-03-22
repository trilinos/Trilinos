// $Id: Stokhos_JacobiBasisImp.hpp,v 1.1.1.1 2010/02/10 20:22:35 kevin Exp $ 
// $Source: /usr/local/cvs/UQ/Ops/Stokhos_JacobiBasisImp.hpp,v $ 
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
Stokhos::JacobiBasis<ordinal_type, value_type>::
JacobiBasis(ordinal_type p,  
  value_type alphaIndex, 
  value_type betaIndex, bool normalize, Stokhos::GrowthPolicy growth) :
  RecurrenceBasis<ordinal_type, value_type>("Jacobi", p, normalize, growth),
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
JacobiBasis(ordinal_type p, const JacobiBasis& basis) :
  RecurrenceBasis<ordinal_type, value_type>(p, basis),
  alphaIndex_(basis.alphaIndex_),
  betaIndex_(basis.betaIndex_)
{
  // Compute coefficients in 3-term recurrsion
  computeRecurrenceCoefficients(p+1, this->alpha, this->beta, this->delta,
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
			      Teuchos::Array<value_type>& alpha,
			      Teuchos::Array<value_type>& beta,
			      Teuchos::Array<value_type>& delta,
			      Teuchos::Array<value_type>& gamma) const
{
  value_type a = alphaIndex_;
  value_type b = betaIndex_;

  if (a==0.0 && b==0.0)
  {
    alpha[0] = 0.0;
    beta[0] = 1.0;
    delta[0] = 1.0;
    gamma[0] = 1.0;
  }
  else
  {
    alpha[0] = getB(0)/getA(0);
    beta[0] = 1.0;
    delta[0] = getC(0)/getA(0);
    gamma[0] = 1.0;
  }
  for (ordinal_type i=1; i<n; i++) 
  {
    alpha[i] = getB(i)/getA(i);
    beta[i] = getD(i)/getA(i);
    delta[i] = getC(i)/getA(i);
    gamma[i] = 1.0;
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
cloneWithOrder(ordinal_type p) const
{
  return 
    Teuchos::rcp(new Stokhos::JacobiBasis<ordinal_type,value_type>(p,*this));
}
