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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
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
