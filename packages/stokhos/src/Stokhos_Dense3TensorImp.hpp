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
Stokhos::Dense3Tensor<ordinal_type, value_type>::
Dense3Tensor(ordinal_type sz) :
  l(sz),
  Cijk_values(l*l*l)
{
}

template <typename ordinal_type, typename value_type>
Stokhos::Dense3Tensor<ordinal_type, value_type>::
~Dense3Tensor()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::Dense3Tensor<ordinal_type, value_type>::
size() const
{
  return l;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::Dense3Tensor<ordinal_type, value_type>::
operator() (ordinal_type i, ordinal_type j, ordinal_type k) const
{
  return Cijk_values[ l*(l*k + j) + i ];
}

template <typename ordinal_type, typename value_type>
value_type&
Stokhos::Dense3Tensor<ordinal_type, value_type>::
operator() (ordinal_type i, ordinal_type j, ordinal_type k)
{
  return Cijk_values[ l*(l*k + j) + i ];
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::Dense3Tensor<ordinal_type, value_type>::
num_values(ordinal_type k) const
{
  return l*l;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::Dense3Tensor<ordinal_type, value_type>::
value(ordinal_type k, ordinal_type ll, 
      ordinal_type& i, ordinal_type& j, value_type& c) const
{
  j = ll/l;
  i = ll-j*l;
  c = triple_value(i,j,k);
}

template <typename ordinal_type, typename value_type>
void 
Stokhos::Dense3Tensor<ordinal_type, value_type>::
print(std::ostream& os) const
{
  for (ordinal_type i=0; i<l; i++)
    for (ordinal_type j=0; j<l; j++)
      for (ordinal_type k=0; k<l; k++)
	os << "i = " << i << ", j = " << j << ", k = " << k << ", Dijk =  " 
	   << operator()(i,j,k) << std::endl;
}
