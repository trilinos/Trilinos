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
