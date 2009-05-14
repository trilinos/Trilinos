// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
Sparse3Tensor(ordinal_type sz) :
  i_indices(sz),
  j_indices(sz),
  Cijk_values(sz)
{
}

template <typename ordinal_type, typename value_type>
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
~Sparse3Tensor()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
size() const
{
  return Cijk_values.size();
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
num_values(ordinal_type k) const
{
  return Cijk_values[k].size();
}

template <typename ordinal_type, typename value_type>
void
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
value(ordinal_type k, ordinal_type l, ordinal_type& i, ordinal_type& j, 
      value_type& c) const
{
  i = i_indices[k][l];
  j = j_indices[k][l];
  c = Cijk_values[k][l];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
add_term(ordinal_type i, ordinal_type j, ordinal_type k, const value_type& c)
{
  i_indices[k].push_back(i);
  j_indices[k].push_back(j);
  Cijk_values[k].push_back(c);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
print(std::ostream& os) const
{
  for (ordinal_type k=0; k<static_cast<ordinal_type>(Cijk_values.size()); k++)
    for (ordinal_type l=0; l<static_cast<ordinal_type>(Cijk_values[k].size()); 
	 l++)
      os << "k = " << k << ", l = " << l 
	 << ", i = " << i_indices[k][l] << ", j = " << j_indices[k][l]
	 << ", Cijk = " << Cijk_values[k][l] << std::endl;
}

