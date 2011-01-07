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
void
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
add_term(ordinal_type i, ordinal_type j, ordinal_type k, const value_type& c)
{
  kji_data[k][j][i] = c;
  ikj_data[i][k][j] = c;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
sum_term(ordinal_type i, ordinal_type j, ordinal_type k, const value_type& c)
{
  kji_data[k][j][i] += c;
  ikj_data[i][k][j] += c;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
print(std::ostream& os) const
{
  for (k_iterator k=k_begin(); k!=k_end(); ++k)
    for (kj_iterator j=j_begin(k); j!=j_end(k); ++j)
      for (kji_iterator i=i_begin(j); i!=i_end(j); ++i)
	os << "k = " << index(k) 
	   << ", j = " << index(j) 
	   << ", i = " << index(i) 
	   << ", Cijk = " << value(i) << std::endl;
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
getValue(ordinal_type i, ordinal_type j, ordinal_type k) const
{
  k_iterator k_it = kji_data.find(k);
  if (k_it == kji_data.end())
    return value_type(0);

  kj_iterator j_it = k_it->second.find(j);
  if (j_it == k_it->second.end())
    return value_type(0);

  kji_iterator i_it = j_it->second.find(i);
  if (i_it == j_it->second.end())
    return value_type(0);

  return i_it->second;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
num_entries() const
{
  ordinal_type num = 0;
  for (k_iterator k = k_begin(); k != k_end(); ++k)
    for (kj_iterator j = j_begin(k); j != j_end(k); ++j)
      for (kji_iterator i = i_begin(j); i != i_end(j); ++i)
	++num;
  return num;
}
