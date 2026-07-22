// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
