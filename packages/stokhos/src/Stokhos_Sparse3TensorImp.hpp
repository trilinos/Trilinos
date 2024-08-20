// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_ConfigDefs.h"
#include "Teuchos_Assert.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
Sparse3Tensor() :
  fill_completed(false)
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
add_term(ordinal_type i, ordinal_type j, ordinal_type k, const value_type& c)
{
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == true, std::logic_error,
		     "You can't call add_term() after calling fillComplete()!");
#endif

  kji_data[k][j][i] = c;
  ikj_data[i][k][j] = c;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
sum_term(ordinal_type i, ordinal_type j, ordinal_type k, const value_type& c)
{
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == true, std::logic_error,
		     "You can't call sum_term() after calling fillComplete()!");
#endif

  kji_data[k][j][i] += c;
  ikj_data[i][k][j] += c;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
fillComplete()
{
  if (fill_completed)
    return;

  ikj_array.resize(ikj_data.size());
  ordinal_type n = 0;
  for (typename ikj_map::const_iterator i_it = ikj_data.begin(); 
       i_it != ikj_data.end(); ++i_it) {
    ikj_array.indices[n] = i_it->first;
    ikj_array.values[n].resize(i_it->second.size());
    ordinal_type l = 0;
    for (typename kj_map::const_iterator k_it = i_it->second.begin(); 
	 k_it != i_it->second.end(); ++k_it) {
      ikj_array.values[n].indices[l] = k_it->first;
      ikj_array.values[n].values[l].resize(k_it->second.size());
      ordinal_type m = 0;
      for (typename j_map::const_iterator j_it = k_it->second.begin(); 
	   j_it != k_it->second.end(); ++j_it) {
	ikj_array.values[n].values[l].indices[m] = j_it->first;
	ikj_array.values[n].values[l].values[m] = j_it->second;
	m++;
      }
      l++;
    }
    n++;
  }

  kji_array.resize(kji_data.size());
  n = 0;
  for (typename kji_map::const_iterator k_it = kji_data.begin(); 
       k_it != kji_data.end(); ++k_it) {
    kji_array.indices[n] = k_it->first;
    kji_array.values[n].resize(k_it->second.size());
    ordinal_type l = 0;
    for (typename ji_map::const_iterator j_it = k_it->second.begin(); 
	 j_it != k_it->second.end(); ++j_it) {
      kji_array.values[n].indices[l] = j_it->first;
      kji_array.values[n].values[l].resize(j_it->second.size());
      ordinal_type m = 0;
      for (typename i_map::const_iterator i_it = j_it->second.begin(); 
	   i_it != j_it->second.end(); ++i_it) {
	kji_array.values[n].values[l].indices[m] = i_it->first;
	kji_array.values[n].values[l].values[m] = i_it->second;
	m++;
      }
      l++;
    }
    n++;
  }

  ikj_data.clear();
  kji_data.clear();

  fill_completed = true;
}

template <typename ordinal_type, typename value_type>
bool
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
fillCompleted() const
{
  return fill_completed;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
print(std::ostream& os) const
{
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling print()!");
#endif

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
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling getValue()!");
#endif

  k_iterator k_it = find_k(k);
  if (k_it == k_end())
    return value_type(0);

  kj_iterator j_it = find_j(k_it, j);
  if (j_it == j_end(k_it))
    return value_type(0);

  kji_iterator i_it = find_i(j_it, i);
  if (i_it == i_end(j_it))
    return value_type(0);

  return i_it.value();
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
num_entries() const
{
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling num_entries()!");
#endif

  ordinal_type num = 0;
  for (k_iterator k = k_begin(); k != k_end(); ++k)
    for (kj_iterator j = j_begin(k); j != j_end(k); ++j)
      for (kji_iterator i = i_begin(j); i != i_end(j); ++i)
	++num;

  return num;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
num_k() const 
{
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling num_k()!");
#endif
  return kji_array.size(); 
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
num_j(const k_iterator& k) const 
{
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling num_j()!");
#endif

  return k.value().size(); 
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
num_i(const kj_iterator& j) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling num_i()!");
#endif
  return j.value().size(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::k_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
find_k(ordinal_type k) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling find_k()!");
#endif
  return kji_array.find(k); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::kj_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
find_j(const k_iterator& k, ordinal_type j) const 
{
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling find_j()!");
#endif
  return k.value().find(j);
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::kji_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
find_i(const kj_iterator& j, ordinal_type i) const 
{
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling find_i()!");
#endif
  return j.value().find(i);
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::k_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
k_begin() const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling k_begin()!");
#endif
  return kji_array.begin(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::k_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
k_end() const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling k_end()!");
#endif
  return kji_array.end(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::k_reverse_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
k_rbegin() const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling k_rbegin()!");
#endif
  return kji_array.rbegin(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::k_reverse_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
k_rend() const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling k_rend()!");
#endif
  return kji_array.rend(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::kj_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
j_begin(const k_iterator& k) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling j_begin()!");
#endif
  return k.value().begin(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::kj_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
j_end(const k_iterator& k) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling j_end()!");
#endif
  return k.value().end(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::kj_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
j_begin(const k_reverse_iterator& k) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling j_begin()!");
#endif
  return k.value().begin(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::kj_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
j_end(const k_reverse_iterator& k) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling j_end()!");
#endif
  return k.value().end(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::kji_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
i_begin(const kj_iterator& j) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling i_begin()!");
#endif
  return j.value().begin(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::kji_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
i_end(const kj_iterator& j) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling i_end()!");
#endif
  return j.value().end(); 
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
num_i() const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling num_i()!");
#endif
  return ikj_array.size(); 
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
num_k(const i_iterator& i) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling num_k()!");
#endif
  return i.value().size(); 
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
num_j(const ik_iterator& k) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling num_j()!");
#endif
  return k.value().size(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::i_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
find_i(ordinal_type i) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling find_i()!");
#endif
  return ikj_array.find(i);
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::ik_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
find_k(const i_iterator& i, ordinal_type k) const 
{
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling find_k()!");
#endif
  return i.value().find(k);
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::ikj_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
find_j(const ik_iterator& k, ordinal_type j) const 
{
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling find_j()!");
#endif
  return k.value().find(j);
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::i_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
i_begin() const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling i_begin()!");
#endif
  return ikj_array.begin(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::i_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
i_end() const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling i_end()!");
#endif
  return ikj_array.end(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::i_reverse_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
i_rbegin() const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling i_rbegin()!");
#endif
  return ikj_array.rbegin();
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::i_reverse_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
i_rend() const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling i_rend()!");
#endif
  return ikj_array.rend(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::ik_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
k_begin(const i_iterator& i) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling k_begin()!");
#endif
  return i.value().begin(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::ik_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
k_end(const i_iterator& i) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling k_end()!");
#endif
  return i.value().end(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::ik_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
k_begin(const i_reverse_iterator& i) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling k_begin()!");
#endif
  return i.value().begin(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::ik_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
k_end(const i_reverse_iterator& i) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling k_end()!");
#endif
  return i.value().end(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::ikj_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
j_begin(const ik_iterator& k) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling j_begin()!");
#endif
  return k.value().begin(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::Sparse3Tensor<ordinal_type, value_type>::ikj_iterator
Stokhos::Sparse3Tensor<ordinal_type, value_type>::
j_end(const ik_iterator& k) const 
{ 
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(fill_completed == false, std::logic_error,
		     "You must call fillComplete() before calling j_end()!");
#endif
  return k.value().end(); 
}
