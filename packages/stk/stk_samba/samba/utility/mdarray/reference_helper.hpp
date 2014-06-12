#ifndef SAMBA_SAMBA_UTILITY_MDARRAY_REFERENCE_HELPER_HPP
#define SAMBA_SAMBA_UTILITY_MDARRAY_REFERENCE_HELPER_HPP

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace samba { namespace detail {

template <typename Array, typename Enable = void>
struct reference_helper;

//TODO -- implement iterators for reference_helpers
template <typename Array, typename Enable = void>
struct reference_iterator {};

template <typename Array>
struct reference_helper< Array
                        ,typename boost::enable_if_c< (boost::rank<Array>::value > 1u) >::type
                       >
{
  typedef Array array;
  typedef typename boost::remove_extent<array>::type sub_array;

  typedef size_t size_type;

  static const size_type static_rank = boost::rank<array>::value;
  static const size_type static_stride = mdarray_length<sub_array>::value;
  static const size_type static_dimension = boost::extent<array>::value;

  typedef typename boost::remove_all_extents<Array>::type data_type;
  typedef reference_helper<sub_array> value_type;

  typedef reference_helper<sub_array> reference;

  static const reference_helper create( data_type * const arg_data, size_type arg_offset=0u)
  { reference_helper ref = {arg_data,arg_offset}; return ref; }

  size_type rank() const { return static_rank; }
  size_type stride() const { return static_stride; }
  size_type dimension() const { return static_dimension; }
  size_type size() const { return dimension(); }

  reference operator[](size_type i) const
  { reference v={m_data,i*static_stride + m_offset}; return v; }

  //data members
  data_type * const m_data;
  const size_type    m_offset;
};

template <typename Array>
struct reference_helper< Array
                        ,typename boost::enable_if_c< (boost::rank<Array>::value == 1u) >::type
                       >
{
  typedef Array array;
  typedef void sub_array;

  typedef size_t size_type;

  static const size_type static_rank = boost::rank<array>::value;
  static const size_type static_stride = 1u;
  static const size_type static_dimension = boost::extent<array>::value;

  typedef typename boost::remove_all_extents<Array>::type data_type;
  typedef data_type value_type;

  typedef value_type & reference;

  static const reference_helper create( data_type * const arg_data, size_type arg_offset=0u)
  { reference_helper ref = {arg_data,arg_offset}; return ref; }

  size_type rank() const { return static_rank; }
  size_type stride() const { return static_stride; }
  size_type dimension() const { return static_dimension; }
  size_type size() const { return dimension(); }

  reference operator[](size_type i) const
  { return m_data[i+m_offset]; }


  //data members
  value_type * const m_data;
  const size_type    m_offset;
};

template <typename T>
struct reference_helper< T
                        ,typename boost::disable_if< boost::is_array<T> >::type
                       >
{
  typedef void array;
  typedef void sub_array;

  typedef size_t size_type;

  static const size_type static_rank = 0;
  static const size_type static_stride = 0;
  static const size_type static_dimension = 0;

  typedef T data_type;
  typedef data_type value_type;

  typedef value_type & reference;

  static const reference_helper create( data_type * const arg_data, size_type /*arg_offset*/)
  { reference_helper ref = {arg_data}; return ref; }

  size_type rank() const { return 0; }
  size_type stride() const { return 0; }
  size_type dimension() const { return 0; }
  size_type size() const { return 1; }

  reference operator[](size_type i) const
  { return *m_data; }

  //data members
  value_type * const m_data;
};

}} //namespace samba::detail


#endif //SAMBA_SAMBA_UTILITY_MDARRAY_REFERENCE_HELPER_HPP

