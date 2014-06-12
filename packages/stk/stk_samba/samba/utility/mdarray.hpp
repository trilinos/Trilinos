#ifndef SAMBA_SAMBA_UTILITY_MDARRAY_HPP
#define SAMBA_SAMBA_UTILITY_MDARRAY_HPP

#include <samba/utility/mdarray/reverse_array.hpp>
#include <samba/utility/mdarray/mdarray_length.hpp>
#include <samba/utility/mdarray/dimension_helper.hpp>
#include <samba/utility/mdarray/reference_helper.hpp>


#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/utility/enable_if.hpp>

#include <algorithm>

namespace samba {


template <typename Array, typename Enable = void>
struct mdarray;


template <typename Array>
struct mdarray<Array, typename boost::enable_if< boost::is_array<Array> >::type >
{
  typedef Array array;
  typedef typename boost::remove_extent<Array>::type sub_array;

  typedef size_t size_type;

  typedef detail::reference_helper<typename boost::remove_cv<array>::type> reference_helper;
  typedef detail::reference_helper<typename boost::add_const<array>::type> const_reference_helper;

  static const size_type static_rank = reference_helper::static_rank;
  static const size_type static_stride = reference_helper::static_stride;
  static const size_type static_dimension = reference_helper::static_dimension;
  static const size_type static_length = detail::mdarray_length<array>::value;

  typedef typename boost::remove_cv<typename reference_helper::data_type>::type data_type;
  typedef typename reference_helper::value_type value_type;

  typedef typename reference_helper::reference reference;
  typedef typename const_reference_helper::reference const_reference;

  size_type rank() const { return static_rank; }
  size_type stride() const { return static_stride; }
  size_type dimension() const { return static_dimension; }

  //size() == dimension so that once we have iterators
  //the array forms a range
  size_type size() const { return dimension(); }

  //how many items of type data_type are allocated
  size_type length() const { return static_length; }

  //prefered fuction to get the dimension
  //allows loops to unroll
  template <size_t Dimension>
  size_type get_dimension() const
  { return boost::extent<array,Dimension>::value;  }

  //convience fuction that only works up to rank 10 arrays
  //loops may not unroll
  size_type dimension(size_type dim) const
  { return detail::dimension_helper<array>(dim); }

  reference operator[](size_type i)
  { return reference_helper::create(m_data,0)[i]; }

  const_reference operator[](size_type i) const
  { return const_reference_helper::create(m_data,0)[i]; }

  mdarray & operator=(const mdarray & b)
  {
    std::copy(b.m_data, b.m_data + static_length, m_data);
    return *this;
  }

  void swap(mdarray &b)
  { std::swap_ranges(b.m_data, b.m_data + static_length, m_data); }

  void assign(data_type const& value)
  { std::fill_n(m_data, static_length, value); }

  //data member
  data_type m_data[static_length];
};


template <typename Scalar>
struct mdarray<Scalar, typename boost::disable_if< boost::is_array<Scalar> >::type >
{
  typedef void array;
  typedef void sub_array;

  typedef size_t size_type;

  static const size_type static_rank = 0u;
  static const size_type static_stride = 0u;
  static const size_type static_dimension = 0u;
  static const size_type static_length = 1u;

  typedef typename boost::remove_cv<typename boost::remove_reference<Scalar>::type>::type data_type;
  typedef data_type value_type;

  typedef value_type &  reference;
  typedef value_type const&  const_reference;

  size_type rank() const { return static_rank; }
  size_type stride() const { return static_stride; }
  size_type dimension() const { return static_dimension; }

  size_type size() const { return 1u; }

  //how many items of type data_type are allocated
  size_type length() const { return static_length; }

  reference operator[](size_type i)
  { return m_data; }

  const_reference operator[](size_type i) const
  { return m_data; }

  mdarray & operator=(const mdarray & b)
  {
    m_data = b.m_data;
    return *this;
  }

  void swap(mdarray &b)
  {
    std::swap(m_data,b.m_data);
  }

  void assign(data_type const& value)
  { m_data = value; }

  //data member
  data_type m_data;
};

} //namespace samba

#endif //SAMBA_SAMBA_UTILITY_MDARRAY_HPP

