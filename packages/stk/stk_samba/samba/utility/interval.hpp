#ifndef SAMBA_SAMBA_UTILITY_INTERVAL_HPP
#define SAMBA_SAMBA_UTILITY_INTERVAL_HPP

#include <samba/utility/is_primitive.hpp>
#include <samba/utility/value_of.hpp>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#ifdef SAMBA_ENABLE_PARALLEL
#include <boost/serialization/utility.hpp>
#endif

#include <functional>
#include <algorithm>


namespace samba {

/**
 * a right-open interval of primitives
 */
template <typename Primitive>
struct interval
{
  BOOST_MPL_ASSERT((is_primitive<Primitive>));

  typedef Primitive value_type;
  typedef typename value_of<Primitive>::type integer_type;

  BOOST_MPL_ASSERT((boost::is_integral<integer_type>));

  typedef interval<value_type> interval_type;
  typedef size_t size_type;

  struct to_value
  {
    typedef value_type result_type;
    value_type operator()(integer_type i) const
    {
      value_type v = value_type::create(i);
      return v;
    }
  };

  typedef boost::transform_iterator<
            to_value,
            boost::counting_iterator<integer_type>
          > iterator;

  typedef iterator const_iterator;


  interval()
    : m_value(0,0)
  {}

  // interval is empty if lower == upper
  interval( value_type f, value_type p)
    : m_value( f(), p() )
  {}

  // interval is empty if lower == upper
  interval( integer_type f, integer_type p)
    : m_value(f,p)
  {}

  value_type lower() const
  {
    value_type v = value_type::create(m_value.first);
    return v;
  }

  value_type upper() const
  {
    value_type v = value_type::create(m_value.second);
    return v;
  }

  template <typename IntegerType>
  value_type operator[](IntegerType i) const
  { value_type v = value_type::create(static_cast<integer_type>(m_value.first+i));
    return v;
  }

  // interval is empty if lower == upper
  bool empty() const
  { return lower() == upper(); }

  size_type size() const
  { return static_cast<size_type>( m_value.second-m_value.first); }

  const_iterator begin() const
  { return const_iterator( m_value.first, to_value() ); }

  const_iterator end() const
  { return const_iterator( m_value.second, to_value() ); }

  // comparsion operators
  // use lexigraphical ordering
  bool operator < (interval_type rhs) const
  { return m_value < rhs.m_value; }
  bool operator <= (interval_type rhs) const
  { return m_value <= rhs.m_value; }

  bool operator > (interval_type rhs) const
  { return m_value > rhs.m_value; }
  bool operator >= (interval_type rhs) const
  { return m_value >= rhs.m_value; }

  bool operator == (interval_type rhs) const
  { return m_value == rhs.m_value; }
  bool operator != (interval_type rhs) const
  { return m_value != rhs.m_value; }

  private:
#ifdef SAMBA_ENABLE_PARALLEL
    friend class boost::serialization::access;
#endif

    template <typename Archive>
    void serialize (Archive & ar, const unsigned version)
    { ar & m_value; }

    std::pair<integer_type,integer_type> m_value;

};

//*****************************************************************************
// Output operator
//*****************************************************************************
template <typename Primitive>
inline std::ostream & print_interval_range (std::ostream& out, interval<Primitive> const& d)
{
  return out << '[' << d.lower()() << ',' <<  d.upper()() << ")";
}

// std::ostream&  << d
// {Tag(): [d.lower(),d.upper()) }
template <typename Primitive>
inline std::ostream & operator << (std::ostream& out, interval<Primitive> const& d)
{
  out << '{' << typename tag_of<Primitive>::type() << ":";
  print_interval_range(out,d);
  out << "}";
  return out;
}

} // namespace samba

#endif  //SAMBA_SAMBA_UTILITY_INTERVAL_HPP

