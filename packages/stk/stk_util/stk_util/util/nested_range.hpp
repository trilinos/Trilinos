#ifndef STK_UTIL_STK_UTIL_UTIL_NESTED_RANGE_HPP
#define STK_UTIL_STK_UTIL_UTIL_NESTED_RANGE_HPP

#include <stk_util/util/nested_iterator.hpp>

#include <boost/range.hpp>
#include <boost/iterator.hpp>
#include <boost/optional.hpp>

#include <boost/mpl/assert.hpp>
#include <boost/type_traits.hpp>

namespace stk {
namespace util {

namespace details {

template<typename T>
struct identity {
  typedef T result_type;

  result_type& operator()(result_type& r) const { return r; }
  const result_type& operator()(const result_type& r) const { return r; }
};

}

/** nested_range is a forward range that flattens iteration over ranges of ranges
 */
template < typename OuterRange,
           typename InnerRange=typename boost::range_value<OuterRange>::type,
           typename OuterToInnerConverter=
             details::identity<
               typename boost::mpl::if_<
                 typename boost::is_same<InnerRange,
                                         typename boost::range_value<OuterRange>::type>,
                 InnerRange,
                 void
               >::type
             >
         >
class nested_range;

template < typename OuterRange, typename InnerRange, typename OuterToInnerConverter >
class nested_range
{
  public:
    typedef OuterRange outer_range;
    typedef InnerRange inner_range;
    BOOST_MPL_ASSERT((boost::has_range_iterator<outer_range>));
    BOOST_MPL_ASSERT((boost::has_range_iterator<inner_range>));

    typedef OuterToInnerConverter converter_type;

    typedef nested_iterator<outer_range,inner_range,converter_type> iterator;
    typedef nested_iterator<typename boost::add_const<outer_range>::type,inner_range,converter_type> const_iterator;

  nested_range() : m_outer(), m_converter() {}

  nested_range(outer_range& outer, converter_type converter=converter_type()) : m_outer(outer), m_converter(converter) {}

  iterator begin() { return iterator(*m_outer, *m_converter); }
  const_iterator begin() const { return const_iterator(*m_outer, *m_converter); }

  iterator end() { return iterator(); }
  const_iterator end() const { return const_iterator(); }

  private:
    boost::optional<outer_range&> m_outer;
    boost::optional<converter_type> m_converter;
};

//template < typename OuterRange>
//class nested_range<OuterRange,typename boost::range_value<OuterRange>::type, void>
//{
//  public:
//    typedef OuterRange outer_range;
//    typedef InnerRange inner_range;
//    BOOST_MPL_ASSERT((boost::has_range_iterator<inner_range>));
//
//    typedef nested_iterator<outer_range,inner_range> iterator;
//    typedef nested_iterator<typename boost::add_const<outer_range>::type,inner_range> const_iterator;
//
//  nested_range() : m_outer() {}
//
//  nested_range(outer_range& outer) : m_outer(outer) {}
//
//  iterator begin() { return iterator(*m_outer); }
//  const_iterator begin() const { return const_iterator(*m_outer); }
//
//  iterator end() { return iterator(); }
//  const_iterator end() const { return const_iterator(); }
//
//  private:
//    boost::optional<outer_range&> m_outer;
//};


} // util
} // stk


#endif //STK_UTIL_STK_UTIL_UTIL_NESTED_RANGE_HPP
