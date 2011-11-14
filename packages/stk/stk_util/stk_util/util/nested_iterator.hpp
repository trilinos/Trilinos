#ifndef STK_UTIL_STK_UTIL_UTIL_NESTED_ITERATOR_HPP
#define STK_UTIL_STK_UTIL_UTIL_NESTED_ITERATOR_HPP

#include <boost/range.hpp>
#include <boost/iterator.hpp>
#include <boost/optional.hpp>

#include <boost/mpl/assert.hpp>
#include <boost/type_traits.hpp>

namespace stk {
namespace util {

/** nested_iterator is a forward iterator that flattens iteration over ranges of ranges
 *
 * Both OuterIterator and InnerIterator must support forward_traversal.
 *
 * OuterValueToInnerRange
 *
 * struct OuterValueToInnerRange {
 *   typedef ... value_type;
 *   typedef ... result_type;
 *
 *   result_type operator()(const value_type & value) {
 *     return ...;
 *   }
 * };
 *
 * iterator_value<OuterIterator>::type is convertible to
 * OuterValueToInnerRange::value_type
 *
 * OuterValueToInnerRange::result_type is convertible to
 * boost::iterator_range<InnerIterator>
 */
template <typename OuterIterator,
          typename InnerIterator,
          typename OuterValueToInnerRange
         >
class nested_iterator
  : public boost::iterator_facade<
        nested_iterator<OuterIterator,InnerIterator,OuterValueToInnerRange>
      , typename boost::iterator_value<InnerIterator>::type //value
      , boost::forward_traversal_tag // traversal tag
    >
{
  public:
    typedef OuterIterator outer_iterator;
    typedef InnerIterator inner_iterator;
    typedef OuterValueToInnerRange outer_value_to_inner_range;

    typedef boost::iterator_range<outer_iterator> outer_range_type;
    typedef boost::iterator_range<inner_iterator> inner_range_type;

  private:
    typedef nested_iterator<outer_iterator,inner_iterator,outer_value_to_inner_range> self;

    typedef typename  boost::iterator_value<outer_iterator>::type outer_value;
    typedef typename outer_value_to_inner_range::value_type converter_value;
    typedef typename outer_value_to_inner_range::result_type converter_result;

    BOOST_MPL_ASSERT_MSG(
        ( boost::is_convertible<
            outer_value,
            converter_value
          >::value
        ),
        NO_CONVERSION_FOUND,
        (
         outer_value,
         converter_value
        )
    );

    BOOST_MPL_ASSERT_MSG(
        ( boost::is_convertible<
            converter_result,
            inner_range_type
          >::value
        ),
        NO_CONVERSION_FOUND,
        ( converter_result,
          inner_range_type
        )
    );


  public:
    nested_iterator()
      : m_outer_current()
      , m_outer_end()
      , m_inner_current()
      , m_inner_end()
      , m_converter()
    {}

    nested_iterator( const outer_range_type & outer_range, const outer_value_to_inner_range converter )
      : m_outer_current(boost::begin(outer_range))
      , m_outer_end(boost::end(outer_range))
      , m_inner_current()
      , m_inner_end()
      , m_converter(converter)
    {
      if ( boost::empty(outer_range) ) {
        m_outer_current = boost::none;
        m_outer_end     = boost::none;
      }
      else {
        find_inner_range_helper();
      }


    }

    //
  private:
    friend class boost::iterator_core_access;

    //functions necessary to implement core operations
    // increment
    // equal
    // dereference

    void increment()
    {
      if (m_inner_current) {
        ++(*m_inner_current);
        //at end of current inner range
        if (m_inner_current == m_inner_end) {
          ++(*m_outer_current);
          find_inner_range_helper();
        }
      }
    }

    bool equal( const self & rhs ) const
    {
      return    (m_outer_current == rhs.m_outer_current)
             && (m_inner_current == rhs.m_inner_current);
    }

    typename boost::iterator_reference<inner_iterator>::type dereference() const
    {
      return **m_inner_current;
    }

    void find_inner_range_helper()
    {
      // find the next none-empty inner_range
      while ( (m_outer_current != m_outer_end) &&
               boost::empty((*m_converter)(**m_outer_current)) )
      {
        ++(*m_outer_current);
      }

      if (m_outer_current != m_outer_end) {
        m_inner_current = boost::begin((*m_converter)(**m_outer_current));
        m_inner_end = boost::end((*m_converter)(**m_outer_current));
      }
      else { //at end of outer range
        m_outer_current = boost::none;
        m_outer_end = boost::none;

        m_inner_current = boost::none;
        m_inner_end = boost::none;
      }
    }


  private:
    boost::optional<outer_iterator> m_outer_current;
    boost::optional<outer_iterator> m_outer_end;

    boost::optional<inner_iterator> m_inner_current;
    boost::optional<inner_iterator> m_inner_end;

    boost::optional<outer_value_to_inner_range> m_converter;
};


} // util
} // stk


#endif //STK_UTIL_STK_UTIL_UTIL_NESTED_ITERATOR_HPP
