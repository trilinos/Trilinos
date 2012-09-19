#ifndef STK_UTIL_STK_UTIL_UTIL_NESTED_ITERATOR_HPP
#define STK_UTIL_STK_UTIL_UTIL_NESTED_ITERATOR_HPP

#include <boost/range.hpp>
#include <boost/iterator.hpp>
#include <boost/optional.hpp>

#include <boost/mpl/assert.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace stk {
namespace util {

/** nested_iterator is a forward iterator that flattens iteration over ranges of ranges
 *
 * Both OuterIterator and InnerIterator must support forward_traversal.
 *
 * iterator_value<OuterIterator>::type is convertible to
 * OuterValueToInnerRange::value_type
 *
 */
template <typename OuterRange,
          typename InnerRange,
          typename OuterToInnerConverter
         >
class nested_iterator
  : public boost::iterator_facade<
        nested_iterator<OuterRange,InnerRange,OuterToInnerConverter>
      , typename boost::range_value<InnerRange>::type //value
      , boost::forward_traversal_tag // traversal tag
    >
{
  public:
    typedef typename boost::range_iterator<OuterRange>::type outer_iterator;
    typedef typename boost::range_iterator<InnerRange>::type inner_iterator;

    typedef OuterRange outer_range_type;
    typedef InnerRange inner_range_type;

    typedef OuterToInnerConverter converter_type;

  private:
    typedef nested_iterator<outer_range_type,inner_range_type,converter_type> self;

    typedef typename  boost::range_value<outer_range_type>::type outer_value;


  public:
    nested_iterator()
      : m_outer_current()
      , m_outer_end()
      , m_inner_current()
      , m_inner_end()
      , m_converter()
    {}

    nested_iterator( outer_range_type & outer_range, converter_type converter = converter_type() )
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

    friend class nested_iterator<typename boost::add_const<outer_range_type>::type,
                                 inner_range_type,converter_type>;


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
      while ( (m_outer_current != m_outer_end) && boost::empty((*m_converter)(**m_outer_current)) )
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

    boost::optional<converter_type> m_converter;
};

template <typename OuterRange,
          typename InnerRange,
          typename OuterToInnerConverter
         >
class nested_iterator<const OuterRange, InnerRange,OuterToInnerConverter>
  : public boost::iterator_facade<
        nested_iterator<const OuterRange,InnerRange,OuterToInnerConverter>
      , typename boost::add_const<typename boost::range_value<const InnerRange>::type>::type //value
      , boost::forward_traversal_tag // traversal tag
    >
{
  public:
    typedef typename boost::range_iterator<const OuterRange>::type outer_iterator;
    typedef typename boost::range_iterator<const InnerRange>::type inner_iterator;

    typedef const OuterRange outer_range_type;
    typedef InnerRange inner_range_type;

    typedef OuterToInnerConverter converter_type;

  private:
    typedef nested_iterator<outer_range_type,inner_range_type,converter_type> self;

    typedef typename  boost::range_value<outer_range_type>::type outer_value;


  public:
    nested_iterator()
      : m_outer_current()
      , m_outer_end()
      , m_inner_current()
      , m_inner_end()
      , m_converter()
    {}

    nested_iterator( outer_range_type & outer_range, converter_type converter = converter_type() )
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

    nested_iterator( nested_iterator<typename boost::remove_const<outer_range_type>::type,
                                     inner_range_type,converter_type> const & itr)
      : m_outer_current(itr.m_outer_current)
      , m_outer_end(itr.m_outer_end)
      , m_inner_current(itr.m_inner_current)
      , m_inner_end(itr.m_inner_end)
      , m_converter(itr.m_converter)
    {}
   
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
      while ( (m_outer_current != m_outer_end) && boost::empty((*m_converter)(**m_outer_current)) )
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

    boost::optional<converter_type> m_converter;
};


} // util
} // stk


#endif //STK_UTIL_STK_UTIL_UTIL_NESTED_ITERATOR_HPP
