#ifndef SIERRA_UTIL_SORTED_VECTOR_HPP
#define SIERRA_UTIL_SORTED_VECTOR_HPP

#include <vector>
#include <algorithm>
#include <functional>

namespace sierra {
namespace util {

struct uniqueS {};
struct not_uniqueS {};

namespace details {

template <class Uniqueness, class Iterator>
struct unique_traits {};

template <class Iterator>
struct unique_traits<uniqueS,Iterator>
{ typedef std::pair<Iterator,bool> type; };

template <class Iterator>
struct unique_traits<not_uniqueS,Iterator>
{ typedef Iterator type; };

} // namespace details


template <
  class ValueType,
  class Compare = std::less<ValueType>,
  class Uniqueness = not_uniqueS,
  class Allocator = std::allocator<ValueType>
  >
class sorted_vector
{
  typedef sorted_vector<ValueType, Compare, Allocator>  self;
  typedef std::vector< ValueType, Allocator >           base_vector;

  struct equal_comp : std::binary_function<ValueType,ValueType,bool> {
    const Compare &m_comp;
    equal_comp( const Compare & comp ) : m_comp(comp) {}
    bool operator() ( const ValueType & lhs, const ValueType & rhs) const {
      return !m_comp(lhs,rhs) && !m_comp(rhs,lhs);
    }
  };

  public:
    // concepts from std::vector
    typedef typename base_vector::const_reference        const_reference;
    typedef typename base_vector::const_reference        reference;

    typedef typename base_vector::const_pointer          const_pointer;
    typedef typename base_vector::const_pointer          pointer;

    typedef typename base_vector::const_iterator         const_iterator;
    typedef typename base_vector::const_iterator         iterator;

    typedef typename base_vector::const_reverse_iterator const_reverse_iterator;
    typedef typename base_vector::const_reverse_iterator reverse_iterator;

    typedef typename base_vector::size_type              size_type;
    typedef typename base_vector::difference_type        difference_type;

    typedef typename base_vector::value_type             value_type;

    typedef typename base_vector::allocator_type         allocator_type;

    typedef Compare                                      value_compare;
    typedef Uniqueness                                   unique_type;


    sorted_vector( const Compare & comp = Compare(), const Allocator & alloc = Allocator())
      : m_vector(alloc)
      , m_comp(comp)
    {}

    sorted_vector(size_t n, const value_type & val = value_type(), const Compare & comp = Compare(), const Allocator & alloc = Allocator())
      : m_vector(n,val,alloc)
      , m_comp(comp)
    {}

    template <class InputIterator>
    sorted_vector( InputIterator first, InputIterator last, const Compare & comp, const Allocator & alloc = Allocator())
      : m_vector()
      , m_comp(comp)
    {
      insert(first,last);
    }

    sorted_vector( const self & rhs)
    : m_vector(rhs.m_vector)
    , m_comp(rhs.m_comp)
    {}

    sorted_vector operator = ( const self & rhs)
    {
      if ( this != &rhs ) {
        self temp(rhs);
        swap(temp);
      }
      return *this;
    }

    const_iterator begin() const { return m_vector.begin(); }
    const_iterator end()   const { return m_vector.end(); }

    const_reverse_iterator rbegin() const { return m_vector.rbegin(); }
    const_reverse_iterator rend()   const { return m_vector.rend(); }

    size_type size()     const { return m_vector.size(); }
    size_type max_size() const { return m_vector.max_size(); }
    size_type capacity() const { return m_vector.capacity(); }

    bool empty() const { return m_vector.empty(); }

    void reserve(size_type n) { m_vector.reserve(n); }

    const_reference operator []( size_type n ) const { return m_vector[n]; }
    const_reference at( size_type n ) const { return m_vector.at(n); }

    const_reference front() const { return m_vector.front(); }
    const_reference back() const { return m_vector.back(); }

    void swap ( self & rhs ) {
      m_vector.swap(rhs.m_vector);
      std::swap(m_comp,rhs.m_comp);
    }

    void clear() { m_vector.clear(); }

    allocator_type get_allocator() const { return m_vector.get_allocator(); }

    iterator insert( const value_type &x)
    {
      return insert( x, unique_type() );
    }

    iterator insert( iterator position, const value_type &x)
    {
      if ( size() == 0 ) {
        return insert(x, unique_type());
      }
      return insert( position, x, unique_type() );
    }

    template <class InputIterator>
    void insert( InputIterator first, InputIterator last)
    {
      insert(first,last,unique_type());
    }

    iterator erase ( iterator pos )
    {
      difference_type n = std::distance( begin(), pos );
      typename base_vector::iterator position = m_vector.begin() + n;
      return m_vector.erase(position);
    }

    iterator erase ( const value_type & x )
    {
      typename base_vector::iterator b = std::lower_bound(m_vector.begin(),m_vector.end(),x, m_comp );
      typename base_vector::iterator e = std::upper_bound(b,m_vector.end(), x, m_comp );
      return   m_vector.erase(b,e);
    }

    iterator erase( iterator first, iterator last)
    {
      difference_type n = std::distance( begin(), first );
      typename base_vector::iterator b = m_vector.begin() + n;

      n = std::distance( begin(), last );
      typename base_vector::iterator e = m_vector.begin() + n;

      return m_vector.erase( b, e );
    }

    void pop_back() { m_vector.pop_back(); }

   // allow serialization
  private:

    //----Begin Insert Helpers-----
    iterator insert( const value_type &x, const uniqueS &)
    {
      typename base_vector::iterator pos = std::lower_bound(m_vector.begin(),m_vector.end(),x, m_comp );
      //bool inserted = false;
      equal_comp eq_comp(m_comp);
      if ( pos == m_vector.end() || !eq_comp(*pos,x) ) {
        //inserted = true;
        pos = m_vector.insert(pos,x);
      }
      //return std::make_pair(iterator(pos),inserted);
      return pos;
    }

    iterator insert( const value_type &x, const not_uniqueS &)
    {
      typename base_vector::iterator pos = std::lower_bound(m_vector.begin(),m_vector.end(),x, m_comp );
      return iterator(m_vector.insert(pos,x));
    }

    template <class InputIterator>
    void insert( InputIterator first, InputIterator last, const uniqueS &)
    {
      size_type n = m_vector.size();
      m_vector.insert(m_vector.end(), first, last);

      typename base_vector::iterator pos = m_vector.begin()+n;
      std::sort(pos,m_vector.end(),m_comp);
      std::inplace_merge(m_vector.begin(),pos,m_vector.end(),m_comp);

      pos = std::unique(m_vector.begin(),m_vector.end(),equal_comp(m_comp));
      m_vector.resize(pos-m_vector.begin());
    }

    template <class InputIterator>
    void insert( InputIterator first, InputIterator last, const not_uniqueS &)
    {
      size_type n = m_vector.size();
      m_vector.insert(m_vector.end(), first, last);

      typename base_vector::iterator pos = m_vector.begin()+n;
      std::sort(pos,m_vector.end(),m_comp);
      std::inplace_merge(m_vector.begin(),pos,m_vector.end(),m_comp);
    }

    iterator insert( iterator position, const value_type &x, const uniqueS &)
    {
      difference_type n = std::distance( begin(), position );
      typename base_vector::iterator pos = m_vector.begin() + n;
      equal_comp eq_comp(m_comp);

      if ( pos != m_vector.end() && m_comp(*(pos-1),x) && m_comp(x,*pos) ) { // correct position
        pos = m_vector.insert(pos,x);
        //return std::make_pair(iterator(pos),true);
        return pos;
      }
      else if ( pos == m_vector.end() && m_comp(*(pos-1),x) ) { // correct but at end
        pos = m_vector.insert(pos,x);
        //return std::make_pair(iterator(pos),true);
        return pos;
      }
      else if ( pos != m_vector.end() && !eq_comp(*pos,x) ) { // correct but not unique
        //return std::make_pair(position,false);
        return pos;
      }

      //incorrect position--do regular insert
      return insert(x, unique_type());
    }

    iterator insert( iterator position, const value_type &x, const not_uniqueS &)
    {
      difference_type n = std::distance( begin(), position );
      typename base_vector::iterator pos = m_vector.begin() + n;

      if ( pos != m_vector.end() && m_comp(*(pos-1),x) && m_comp(x,*pos) ) { // correct position
        pos = m_vector.insert(pos,x);
        //return std::make_pair(iterator(pos),true);
        return pos;
      }
      else if ( pos == m_vector.end() && m_comp(*(pos-1),x) ) { // correct but at end
        pos = m_vector.insert(pos,x);
        //return std::make_pair(iterator(pos),true);
        return pos;
      }

      //incorrect position--do regular insert
      return insert(x, unique_type());
    }

    //----End Insert Helpers-----


    base_vector   m_vector;
    value_compare m_comp;

};

} // namespace util
} // namespace sierra

#endif //SIERRA_UTIL_SORTED_VECTOR_HPP
