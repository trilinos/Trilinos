#ifndef STK_UTIL_BLOCK_VECTOR_HPP
#define STK_UTIL_BLOCK_VECTOR_HPP

#include <stk_util/util/CacheAlignedAllocator.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <algorithm>

namespace stk {

namespace details {

template <size_t N>
struct power_of_2
{
  static const size_t value = 1 + power_of_2<(N>>1)>::value;
};

template <>
struct power_of_2<2>
{
  static const size_t value = 1;
};

template <>
struct power_of_2<1>
{
  static const size_t value = 0;
};

template <>
struct power_of_2<0>
{
  static const size_t value = 0;
};

} // namespace details

template <typename Type, size_t BlockSize = 512, typename Allocator = stk::cache_aligned_allocator<Type> >
class block_vector
{
public:
  // require that BlockSize be a power of 2
  BOOST_STATIC_ASSERT(( BlockSize != 0u && !( BlockSize & (BlockSize-1u)) ));

  typedef Type         value_type;
  typedef Type       * pointer;
  typedef Type const * const_pointer;
  typedef Type       & reference;
  typedef Type const & const_reference;
  typedef ptrdiff_t    difference_type;
  typedef size_t       size_type;

  typedef Allocator    block_allocator;
  typedef typename block_allocator::template rebind<pointer>::other allocator;

  static const size_type block_size = BlockSize;
  static const size_type mask = BlockSize - 1;
  static const size_type shift = details::power_of_2<BlockSize>::value;

  struct const_iterator
  {
    typedef std::random_access_iterator_tag iterator_category;

    typedef Type         value_type;
    typedef Type       * pointer;
    typedef Type const * const_pointer;
    typedef Type       & reference;
    typedef Type const & const_reference;
    typedef ptrdiff_t    difference_type;
    typedef size_t       size_type;


    //operator *
    const_reference operator*() const
    { return m_blocks[m_offset >> shift][m_offset & mask]; }

    const_pointer operator->() const
    { return m_blocks[m_offset >> shift] + (m_offset & mask); }

    //operator[n]
    template <typename IntType>
    const_reference operator[](IntType n) const
    { return m_blocks[ (m_offset + n) >> shift][ (m_offset + n) & mask]; }

    //operator++
    const_iterator & operator++()
    { ++m_offset; return *this; }

    //operator++(int)
    const_iterator operator++(int)
    { const_iterator tmp(*this); ++m_offset; return tmp; }

    //operator--
    const_iterator & operator--()
    { --m_offset; return *this; }

    //operator--(int)
    const_iterator operator--(int)
    { const_iterator tmp(*this); --m_offset; return tmp; }

    //operator += n
    template <typename IntType>
    const_iterator & operator+=(IntType n)
    { m_offset += n; return *this; }

    //operator -= n
    template <typename IntType>
    const_iterator & operator-=(IntType n)
    { m_offset -= n; return *this; }

    //const_iterator + n
    template <typename IntType>
    friend inline const_iterator operator+(const_iterator tmp, IntType n)
    { return tmp += n; }

    //n + const_iterator
    template <typename IntType>
    friend inline const_iterator operator+(IntType n, const_iterator tmp)
    { return tmp += n; }

    //const_iterator - n
    template <typename IntType>
    friend inline const_iterator operator-(const_iterator tmp, IntType n)
    { return tmp -= n; }

    //difference
    friend inline difference_type operator - (const_iterator const & lhs, const_iterator const & rhs)
    { return lhs.m_offset - rhs.m_offset; }

    // comparisons
    friend inline bool operator == (const_iterator const& lhs, const_iterator const& rhs)
    { return lhs.m_offset == rhs.m_offset; }

    friend inline bool operator != (const_iterator const& lhs, const_iterator const& rhs)
    { return lhs.m_offset != rhs.m_offset; }

    friend inline bool operator < (const_iterator const& lhs, const_iterator const& rhs)
    { return lhs.m_offset < rhs.m_offset; }

    friend inline bool operator <= (const_iterator const& lhs, const_iterator const& rhs)
    { return lhs.m_offset <= rhs.m_offset; }

    friend inline bool operator > (const_iterator const& lhs, const_iterator const& rhs)
    { return lhs.m_offset > rhs.m_offset; }

    friend inline bool operator >= (const_iterator const& lhs, const_iterator const& rhs)
    { return lhs.m_offset >= rhs.m_offset; }


    //members
    value_type const * const * m_blocks;
    size_type                  m_offset;
  };

  struct iterator
  {
    typedef std::random_access_iterator_tag iterator_category;

    typedef Type         value_type;
    typedef Type       * pointer;
    typedef Type const * const_pointer;
    typedef Type       & reference;
    typedef Type const & const_reference;
    typedef ptrdiff_t    difference_type;
    typedef size_t       size_type;


    //operator *
    reference operator*() const
    { return m_blocks[m_offset >> shift][m_offset & mask]; }

    pointer operator->() const
    { return m_blocks[m_offset >> shift] + (m_offset & mask); }

    //operator[n]
    template <typename IntType>
    reference operator[](IntType n) const
    { return m_blocks[ (m_offset + n) >> shift][ (m_offset + n) & mask]; }

    //operator++
    iterator & operator++()
    { ++m_offset; return *this; }

    //operator++(int)
    iterator operator++(int)
    { iterator tmp(*this); ++m_offset; return tmp; }

    //operator--
    iterator & operator--()
    { --m_offset; return *this; }

    //operator--(int)
    iterator operator--(int)
    { iterator tmp(*this); --m_offset; return tmp; }

    //operator += n
    template <typename IntType>
    iterator & operator+=(IntType n)
    { m_offset += n; return *this; }

    //operator -= n
    template <typename IntType>
    iterator & operator-=(IntType n)
    { m_offset -= n; return *this; }

    //iterator + n
    template <typename IntType>
    friend inline iterator operator+(iterator tmp, IntType n)
    { return tmp += n; }

    //n + iterator
    template <typename IntType>
    friend inline iterator operator+(IntType n, iterator tmp)
    { return tmp += n; }

    //iterator - n
    template <typename IntType>
    friend inline iterator operator-(iterator tmp, IntType n)
    { return tmp -= n; }

    //difference
    friend inline difference_type operator - (iterator const & lhs, iterator const & rhs)
    { return lhs.m_offset - rhs.m_offset; }

    // comparisons
    friend inline bool operator == (iterator const& lhs, iterator const& rhs)
    { return lhs.m_offset == rhs.m_offset; }

    friend inline bool operator != (iterator const& lhs, iterator const& rhs)
    { return lhs.m_offset != rhs.m_offset; }

    friend inline bool operator < (iterator const& lhs, iterator const& rhs)
    { return lhs.m_offset < rhs.m_offset; }

    friend inline bool operator <= (iterator const& lhs, iterator const& rhs)
    { return lhs.m_offset <= rhs.m_offset; }

    friend inline bool operator > (iterator const& lhs, iterator const& rhs)
    { return lhs.m_offset > rhs.m_offset; }

    friend inline bool operator >= (iterator const& lhs, iterator const& rhs)
    { return lhs.m_offset >= rhs.m_offset; }

    operator const_iterator() const
    { const_iterator itr = { m_blocks, m_offset }; return itr; }

    //members
    value_type       ** m_blocks;
    size_type           m_offset;
  };

  struct const_reverse_iterator
  {
    typedef std::random_access_iterator_tag iterator_category;

    typedef Type         value_type;
    typedef Type       * pointer;
    typedef Type const * const_pointer;
    typedef Type       & reference;
    typedef Type const & const_reference;
    typedef ptrdiff_t    difference_type;
    typedef size_t       size_type;


    //operator *
    const_reference operator*() const
    { return m_blocks[ (m_offset-1u) >> shift][ (m_offset-1u) & mask]; }

    const_pointer operator->() const
    { return m_blocks[(m_offset -1u) >> shift] + ((m_offset -1u) & mask); }

    //operator[n]
    template <typename IntType>
    const_reference operator[](IntType n) const
    { return m_blocks[ ((m_offset -1u) + n) >> shift][ ((m_offset -1u) + n) & mask]; }

    //operator++
    const_reverse_iterator & operator++()
    { --m_offset; return *this; }

    //operator++(int)
    const_reverse_iterator operator++(int)
    { const_reverse_iterator tmp(*this); --m_offset; return tmp; }

    //operator--
    const_reverse_iterator & operator--()
    { ++m_offset; return *this; }

    //operator--(int)
    const_reverse_iterator operator--(int)
    { const_reverse_iterator tmp(*this); ++m_offset; return tmp; }

    //operator += n
    template <typename IntType>
    const_reverse_iterator & operator+=(IntType n)
    { m_offset -= n; return *this; }

    //operator -= n
    template <typename IntType>
    const_reverse_iterator & operator-=(IntType n)
    { m_offset += n; return *this; }

    //const_reverse_iterator + n
    template <typename IntType>
    friend inline const_reverse_iterator operator+(const_reverse_iterator tmp, IntType n)
    { return tmp += n; }

    //n + const_reverse_iterator
    template <typename IntType>
    friend inline const_reverse_iterator operator+(IntType n, const_reverse_iterator tmp)
    { return tmp += n; }

    //const_reverse_iterator - n
    template <typename IntType>
    friend inline const_reverse_iterator operator-(const_reverse_iterator tmp, IntType n)
    { return tmp -= n; }

    //difference
    friend inline difference_type operator - (const_reverse_iterator const & lhs, const_reverse_iterator const & rhs)
    { return rhs.m_offset - lhs.m_offset; }

    // comparisons
    friend inline bool operator == (const_reverse_iterator const& lhs, const_reverse_iterator const& rhs)
    { return lhs.m_offset == rhs.m_offset; }

    friend inline bool operator != (const_reverse_iterator const& lhs, const_reverse_iterator const& rhs)
    { return lhs.m_offset != rhs.m_offset; }

    friend inline bool operator < (const_reverse_iterator const& lhs, const_reverse_iterator const& rhs)
    { return lhs.m_offset > rhs.m_offset; }

    friend inline bool operator <= (const_reverse_iterator const& lhs, const_reverse_iterator const& rhs)
    { return lhs.m_offset >= rhs.m_offset; }

    friend inline bool operator > (const_reverse_iterator const& lhs, const_reverse_iterator const& rhs)
    { return lhs.m_offset < rhs.m_offset; }

    friend inline bool operator >= (const_reverse_iterator const& lhs, const_reverse_iterator const& rhs)
    { return lhs.m_offset <= rhs.m_offset; }


    //members
    value_type const * const * m_blocks;
    size_type                  m_offset;
  };

  struct reverse_iterator
  {
    typedef std::random_access_iterator_tag iterator_category;

    typedef Type         value_type;
    typedef Type       * pointer;
    typedef Type const * const_pointer;
    typedef Type       & reference;
    typedef Type const & const_reference;
    typedef ptrdiff_t    difference_type;
    typedef size_t       size_type;


    //operator *
    reference operator*() const
    { return m_blocks[(m_offset -1u) >> shift][(m_offset -1u) & mask]; }

    pointer operator->() const
    { return m_blocks[(m_offset -1u) >> shift] + ((m_offset -1u) & mask); }

    //operator[n]
    template <typename IntType>
    reference operator[](IntType n) const
    { return m_blocks[ ((m_offset -1u) + n) >> shift][ ((m_offset -1u) + n) & mask]; }

    //operator++
    reverse_iterator & operator++()
    { --m_offset; return *this; }

    //operator++(int)
    reverse_iterator operator++(int)
    { reverse_iterator tmp(*this); --m_offset; return tmp; }

    //operator--
    reverse_iterator & operator--()
    { ++m_offset; return *this; }

    //operator--(int)
    reverse_iterator operator--(int)
    { reverse_iterator tmp(*this); ++m_offset; return tmp; }

    //operator += n
    template <typename IntType>
    reverse_iterator & operator+=(IntType n)
    { m_offset -= n; return *this; }

    //operator -= n
    template <typename IntType>
    reverse_iterator & operator-=(IntType n)
    { m_offset += n; return *this; }

    //reverse_iterator + n
    template <typename IntType>
    friend inline reverse_iterator operator+(reverse_iterator tmp, IntType n)
    { return tmp += n; }

    //n + reverse_iterator
    template <typename IntType>
    friend inline reverse_iterator operator+(IntType n, reverse_iterator tmp)
    { return tmp += n; }

    //reverse_iterator - n
    template <typename IntType>
    friend inline reverse_iterator operator-(reverse_iterator tmp, IntType n)
    { return tmp -= n; }

    //difference
    friend inline difference_type operator - (reverse_iterator const & lhs, reverse_iterator const & rhs)
    { return rhs.m_offset - lhs.m_offset; }

    // comparisons
    friend inline bool operator == (reverse_iterator const& lhs, reverse_iterator const& rhs)
    { return lhs.m_offset == rhs.m_offset; }

    friend inline bool operator != (reverse_iterator const& lhs, reverse_iterator const& rhs)
    { return lhs.m_offset != rhs.m_offset; }

    friend inline bool operator < (reverse_iterator const& lhs, reverse_iterator const& rhs)
    { return lhs.m_offset > rhs.m_offset; }

    friend inline bool operator <= (reverse_iterator const& lhs, reverse_iterator const& rhs)
    { return lhs.m_offset >= rhs.m_offset; }

    friend inline bool operator > (reverse_iterator const& lhs, reverse_iterator const& rhs)
    { return lhs.m_offset < rhs.m_offset; }

    friend inline bool operator >= (reverse_iterator const& lhs, reverse_iterator const& rhs)
    { return lhs.m_offset <= rhs.m_offset; }

    operator const_reverse_iterator() const
    { const_reverse_iterator itr = { m_blocks, m_offset }; return itr; }

    //members
    value_type       ** m_blocks;
    size_type           m_offset;
  };


  block_vector()
    : m_size(0)
    , m_blocks_size(0)
    , m_blocks_capacity(0)
    , m_blocks(NULL)
  {}

  block_vector(size_type n, const_reference val = value_type())
    : m_size(0)
    , m_blocks_size(0)
    , m_blocks_capacity(0)
    , m_blocks(NULL)
  {
    allocate(n);

    for (size_type i=0; i<m_size; ++i) {
      pointer buf = m_blocks[i >> shift] + (i & mask);
      block_allocator::construct(buf, val);
    }
  }

  block_vector( const block_vector & rhs )
    : m_size(0)
    , m_blocks_size(0)
    , m_blocks_capacity(0)
    , m_blocks(NULL)
  {
    allocate(rhs.m_size);
    for (size_type i=0; i<m_size; ++i) {
      pointer buf = m_blocks[i >> shift] + (i & mask);
      block_allocator::construct(buf, rhs[i]);
    }
  }

  template <typename InputIterator>
  block_vector( InputIterator first, InputIterator second
               ,typename boost::enable_if_c< !(boost::is_integral<InputIterator>::value || boost::is_enum<InputIterator>::value), int>::type = 0
              )
    : m_size(0)
    , m_blocks_size(0)
    , m_blocks_capacity(0)
    , m_blocks(NULL)
  {
    assign(first, second);
  }

  block_vector & operator=( const block_vector & rhs)
  {
    if (&rhs != this) {
      block_vector tmp(rhs);
      this->swap(tmp);
    }
    return *this;
  }

  ~block_vector()
  {
    //delete objects
    for (size_type i=0; i<m_size; ++i) {
      pointer p = m_blocks[i >> shift] + (i & mask);
      block_allocator::destroy(p);
    }

    //deallocate blocks
    for (size_type block=0; block<m_blocks_size; ++block) {
      if (m_blocks[block]) {
        block_allocator::deallocate(m_blocks[block], block_size);
      }
    }
    if (m_blocks) {
      allocator::deallocate(m_blocks, m_blocks_capacity);
    }
  }


        iterator  begin()       {       iterator itr = {m_blocks,0}; return itr; }
  const_iterator  begin() const { const_iterator itr = {m_blocks,0}; return itr; }
  const_iterator cbegin() const { const_iterator itr = {m_blocks,0}; return itr; }

        iterator  end()       {       iterator itr = {m_blocks,m_size}; return itr; }
  const_iterator  end() const { const_iterator itr = {m_blocks,m_size}; return itr; }
  const_iterator cend() const { const_iterator itr = {m_blocks,m_size}; return itr; }

        reverse_iterator  rbegin()       {       reverse_iterator itr = {m_blocks,m_size}; return itr; }
  const_reverse_iterator  rbegin() const { const_reverse_iterator itr = {m_blocks,m_size}; return itr; }
  const_reverse_iterator crbegin() const { const_reverse_iterator itr = {m_blocks,m_size}; return itr; }

        reverse_iterator  rend()       {       reverse_iterator itr = {m_blocks,0u}; return itr; }
  const_reverse_iterator  rend() const { const_reverse_iterator itr = {m_blocks,0u}; return itr; }
  const_reverse_iterator crend() const { const_reverse_iterator itr = {m_blocks,0u}; return itr; }


  size_type size() const { return m_size; }
  size_type capacity() const { return m_blocks_size * block_size; }
  bool empty() const { return m_size == 0u; };

  void resize(size_type n, const_reference val = value_type())
  {
    if (n > m_size) { //grow array
      const size_t curr_size = m_size;
      allocate(n);
      for (size_type i=curr_size; i<m_size; ++i) {
        pointer buf = m_blocks[i >> shift] + (i & mask);
        block_allocator::construct(buf, val);
      }
    }
    else if (n < m_size) { //shrink array
      for (size_type i=n; i<m_size; ++i) {
        pointer p = m_blocks[i >> shift] + (i & mask);
        block_allocator::destroy(p);
      }
    }

    m_size = n;
  }

  void shrink_to_fit()
  {
    size_type required_blocks = (m_size + block_size - 1) / block_size;

    if ( required_blocks < m_blocks_capacity ) {

      //free unused blocks
      for (size_type block = required_blocks; block < m_blocks_size; ++block) {
        if (m_blocks[block]) {
          block_allocator::deallocate(m_blocks[block], block_size);
          m_blocks[block] = NULL;
        }
      }

      size_type old_capacity = m_blocks_capacity;

      m_blocks_size = required_blocks;
      m_blocks_capacity = required_blocks;

      //copy over remaining blocks blocks
      value_type **tmpdata = allocator::allocate(m_blocks_capacity);

      memmove ( reinterpret_cast<void*>(tmpdata), reinterpret_cast<void*>(m_blocks), m_blocks_size * sizeof(value_type**));
      if (old_capacity) {
        allocator::deallocate(m_blocks, old_capacity);
      }
      m_blocks = tmpdata;
    }
  }


  template <typename IntType>       reference operator[](IntType n)       { return m_blocks[n >> shift][n & mask]; }
  template <typename IntType> const_reference operator[](IntType n) const { return m_blocks[n >> shift][n & mask]; }

  // TODO: check bounds in debug
  template <typename IntType>       reference at(IntType n)       { return m_blocks[n >> shift][n & mask]; }
  template <typename IntType> const_reference at(IntType n) const { return m_blocks[n >> shift][n & mask]; }

        reference front()       { return (*this)[0]; }
  const_reference front() const { return (*this)[0]; }

        reference back()       { return (*this)[m_size-1];}
  const_reference back() const { return (*this)[m_size-1];}

  template <typename InputIterator>
  typename boost::enable_if_c< !(boost::is_integral<InputIterator>::value || boost::is_enum<InputIterator>::value), void>::type
  assign( InputIterator first, InputIterator second)
  {
    clear();
    push_back(first, second);
  }

  void assign( size_type n, const_reference val)
  {
    clear();
    resize(n, val);
  }

  void push_back( const_reference val)
  {
    allocate(m_size+1);
    size_type i = m_size - 1;
    pointer buf = m_blocks[i >> shift] + (i & mask);
    block_allocator::construct(buf, val);
  }

  template <typename InputIterator>
  void push_back( InputIterator first, InputIterator second)
  {
    const size_type n = std::distance(first, second);
    const size_type curr_size = m_size;
    allocate(m_size+n);
    for (size_t i=curr_size; i<m_size; ++i, ++first) {
      pointer buf = m_blocks[ i >> shift ] + (i & mask);
      block_allocator::construct(buf, *first);
    }
  }

  void pop_back()
  {
    if (m_size > 0u) {
      --m_size;
      pointer p = m_blocks[m_size >> shift] + (m_size & mask);
      block_allocator::destroy(p);
    }
  }

  void pop_back(size_type n)
  {
    if (m_size > n ) {
      resize(m_size-n);
    }
    else {
      clear();
    }
  }

  iterator insert( iterator position, const_reference value)
  {
    size_type pos = position.m_offset;
    resize(size() + 1);
    // shift up
    for (size_type i = size()-1; i > pos; --i) {
      operator[](i) = operator[](i-1);
    }
    operator[](pos) = value;
    return position;
  }

  void insert(iterator position, size_type n, const_reference value)
  {
    size_type pos = position.m_offset;
    resize( size() + n );
    // shift up
    for (size_type i = size()-1; i > pos+n-1u; --i) {
      operator[](i) = operator[](i-n);
    }
    for (size_type i=0; i<n; ++i) {
      operator[](pos+i) = value;
    }
  }

  template <typename InputIterator>
  typename boost::enable_if_c< !(boost::is_integral<InputIterator>::value || boost::is_enum<InputIterator>::value), void>::type
  insert(iterator position, InputIterator first, InputIterator second)
  {
    const size_type n = std::distance(first, second);
    size_type pos = position.m_offset;
    resize( size() + n );
    // shift up
    for (size_type i = size()-1; i > pos+n-1u; --i) {
      operator[](i) = operator[](i-n);
    }
    for (size_type i=0; i < n; ++i, ++first) {
      operator[](pos+i) = *first;
    }
  }

  iterator erase(iterator position)
  {
    // shift down
    for (size_type pos = position.m_offset, pend = size()-1; pos < pend; ++pos) {
      operator[](pos) = operator[](pos+1);
    }
    pop_back();
    return position;
  }

  iterator erase(iterator first, iterator second)
  {
    const size_type n = std::distance(first, second);
    // shift down
    for (size_type pos = first.m_offset, pend = size()-n; pos < pend; ++pos) {
      operator[](pos) = operator[](pos+n);
    }
    pop_back(n);
    return first;
  }

  void swap(block_vector & rhs)
  {
    std::swap(m_size, rhs.m_size);
    std::swap(m_blocks_size, rhs.m_blocks_size);
    std::swap(m_blocks_capacity, rhs.m_blocks_capacity);
    std::swap(m_blocks, rhs.m_blocks);
  }

  void clear()
  {
    for (size_type i=0; i<m_size; ++i) {
      pointer p = m_blocks[i >> shift] + (i & mask);
      block_allocator::destroy(p);
    }
    m_size = 0;
  }

  template <typename TType, size_t BBlockSize, typename AAllocator>
  friend void std::swap( block_vector<TType, BBlockSize, AAllocator> & lhs
                        ,block_vector<TType, BBlockSize, AAllocator> & rhs );

private:

  void allocate(size_type num_entities)
  {
    const size_type required_blocks = (num_entities + block_size - 1) / block_size;


    // allocate new blocks
    if ( required_blocks > m_blocks_size ) {

      // grow block array
      if ( required_blocks > m_blocks_capacity ) {

        size_type old_capacity = m_blocks_capacity;

        m_blocks_capacity = m_blocks_capacity == 0 ? 64 : m_blocks_capacity << 1;
        m_blocks_capacity = m_blocks_capacity < required_blocks ? required_blocks : m_blocks_capacity;

        //copy over remaining blocks blocks
        value_type **tmpdata = allocator::allocate(m_blocks_capacity);

        memmove ( reinterpret_cast<void*>(tmpdata), reinterpret_cast<void*>(m_blocks), m_blocks_size * sizeof(value_type**));
        if (old_capacity) {
          allocator::deallocate(m_blocks, old_capacity);
        }
        m_blocks = tmpdata;
      }

      // allocate new blocks
      for (size_type block = m_blocks_size; block < required_blocks; ++block) {
        m_blocks[block] = block_allocator::allocate(block_size);
      }
      m_blocks_size = required_blocks;
    }
    m_size = num_entities;
  }

  //members
  size_type     m_size;
  size_type     m_blocks_size;
  size_type     m_blocks_capacity;
  value_type ** m_blocks;
};


template <typename Type, size_t BlockSize , typename Allocator>
bool operator==( const block_vector<Type,BlockSize,Allocator> & lhs
                ,const block_vector<Type,BlockSize,Allocator> & rhs )
{
  return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

template <typename Type, size_t BlockSize , typename Allocator>
bool operator!=( const block_vector<Type,BlockSize,Allocator> & lhs
                ,const block_vector<Type,BlockSize,Allocator> & rhs)
{
  return lhs.size() != rhs.size() || !std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

template <typename Type, size_t BlockSize , typename Allocator>
bool operator<( const block_vector<Type,BlockSize,Allocator> & lhs
               ,const block_vector<Type,BlockSize,Allocator> & rhs)
{
  return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

template <typename Type, size_t BlockSize , typename Allocator>
bool operator<=( const block_vector<Type,BlockSize,Allocator> & lhs
                ,const block_vector<Type,BlockSize,Allocator> & rhs)
{
  return (lhs == rhs) || (lhs < rhs);
}

template <typename Type, size_t BlockSize , typename Allocator>
bool operator>( const block_vector<Type,BlockSize,Allocator> & lhs
               ,const block_vector<Type,BlockSize,Allocator> & rhs)
{
  return rhs < lhs;
}

template <typename Type, size_t BlockSize , typename Allocator>
bool operator>=( const block_vector<Type,BlockSize,Allocator> & lhs
                ,const block_vector<Type,BlockSize,Allocator> & rhs)
{
  return (lhs == rhs) || (rhs < lhs);
}




} // namespace stk

namespace std {

template <typename Type, size_t BlockSize, typename Allocator>
inline void swap(   stk::block_vector<Type,BlockSize,Allocator> & lhs
                  , stk::block_vector<Type,BlockSize,Allocator> & rhs )
{ lhs.swap(rhs); }

} //namespace std


#endif //STK_UTIL_BLOCK_VECTOR_HPP
