#ifndef stk_percept_NoMallocArray_hpp
#define stk_percept_NoMallocArray_hpp

#include <algorithm>
#include <stk_percept/Util.hpp>

namespace stk {
  namespace percept {

    // patterned after boos::array

    template<class T, std::size_t N>
    class NoMallocArray {

    public:
      //enum { MAX_CAPACITY = N };

      // type definitions
      typedef T              value_type;
      typedef T*             iterator;
      typedef const T*       const_iterator;
      typedef T&             reference;
      typedef const T&       const_reference;
      typedef std::size_t    size_type;
      typedef std::ptrdiff_t difference_type;

    private:
      size_type m_current_max_size;
      size_type m_size;
      T m_data[N];    // fixed-size array of elements of type T

    public:

      NoMallocArray(size_type n=N) : m_current_max_size(n), m_size(0) 
      {
        if (n > N) throw std::runtime_error("error n > N");
        for (size_type i = 0; i < N; i++)
          {
            m_data[i] = 0.0;
          }
      }

      void clear() { m_size=0; }

      // iterator support
      iterator begin() { return m_data; }
      const_iterator begin() const { return m_data; }
      iterator end() { return m_data+m_size; }
      const_iterator end() const { return m_data+m_size; }

      // reverse iterator support
      typedef std::reverse_iterator<iterator> reverse_iterator;
      typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

      reverse_iterator rbegin() { return reverse_iterator(end()); }
      const_reverse_iterator rbegin() const {
        return const_reverse_iterator(end());
      }
      reverse_iterator rend() { return reverse_iterator(begin()); }
      const_reverse_iterator rend() const {
        return const_reverse_iterator(begin());
      }

      // operator[]
      reference operator[](size_type i) 
      { 
        VERIFY_OP( i, <,  N,  "out of range" ); 
        return m_data[i];
      }
        
      const_reference operator[](size_type i) const 
      {     
        VERIFY_OP( i, <,  N,  "out of range" ); 
        return m_data[i]; 
      }

      // at() with range check
      reference at(size_type i) { rangecheck(i); return m_data[i]; }
      const_reference at(size_type i) const { rangecheck(i); return m_data[i]; }
    
      // front() and back()
      reference front() 
      { 
        return m_data[0]; 
      }
        
      const_reference front() const 
      {
        return m_data[0];
      }
        
      reference back() 
      { 
        return m_data[m_size-1]; 
      }
        
      const_reference back() const 
      { 
        return m_data[m_size-1]; 
      }

      size_type size() const { return m_size;}

      void insert(T val) 
      {
        VERIFY_OP(m_size, < ,  m_current_max_size, "out of bounds");
        (*this)[m_size] = val;
        m_size++;
        //if (m_size > m_current_max_size) throw std::runtime_error("m_size");
        //std::sort( begin(), end() );
      }

      bool empty() { return m_size == 0; }
      size_type max_size() { return m_current_max_size; }
      size_type max_capacity() { return N; }
      

      bool contains(T val)
      {
        for (size_type i = 0; i < m_size; i++)
          {
            if (val == (*this)[i])
              {
                return true;
              }
          }
        return false;
      }

      // swap (note: linear complexity)
      void swap (NoMallocArray<T,N>& y) {
        std::swap_ranges(begin(),end(),y.begin());
      }

      // direct access to data (read-only)
      const T* data() const { return m_data; }
      T* data() { return m_data; }

      // use array as C array (direct read/write access to data)
      T* c_array() { return m_data; }

      // assignment with type conversion
      template <typename T2>
      NoMallocArray<T,N>& operator= (const NoMallocArray<T2,N>& rhs) {
        std::copy(rhs.begin(),rhs.end(), begin());
        return *this;
      }

      // assign one value to all elements
      void assign (const T& value)
      {
        std::fill_n(begin(),size(),value);
      }

      // check range (may be private because it is static)
       void rangecheck (size_type i) {
        if (i >= size()) {
          throw std::out_of_range("NoMallocArray<>: index out of range");
        }
      }
    };

  }
}
#endif
