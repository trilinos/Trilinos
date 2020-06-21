// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_NoMallocArray_hpp
#define percept_NoMallocArray_hpp

#include <algorithm>
#include <percept/Util.hpp>

  namespace percept {

    // patterned after std::array

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

    public:
      size_type m_size;  // 0...N-1
      T m_data[N];    // fixed-size array of elements of type T

    public:

      explicit NoMallocArray() : m_size(0u) {}

      NoMallocArray(size_type sz, const T& val) :  m_size(sz)
      {
#ifdef NOMALLOC_ARRAY_CHECK_SIZES
#endif
#if 1
        if (sz > N) throw std::runtime_error("error n > N");
        for (size_type i = 0; i < sz; i++)
          {
            m_data[i] = val;
          }
#endif
      }

#if 0
      NoMallocArray(const NoMallocArray& nma) : m_size(nma.m_size) {
        for (size_type i = 0; i < m_size; i++)
          {
            m_data[i] = nma.m_data[i];
          }
      }
#endif

      void clear() { m_size=0; }

      // iterator support
      iterator begin() { return m_data; }
      const_iterator begin() const { return m_data; }
      iterator end() { return m_data+m_size; }
      const_iterator end() const { return m_data+m_size; }

      void resize(size_type n) {
        if (n > N) throw std::runtime_error("error n > N");
        m_size=n;
      }
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

      inline size_type size() const { return m_size;}

      void insert(T val)
      {
        VERIFY_OP(m_size, < ,  N, "out of bounds");
        (*this)[m_size] = val;
        m_size++;
        //if (m_size > m_capacity) throw std::runtime_error("m_size");
        //std::sort( begin(), end() );
      }

      bool empty() { return m_size == 0; }
      size_type max_size() { return N; }
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

      // assignment with type conversion
      template <typename T2>
      NoMallocArray<T,N>& operator= (const NoMallocArray<T2,N>& rhs) {
        std::copy(rhs.begin(),rhs.end(), begin());
        return *this;
      }

#if 0
      bool operator== (const NoMallocArray<T,N>& rhs) const {
        for (size_type i = 0; i < size(); i++)
          {
            if (m_data[i] != rhs.m_data[i])
              return false;
          }
        return true;
      }
#endif

      // assign one value to all elements
      void assign (const T& value)
      {
        std::fill_n(begin(),size(),value);
      }

      // check range (may be private because it is static)
       void rangecheck (size_type i) const {
        if (i >= size()) {
          throw std::out_of_range("NoMallocArray<>: index out of range");
        }
      }
    };

    template<class T, std::size_t N>
    inline std::ostream &operator<<(std::ostream& out, const NoMallocArray<T,N>& arr)
    {
      out << arr[0];
      return out;
    }
  }
#endif
