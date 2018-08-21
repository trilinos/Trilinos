// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef STK_UTIL_UTIL_vecset_hpp
#define STK_UTIL_UTIL_vecset_hpp

#include <utility>
#include <vector>
#include <functional>

namespace sierra {

/**
 * @class vecset
 * @brief Vector-based std::set functionality
 *
 * @par Purpose: Mimic the 'std::set' interface
 *
 *  This template class mimics the 'std::set'
 *  associative container interface; however,
 *  its semantics are significantly different.
 *  Storage for the set-lite class is provided
 *  by the std::vector class where the entries
 *  are sorted by key value.
 *
 * @par Domain of Applicability
 *
 *  Light weight associative container functionality
 *  for small keys and values, e.g. key = integer
 *  and value = pointer.
 *
 * @par Associative container violations
 *
 *  Modifications to the vecset contents are linear
 *  complexity in violation of the associative
 *  container requirement for logarithmic complexity.
 *  Furthermore, modification operations are guaranteed
 *  to invalidate all iterators after the insert/erase
 *  point.  Insertion operations may also invalidate
 *  all iterators if the storage is reallocated.
 *
 * @par Associative container compliance
 *
 *  All non-modifying query operations conform to
 *  either the constant or logarithmic complexity.
 */
template <class Key, class Compare = std::less<Key> >
class vecset {
private:
  template <typename T>
  struct const_key_type_meta_func
  {
    typedef const T type;
  };

  template <typename T>
  struct const_key_type_meta_func<T*>
  {
    typedef const T* const type;
  };

public:

  typedef Key      key_type;
  typedef Key      value_type;
  typedef Compare  key_compare;
  typedef Compare  value_compare;
  typedef typename const_key_type_meta_func<Key>::type const_key_type;

private:

  typedef std::vector<key_type> storage ;

public:

  typedef typename storage::allocator_type         allocator_type ;
  typedef typename allocator_type::reference       reference ;
  typedef typename allocator_type::const_reference const_reference ;
  typedef typename allocator_type::pointer         pointer ;
  typedef typename allocator_type::const_pointer   const_pointer ;
  typedef typename storage::size_type              size_type ;
  typedef typename storage::difference_type        difference_type ;
  typedef typename storage::iterator		   iterator ;
  typedef typename storage::const_iterator         const_iterator ;
  typedef typename storage::reverse_iterator       reverse_iterator ;
  typedef typename storage::const_reverse_iterator const_reverse_iterator ;

  //--------------------------------------------------------------------
private:

  key_compare KeyComp ;
  storage     Storage ;

public:
  //--------------------------------------------------------------------

  ~vecset()
  {}

  vecset() : Storage()
  {}

  vecset( const vecset<Key,Compare> & rhs ) : Storage(rhs.Storage)
  {}

  vecset<Key,Compare> & operator = ( const vecset<Key,Compare> & rhs ) {
    Storage = rhs.Storage ;
    return *this ;
  }

  void swap( vecset<Key,Compare> & v ) {
    Storage.swap( v.Storage );
  }

  iterator               begin()          { return Storage.begin(); }
  iterator               end()            { return Storage.end(); }
  const_iterator         begin() const    { return Storage.begin(); }
  const_iterator         end() const      { return Storage.end(); }
  reverse_iterator       rbegin()         { return Storage.rbegin(); }
  reverse_iterator       rend()           { return Storage.rend(); }
  const_reverse_iterator rbegin() const   { return Storage.rbegin(); }
  const_reverse_iterator rend() const     { return Storage.rend(); }
  bool                   empty() const    { return Storage.empty(); }
  size_type              size() const     { return Storage.size(); }

  iterator lower_bound( const_key_type & k ) {
    iterator __first = begin();
    iterator __last = end();

    difference_type __len = __last - __first;
    difference_type __half;
    iterator __middle;

    while (__len > 0) {
      __half = __len >> 1;
      __middle = __first;
      __middle += __half;
      if (KeyComp(*__middle, k)) {
	__first = __middle;
	++__first;
	__len = __len - __half - 1;
      }
      else
	__len = __half;
    }
    return __first;
  }

  const_iterator lower_bound( const_key_type & k ) const {
    const_iterator __first = begin();
    const_iterator __last = end();

    difference_type __len = __last - __first;
    difference_type __half;
    const_iterator __middle;

    while (__len > 0) {
      __half = __len >> 1;
      __middle = __first;
      __middle += __half;
      if (KeyComp(*__middle, k)) {
	__first = __middle;
	++__first;
	__len = __len - __half - 1;
      }
      else
	__len = __half;
    }
    return __first;
  }

  iterator upper_bound( const_key_type & k ) {
    iterator __first = begin();
    iterator __last = end();

    difference_type __len = __last - __first;
    difference_type __half;
    iterator __middle;

    while (__len > 0) {
      __half = __len >> 1;
      __middle = __first;
      __middle += __half;
      if (__comp(k, *__middle))
	__len = __half;
      else {
	__first = __middle;
	++__first;
	__len = __len - __half - 1;
      }
    }
    return __first;
  }


  const_iterator upper_bound( const_key_type & k ) const {
    const_iterator __first = begin();
    const_iterator __last = end();

    difference_type __len = __last - __first;
    difference_type __half;
    const_iterator __middle;

    while (__len > 0) {
      __half = __len >> 1;
      __middle = __first;
      __middle += __half;
      if (__comp(k, *__middle))
	__len = __half;
      else {
	__first = __middle;
	++__first;
	__len = __len - __half - 1;
      }
    }
    return __first;
  }

  std::pair<iterator,bool> insert( const value_type & v ) {
    typename storage::iterator ip = lower_bound(v);
    // (ip-1)->first < v <= ip->first
    const bool b = Storage.end() == ip || KeyComp( v , *ip );
    // b = v != *ip
    if ( b ) ip = Storage.insert(ip,v);
    return std::pair<iterator,bool>( ip, b );
  }

  template< typename FORWARD_ITERATOR >
  void insert( FORWARD_ITERATOR b, FORWARD_ITERATOR e) {
    for (FORWARD_ITERATOR it = b ; it != e ; ++it) {
        insert(*it);
    }
  }

  void erase( iterator i ) {
    Storage.erase( Storage.begin() + ( i - begin() ) );
  }

  void erase( iterator first , iterator last ) {
    Storage.erase( Storage.begin() + ( first - begin() ) ,
		   Storage.begin() + ( last  - begin() ) );
  }

  size_type erase( const_key_type & k ) {
    typename storage::iterator i = lower_bound(k );
    return KeyComp( k , *i ) ? 0 : ( Storage.erase(i) , 1 );
  }

  void clear() {
    Storage.clear();
  }

  key_compare   key_comp() const   { return KeyComp ; }
  value_compare value_comp() const { return KeyComp ; }

  iterator find( const_key_type & k ) {
    const iterator i = lower_bound(k);
    return end() == i || KeyComp( k , *i ) ? end() : i ;
  }

  const_iterator find( const_key_type & k ) const {
    const const_iterator i = lower_bound(k);
    return end() == i || KeyComp( k , *i ) ? end() : i ;
  }

  size_type count( const_key_type & k ) const {
    const const_iterator i = lower_bound(k);
    return end() == i || KeyComp( k , *i ) ? 0 : 1 ;
  }

  //--------------------------------------------------------------------
  // Allow conversion to constant vector

  operator const storage & () const {
    return Storage ;
  }

  void reserve( size_type n ) {
    Storage.reserve( n );
  }

  bool operator == ( const vecset<Key,Compare> & rhs ) const {
    return Storage == rhs.Storage ;
  }

  bool operator != ( const vecset<Key,Compare> & rhs ) const {
    return Storage != rhs.Storage ;
  }
};

}

#endif // STK_UTIL_UTIL_vecset_hpp
