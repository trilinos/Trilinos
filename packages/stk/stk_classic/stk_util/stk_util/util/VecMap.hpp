/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
/**
 * @file
 * @author H. Carter Edwards
 * @date   October 2002
 */

#ifndef STK_UTIL_UTIL_vecmap_hpp
#define STK_UTIL_UTIL_vecmap_hpp

#include <utility>
#include <functional>
#include <vector>
#include <algorithm>

namespace sierra {

/**
 * @class vecmap
 * @brief Vector-based std::map functionality
 *
 * @par Purpose: Mimic the 'std::map' interface
 *
 *  This template class mimics the 'std::map'
 *  associative container interface; however,
 *  its semantics are significantly different.
 *  Storage for the map-lite class is provided
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
 *  Modifications to the vecmap contents are linear
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

template <class Key, class T, class Compare = std::less<Key> >
class vecmap {
public:
  typedef Key      key_type ;
  typedef T        mapped_type ;
  typedef Compare  key_compare ;

private: // Hidden storage type
  typedef std::vector< std::pair<const key_type,mapped_type> > storage ;

public:
  typedef typename storage::value_type             value_type ;
  typedef typename storage::allocator_type         allocator_type ;
  typedef typename allocator_type::reference       reference ;
  typedef typename allocator_type::const_reference const_reference ;
  typedef typename allocator_type::pointer         pointer ;
  typedef typename allocator_type::const_pointer   const_pointer ;
  typedef typename storage::size_type              size_type ;
  typedef typename storage::difference_type        difference_type ;
  typedef typename storage::iterator               iterator ;
  typedef typename storage::const_iterator         const_iterator ;
  typedef typename storage::reverse_iterator       reverse_iterator ;
  typedef typename storage::const_reverse_iterator const_reverse_iterator ;

  typedef std::pair<key_type,mapped_type> value_type_unconst_key ;

private: // key compare functors
  class value_compare
    : public std::binary_function<value_type,value_type,bool> {
  private:
    key_compare comp ;
  public:
    bool operator ()( const value_type & RHS , const value_type & LHS ) const
    { return comp( RHS.first , LHS.first ); }
  };

  class value_compare_key
    : public std::binary_function<value_type,key_type,bool> {
  private:
    key_compare comp ;
  public:
    bool operator()( const value_type & RHS , const key_type & LHS ) const
    { return comp( RHS.first , LHS ); }
  };

  class value_compare_unconst_key
    : public std::binary_function<value_type_unconst_key,key_type,bool> {
  private:
    key_compare comp ;
  public:
    bool operator()( const value_type_unconst_key & RHS , const key_type & LHS ) const
    { return comp( RHS.first , LHS ); }
  };

private: // Object data
  std::vector<value_type_unconst_key>	Storage;
  value_compare			ValueComp ;
  value_compare_key			ValueKeyComp ;
  value_compare_unconst_key		UnconstValueKeyComp ;
  key_compare				KeyComp ;

private:
  storage & pub() {
    char &i = reinterpret_cast<char&>(Storage);
    return reinterpret_cast<storage&>(i);
  }

  const storage & pub() const {
    const char &i = reinterpret_cast<const char&>(Storage);
    return reinterpret_cast<const storage&>(i);
  }

  typedef typename std::vector<value_type_unconst_key>::iterator iterator_private ;

  static iterator_private & castit( iterator & i )
  { return reinterpret_cast<iterator_private &>(i); }

  static iterator & castit( iterator_private & i )
  { return reinterpret_cast<iterator &>(i); }

public:
  ~vecmap() {}

  vecmap() : Storage() {}

  vecmap( const vecmap<Key,T,Compare> & rhs ) : Storage(rhs.Storage) {}

  vecmap<Key,T,Compare> & operator = ( const vecmap<Key,T,Compare> & rhs ) {
    Storage = rhs.Storage ;
    return *this ;
  }

  void swap( vecmap<Key,T,Compare> & v ) {
    Storage.swap( v.Storage );
  }

  iterator               begin()          { return pub().begin(); }
  iterator               end()            { return pub().end(); }
  const_iterator         begin() const    { return pub().begin(); }
  const_iterator         end() const      { return pub().end(); }
  reverse_iterator       rbegin()         { return pub().rbegin(); }
  reverse_iterator       rend()           { return pub().rend(); }
  const_reverse_iterator rbegin() const   { return pub().rbegin(); }
  const_reverse_iterator rend() const     { return pub().rend(); }
  bool                   empty() const    { return Storage.empty(); }
  size_type              size() const     { return Storage.size(); }
  size_type              max_size() const { return Storage.max_size(); }


  typename std::vector<value_type_unconst_key>::iterator lower_bound_private_comp_( const key_type & k ) {
    typename std::vector<value_type_unconst_key>::iterator __first = Storage.begin();
    typename std::vector<value_type_unconst_key>::iterator __last = Storage.end();

//    difference_type __len = std::distance(__first, __last);
    difference_type __len = __last - __first;
    difference_type __half;
    typename std::vector<value_type_unconst_key>::iterator __middle;

    while (__len > 0) {
      __half = __len >> 1;
      __middle = __first;
//      std::advance(__middle, __half);
      __middle += __half;
      if (UnconstValueKeyComp(*__middle, k)) {
	__first = __middle;
	++__first;
	__len = __len - __half - 1;
      }
      else
	__len = __half;
    }
    return __first;
  }

  std::pair<iterator,bool> insert( const value_type & v ) {
    typename std::vector<value_type_unconst_key>::iterator ip =
      lower_bound_private_comp_(v.first);
    // (ip-1)->first < v.first <= ip->first
    const bool b = Storage.end() == ip || KeyComp( v.first , ip->first );
    // b = v.first != ip->first
    if ( b ) ip = Storage.insert(ip,value_type_unconst_key(v.first,v.second));
    return std::pair<iterator,bool>( castit(ip), b );
  }

  mapped_type & operator[]( const key_type & k ) {
    typename std::vector<value_type_unconst_key>::iterator ip =
      lower_bound_private_comp_(k);

    // (ip-1)->first < k <= ip->first
    if ( Storage.end() == ip || KeyComp(k,ip->first) ) {
      // k != ip->first => insert
      ( ip = Storage.insert(ip,value_type_unconst_key()) )->first = k ;
    }
    return ip->second ;
  }

  void erase( iterator i ) {
    Storage.erase( castit(i) );
  }

  void erase( iterator first , iterator last ) {
    Storage.erase( castit(first), castit(last) );
  }

  size_type erase( const key_type & k ) {
    const typename std::vector<value_type_unconst_key>::iterator i =
      lower_bound_private_comp_(k);
    return KeyComp( k , i->first ) ? 0 : ( Storage.erase(i) , 1 );
  }

  void clear() {
    Storage.clear();
  }

  key_compare   key_comp() const   { return KeyComp ; }

  value_compare value_comp() const { return ValueComp ; }

  iterator lower_bound( const key_type & k ) {
    iterator __first = begin();
    iterator __last = end();

//    difference_type __len = std::distance(__first, __last);
    difference_type __len = __last - __first;
    difference_type __half;
    iterator __middle;

    while (__len > 0) {
      __half = __len >> 1;
      __middle = __first;
//      std::advance(__middle, __half);
      __middle += __half;
      if (ValueKeyComp(*__middle, k)) {
	__first = __middle;
	++__first;
	__len = __len - __half - 1;
      }
      else
	__len = __half;
    }
    return __first;
  }

  const_iterator lower_bound( const key_type & k ) const {
    const_iterator __first = begin();
    const_iterator __last = end();

//    difference_type __len = std::distance(__first, __last);
    difference_type __len = __last - __first;
    difference_type __half;
    const_iterator __middle;

    while (__len > 0) {
      __half = __len >> 1;
      __middle = __first;
//      std::advance(__middle, __half);
      __middle += __half;
      if (ValueKeyComp(*__middle, k)) {
	__first = __middle;
	++__first;
	__len = __len - __half - 1;
      }
      else
	__len = __half;
    }
    return __first;
  }

  iterator upper_bound( const key_type & k ) {
    iterator __first = begin();
    iterator __last = end();

//    difference_type __len = std::distance(__first, __last);
    difference_type __len = __last - __first;
    difference_type __half;
    iterator __middle;

    while (__len > 0) {
      __half = __len >> 1;
      __middle = __first;
//      std::advance(__middle, __half);
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


  const_iterator upper_bound( const key_type & k ) const {
    const_iterator __first = begin();
    const_iterator __last = end();

//    difference_type __len = std::distance(__first, __last);
    difference_type __len = __last - __first;
    difference_type __half;
    const_iterator __middle;

    while (__len > 0) {
      __half = __len >> 1;
      __middle = __first;
//      std::advance(__middle, __half);
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


  iterator find( const key_type & k ) {
    const iterator i = lower_bound(k);
    return end() == i || KeyComp( k , i->first ) ? end() : i ;
  }

  const_iterator find( const key_type & k ) const {
    const const_iterator i = lower_bound(k);
    return end() == i || KeyComp( k , i->first ) ? end() : i ;
  }

  size_type count( const key_type & k ) const {
    const const_iterator i = lower_bound(k);
    return end() == i || KeyComp( k , i->first ) ? 0 : 1 ;
  }

  //--------------------------------------------------------------------
  // Is ok to get constant versions of the underlyingstorage

  operator const std::vector<value_type> & () const {
    return * reinterpret_cast<const std::vector<value_type>*>( &Storage );
  }

  operator const std::vector< std::pair<key_type,mapped_type> > & () const {
    return Storage ;
  }

  void reserve( size_type n ) {
    Storage.reserve(n);
  }

  bool operator == ( const vecmap<Key,T,Compare> & rhs ) const {
    return Storage == rhs.Storage ;
  }

  bool operator != ( const vecmap<Key,T,Compare> & rhs ) const {
    return Storage != rhs.Storage ;
  }
};

} // namespace sierra

#endif // STK_UTIL_UTIL_vecmap_hpp
