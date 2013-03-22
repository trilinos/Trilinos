#ifndef STK_UTIL_DIAG_Mapv_h
#define STK_UTIL_DIAG_Mapv_h

#include <stk_util/util/FeatureTest.hpp>
#include <cstddef>

// ---------------------------------------------------------------------
// Author: H. Carter Edwards
//
// Purpose:
//   A map-style associative container for large allocated data objects.
// ---------------------------------------------------------------------
// Usage / Instantiation:
//
// Each node in a "virtual map" must be derived from the 'MapvNode'
// base class.  Their should be only one 'MapvNode' base class.
//
//  class MY_CLASS : public MapvNode<Key_Type> , other_base_classes {
//    // 'MY_CLASS' internal data & methods
//  };
//
//  class MY_CLASS2 : MY_CLASS {
//
//  };
//
// The Key_Type template argument is the data type of the ordering key
// for the virtual map.  The Key_Type must support the 'less<Key_Type>'
// template function object.
// ---------------------------------------------------------------------
// Features and Limitations:
//
// The 'MapvNode' derived objects may be inserted into a 'Mapv' container
// without being copied.
//
// When a 'MapvNode' derived object is destroyed it is automatically
// removed from the 'Mapv' container.  Destruction of a 'Mapv' container
// automatically invokes the 'delete' operator for each 'MapvNode'
// derived object that resides in the container at the time of its
// destruction.
//
// The 'insert' and 'remove' operations cause the 'Mapv' container to
// rebalance its binary tree.  These rebalancing algorithms are lengthy.
// As such the majority of the 'insert' and 'remove' operations are
// implemented in the seperately compiled 'Mapv.C' file.  The subprograms
// in 'Mapv.C' are shared by all instantiations of 'Mapv', thus
// "code bloat" is kept to a minimum.
//
// Alteration of a 'mavnode' derived objects 'key' value is likely
// to invalidate its position in the binary.  Thus alteration of the
// 'key' value removes the object from its 'Mapv' container.
//
// ---------------------------------------------------------------------
// Acknowledgements:
//
//   Most all of the algorithms in this class were obtained from
// the Hewlett-Packard source for the Standard Template Library,
// thus the inclusion of Hewlett-Packard's copyright notice.
// ---------------------------------------------------------------------
/*
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Hewlett-Packard Company makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 */
/*
Red-black tree class, designed for use in implementing STL
associative containers (set, multiset, map, and multimap).
The insertion and deletion algorithms are based on those in Cormen,
Leiserson, and Rivest, Introduction to Algorithms (MIT Press, 1990),
except that

(1) the header cell is maintained with links not only to the root
but also to the leftmost node of the tree, to enable constant time
begin(), and to the rightmost node of the tree, to enable linear time
performance when used with the generic set algorithms (set_union,
etc.);

(2) when a node being deleted has two children its successor node is
relinked into its place, rather than copied, so that the only
iterators invalidated are those referring to the deleted node.
*/
// ---------------------------------------------------------------------
// Don't force usage of 'assert', but if available then use it

#ifdef MAPV_ASSERT_H
#define mapv_assert( expr )	assert( expr )
#else
#define mapv_assert( expr )	/* do nothing */
#endif

// STL includes:

#include <utility>
#include <iterator>
#include <functional>

#ifdef SIERRA_IA64_OPTIMIZER_FIX
#pragma optimize("", off)
#endif

namespace sierra {

// class MyType : public MapvNode<MyKey> { ... };

template < typename Key_Type , class Key_Compare > class MapvNode ;

template < class Derived_Type , class Memory_Policy > class Mapv ;

template < class Derived_Type , class Next , class Prev > class MapvIterator ;

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
/** A very private base class for Mapv nodes */

class MapvIterNext ;
class MapvIterPrev ;

class MapvNodeBase {
private:

  enum { black = 0 , red = 1 };

  MapvNodeBase * parent ; // Binary tree node
  MapvNodeBase * left ;   // Binary tree node
  MapvNodeBase * right ;  // Binary tree node
  unsigned       color ;  // Red-black color

  inline void remove_from_container();

  /** Destructor removes the node from its container */
  virtual ~MapvNodeBase() { remove_from_container(); }

  MapvNodeBase() : parent(0), left(0), right(0), color(black) {}

  // Friends of this very private class

  friend class MapvBase ;
  friend class MapvIterNext ;
  friend class MapvIterPrev ;
  template < class dType , class N, class P > friend class MapvIterator ;
  template < class dType , class mPolicy > friend class Mapv ;
  template < class kType , class kCompare > friend class MapvNode ;
};

// ---------------------------------------------------------------------
/** Base class for Mapv nodes.
 *  Objects stored in a Mapv container are derived from this
 *  template base class as follows.
 *
 *    class Derived_Type : [ other classes , ]
 *                         public MapvNode<Key_Type>
 *                         [ , other classes ] { ... };
 *
 * OR
 *    class Derived_Type : [ other classes , ]
 *                         public MapvNode<Key_Type,Key_Less>
 *                         [ , other classes ] { ... };
 */

template < typename Key_Type , class Key_Compare = std::less<Key_Type> >
class MapvNode : private MapvNodeBase {
public:

  typedef Key_Type    key_type ;
  typedef Key_Compare key_compare ;

  bool mapv_valid() const { return this && MapvNodeBase::parent ; }
  void mapv_remove()      { MapvNodeBase::remove_from_container(); }

  /** Query key */
  const key_type & mapv_key() const { return Key ; }

  /** Set key, => forces removal of the node from its 'Mapv' container */
  const Key_Type & mapv_key( const key_type & K )
    { remove_from_container(); return Key = K  ; }

  /** A destroyed node is automatically removed from its container */
  virtual ~MapvNode() {}
  MapvNode() {}
  MapvNode( const key_type & K ) : Key(K) {}

private:

  Key_Type Key ;

  // Disallow copy constructor and assignment operator

  MapvNode( const MapvNode<Key_Type,Key_Compare> & );

  MapvNode<Key_Type,Key_Compare> &
    operator = ( const MapvNode<Key_Type,Key_Compare> & );

  // friends to access the MapvNodeBase base class and Key member
  template< class dType , class N , class P > friend class MapvIterator ;
  template< class dType , class mPolicy > friend class Mapv ;
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
/** Iterators for Mapv containers are bidirectional */

class MapvIterNext {
public:
  typedef       MapvNodeBase *  ptr ;
  typedef const MapvNodeBase * cptr ;
  static ptr op( cptr x )
    {
      ptr y = x->right;
      if ( y ) {
	while ( y->left ) y = y->left ;
      }
      else {
	mapv_assert( x->parent /* incrementing 'end()' iterator */ );
	y = x->parent ;
	while ( x == y->right ) y = ( x = y )->parent ;
	if ( x == y->parent ) y = y->right ;
      }
      return y ;
    }
};

class MapvIterPrev {
public:
  typedef       MapvNodeBase *  ptr ;
  typedef const MapvNodeBase * cptr ;
  static ptr op( cptr x )
    {
      ptr y = x->left;
      if ( y ) {
	while ( y->right ) y = y->right ;
      }
      else {
	mapv_assert( x->parent /* decrementing 'end()' iterator */ );
	y = x->parent ;
	while ( x == y->left ) y = ( x = y )->parent ;
	if ( x == y->parent ) y = y->left ;
      }
      return y ;
    }
};

template < class Type , class Next , class Prev >
class MapvIterator : public std::iterator<std::bidirectional_iterator_tag,Type>
{
private:
  template < class T , class M > friend class Mapv ;
  template < class T , class N, class P > friend class MapvIterator ;

  Type * n ;

  /** Dereference of 'end()' or 'rend()' iterators
   *  is erroneous and can lead to corruption of the container
   */
  Type * ptr() const { return n->MapvNodeBase::parent == 0 ? (Type*) 0 : n ; }

  MapvIterator( MapvNodeBase * x ) : n( static_cast<Type*>(x) )  {}
  MapvIterator( Type * x ) : n( x )  {}
  MapvIterator( Type & x ) : n( &x ) {}

public:

  typedef MapvIterator<Type,Next,Prev> SelfType ;

  /** Construct iterator from other similar iterators.
   *  The contained pointer types must be compile-time compatible
   *  for correctness of the pointer copy-constructor.
   */
  template<class T,class N, class P>
    MapvIterator( const MapvIterator<T,N,P> & x ) : n( x.n ) {}

  /** Assign iterator from other similar iterators.
   *  The contained pointer types must be compile-time compatible
   *  for correctness of the pointer assignment operator.
   */
  template<class T,class N,class P>
    SelfType & operator = ( const MapvIterator<T,N,P> & x )
      { n = x.n ; return *this ; }

  /** Query if the iterator is valid to dereference
   *  1) Is not NULL
   *  2) Is not an 'end()' or 'rend()' iterator.
   */
  bool valid_to_dereference() const { return n && n->MapvNodeBase::parent ; }

//   /** Conversion to bool: true for dereferenceable iterator */
//   operator bool () const { return valid_to_dereference(); }

  //--------------------------------------------------------------------

  MapvIterator() : n(0) {}
  MapvIterator( const SelfType & x ) : n( x.n ) {}

  SelfType & operator = ( const SelfType & x ) { n = x.n ; return *this ; }

  Type & operator * () const
    { mapv_assert( valid_to_dereference() ); return *n ; }
  Type * operator ->() const
    { mapv_assert( valid_to_dereference() ); return n ; }

  SelfType & operator++(){ n = static_cast<Type*>(Next::op(n)); return *this; }
  SelfType & operator--(){ n = static_cast<Type*>(Prev::op(n)); return *this; }

  SelfType operator++(int)
    { Type * t = n ; n = static_cast<Type*>(Next::op(n)); return SelfType(t); }

  SelfType operator--(int)
    { Type * t = n ; n = static_cast<Type*>(Prev::op(n)); return SelfType(t); }

  template<class T,class N, class P>
    bool operator == ( const MapvIterator<T,N,P> & y ) const
      { return n == y.n ; }

  template<class T,class N, class P>
    bool operator != ( const MapvIterator<T,N,P> & y ) const
      { return n != y.n ; }
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
/** Very private base class for Mapv container */

class MapvBase : private MapvNodeBase {
private:

  MapvNodeBase left_end ;  // 'rend()' node
  MapvNodeBase right_end ; // 'end()' node

  size_t Count ;

  //------------------------------------

  typedef MapvNodeBase   nType ;
  typedef MapvNodeBase * pType ;

  static MapvBase * container( const nType * x )
    {
      MapvBase * h = 0 ;
      if ( x && x->parent ) {
	// Search for root node, while loop breaks at root node
	while ( x != x->parent->parent ) x = x->parent ;
	// the root node's parent is the header
	h = static_cast<MapvBase*>( x->parent );
      }
      return h ;
    }

// ---------------------------------------------------------------------

  static nType * minimum( nType * );
  static nType * maximum( nType * );

  void rotate_left(  nType * x );
  void rotate_right( nType * x );

  nType * header()    const { return const_cast<nType*>((const nType*)this); }
  nType * rightmost() const { return right_end.left ; }
  nType * leftmost()  const { return left_end.right ; }

  void rightmost( nType * N ) { right_end.left = N ; }
  void leftmost(  nType * N ) { left_end.right = N ; }
  void root(      nType * N ) { header()->parent = N ; }

  //------------------------------------

  virtual ~MapvBase();

  void remove( MapvNodeBase * );

  nType * nRoot() const { return header()->parent ; }
  nType * nEnd()  const { return const_cast<nType*>( & right_end ); }
  nType * nREnd() const { return const_cast<nType*>( & left_end ); }

  nType * nBegin() const
    { return ( left_end.right != nREnd() ) ? left_end.right : nEnd(); }

  nType * nRBegin() const
    { return ( right_end.left != nEnd() ) ? right_end.left : nREnd(); }

  MapvBase() : MapvNodeBase(), left_end(), right_end(), Count(0)
    {
      header()->color = red ;  /* Color the header node red */
      leftmost(  header()->left  = nREnd() ); // left end of the tree
      rightmost( header()->right = nEnd() );  // right end of the tree
    }

  void insert( nType * y , nType * z , bool z_lt_y );

  nType * unbalancing_removal( nType ** n );
  static void WarnOptimize();

  friend class MapvNodeBase ;
  template< class dType , class mPolicy > friend class Mapv ;
};

inline void MapvNodeBase::remove_from_container()
{
  MapvBase * const c = MapvBase::container(this);
  if ( c ) c->remove( this );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template <class Derived_Type>
struct MapvNewDeletePolicy {
  typedef typename Derived_Type::key_type key_type ;

  Derived_Type * create( const key_type & K ) { return new Derived_Type(K); }

  void destroy( Derived_Type * p ) { delete p ; }
};

template <class Derived_Type>
struct MapvDeleteOnlyPolicy {
  typedef typename Derived_Type::key_type key_type ;

  Derived_Type * create( const key_type & K ) { return (Derived_Type*) 0 ; }

  void destroy( Derived_Type * p ) { delete p ; }
};

template <class Derived_Type>
struct MapvNullPolicy {
  typedef typename Derived_Type::key_type key_type ;

  Derived_Type * create( const key_type & K ) { return (Derived_Type*) 0 ; }

  void destroy( Derived_Type * p ) {}
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template < class Derived_Type ,
	   class Memory_Policy = MapvDeleteOnlyPolicy<Derived_Type> >
class Mapv : public MapvBase {
public:

  // 'self' type
  typedef Mapv<Derived_Type,Memory_Policy> SelfType ;

  // NOTE: value_type == Derived_Type, the chain of
  //       typedefs insures the correct derivation has occured.

  typedef typename Derived_Type::key_type     key_type ;
  typedef typename Derived_Type::key_compare  key_compare ;
  typedef Derived_Type                        value_type ;
  typedef size_t                              size_type ;
  typedef ptrdiff_t                           difference_type ;

  typedef       value_type * pointer ;
  typedef const value_type * const_pointer ;
  typedef       value_type & reference ;
  typedef const value_type & const_reference ;

  struct value_compare
    : public std::binary_function<value_type,value_type,bool>
  {
  protected:
    key_compare comp ;
  public:
    bool operator()(const value_type& x, const value_type& y) const
      { return comp( x.mapv_key() , y.mapv_key() ); }
  };

  typedef MapvIterNext Next ;
  typedef MapvIterPrev Prev ;

  typedef MapvIterator<      Derived_Type,Next,Prev> iterator ;
  typedef MapvIterator<const Derived_Type,Next,Prev> const_iterator ;
  typedef MapvIterator<      Derived_Type,Prev,Next> reverse_iterator ;
  typedef MapvIterator<const Derived_Type,Prev,Next> const_reverse_iterator ;

  typedef std::pair<iterator, bool> pair_iterator_bool ;

private:

  // Disallow copy and assignment

  Mapv( const SelfType & );

  SelfType & operator = ( const SelfType & );

  Memory_Policy memoryPolicy ;
  key_compare   key_less ;      // An "abstract" object
  value_compare value_less ;    // An "abstract" object

  static const key_type & key( nType * const n )
    { return static_cast<MapvNode<key_type,key_compare>*>(n)->mapv_key(); }

  static value_type * cast( nType * n )
    { return static_cast<value_type*>(n); }

// WORKAROUND: 4/7/2010 -- This works for our SGI Itanium Intel compiler on sasg1132 
#ifdef SIERRA_IA64_OPTIMIZER_FIX
   Derived_Type * lb( const key_type & k ) const
     {
      volatile pType y = nEnd();  // [DGB] -- Altix intel-10.1 optimizer hack
      volatile pType x = nRoot(); // [DGB] -- Altix intel-10.1 optimizer hack
      while ( x ) {
        bool go_right = key_less(key(x), k);
        
        if (go_right) {
          x = x->right;
        }
        else {
          y = x;
          x = x->left;
        }
      }      
      return cast(y);
    }

  Derived_Type * ub( const key_type & k ) const
    {
      pType y = nEnd();
      pType x = nRoot();
      while ( x ) x = key_less(k,key(x)) ? ( y = x )->left : x->right ;
      return cast(y);
    }

#else
  Derived_Type * lb( const key_type & k ) const
    {
      pType y = nEnd();
      pType x = nRoot();
      while ( x ) x = key_less(key(x),k) ? x->right : ( y = x )->left ;
      return cast(y);
    }

  Derived_Type * ub( const key_type & k ) const
    {
      pType y = nEnd();
      pType x = nRoot();
      while ( x ) x = key_less(k,key(x)) ? ( y = x )->left : x->right ;
      return cast(y);
    }
#endif // SIERRA_IA64_OPTIMIZER_FIX

  Derived_Type * f( const key_type & k ) const
    {
      nType * const e = nEnd();
      nType * const y = lb(k); // k <= y->mapv_key()

      // If 'end()' or k < y->mapv_key() then not found
      return cast( ( y == e || key_less(k,key(y)) ) ? e : y );
    }

public:

  static SelfType * container( const Derived_Type & n )
    {
      MapvBase * const c = MapvBase::container(&n);
      return c ? static_cast<SelfType*>( c ) : (SelfType*) 0 ;
    }

  static SelfType * container( const Derived_Type * n )
    { return n ? container( *n ) : (SelfType*) 0 ; }

  // -------------------------------------------------------------------

  key_compare   key_comp()   const { return key_less ; }
  value_compare value_comp() const { return value_less ; }

  const_iterator begin()  const { return const_iterator(cast(nBegin())); }
  const_iterator end()    const { return const_iterator(cast(nEnd())); }
  const_iterator rbegin() const { return const_iterator(cast(nRBegin())); }
  const_iterator rend()   const { return const_iterator(cast(nREnd())); }

  iterator begin()  { return iterator(cast(nBegin())); }
  iterator end()    { return iterator(cast(nEnd())); }
  iterator rbegin() { return iterator(cast(nRBegin())); }
  iterator rend()   { return iterator(cast(nREnd())); }

  bool      empty() const { return MapvBase::Count == 0 ; }
  size_type size()  const { return MapvBase::Count ; }

  // -------------------------------------------------------------------
  // search operations:

  iterator lower_bound( const key_type & k ) { return iterator( lb(k) ); }
  iterator upper_bound( const key_type & k ) { return iterator( ub(k) ); }

  const_iterator lower_bound( const key_type & k ) const
    { return const_iterator( lb(k) ); }

  const_iterator upper_bound( const key_type & k ) const
    { return const_iterator( ub(k) ); }

  iterator find( const key_type & k ) { return iterator( f(k) ); }

  const_iterator find( const key_type & k ) const
    { return const_iterator( f(k) ); }

  // -------------------------------------------------------------------

  /** Return member with the input key value.
   *  If does not already existing then create and insert.
   */
  value_type & operator[]( const key_type & k )
    {
      pType y = nEnd();
      pType x = nRoot();

      bool flag = true ;

      while ( x )
	{ y = x ; x = ( flag = key_less(k, key(x)) ) ? x->left : x->right ; }

      /* flag = k < y , check previous value if exists */

      const bool k_lt_y = flag ;

      x = flag && y != nBegin() ? ( flag = false , MapvIterPrev::op(y) ) : y ;

      if ( flag || ( flag = key_less( key(x) , k ) ) ) {
	x = memoryPolicy.create(k);
	MapvBase::insert( y , x , k_lt_y );
      }
      return *static_cast<pointer>(x);
    }

  /** Insert an object into this container.
   *  The input object is itself inserted, no copy is made.
   *  If the input object is already a member of a container
   *  then it is removed from that container.
   *  Return iterator to member with the key value,
   *  if already existing the (iterator,false) is returned.
   */

  pair_iterator_bool insert( const pointer v )
    {
      WarnOptimize();
      pType y = nEnd();
      pType x = nRoot();

      bool flag = true ;

      while ( x ) {
	y = x ;
	x = ( flag = key_less(v->mapv_key(), key(x)) ) ? x->left : x->right ;
      }

      /* flag = k < y , check previous value if exists */

      const bool k_lt_y = flag ;

      x = flag && y != nBegin() ? ( flag = false , MapvIterPrev::op(y) ) : y ;

      if ( flag || ( flag = key_less( key(x) , v->mapv_key() ) ) ) {
	x = v ;
	MapvBase::insert( y , x , k_lt_y );
      }
      return pair_iterator_bool( iterator(x), flag );
    }

  pair_iterator_bool insert( reference v ) { return insert( &v ); }

  // -------------------------------------------------------------------
  // Remove & erase operations

  pointer remove( const key_type & k )
    {
      pointer v = f(k);
      return ( nEnd() != v ) ? ( MapvBase::remove(v) , v ) : (pointer) 0 ;
    }

  pointer remove( iterator i )
    { MapvBase::remove( i.n ); return i.n ; }

  void erase( const key_type & K )
    { pointer v = remove(K); if ( v ) memoryPolicy.destroy(v); }

  void erase( iterator i )
    { pointer v = remove(i); if ( v ) memoryPolicy.destroy(v); }

  // -------------------------------------------------------------------
  // construction / destruction
  // Destruction will apply the 'Derived_Type::destroy'
  // method to each node in the tree.

  Mapv() {}

  void clear()
    {
      if ( Count ) {
	nType * n = nBegin();
	nType * t ;
	while ( ( t = unbalancing_removal( &n ) ) ) {
	  memoryPolicy.destroy( cast( t ) );
	}
      }
    }

  virtual ~Mapv() { clear(); }

  // -------------------------------------------------------------------
  // Verify the integrity of the tree

  bool verify() const
    {
      size_type count = 0 ;
      const_iterator  i = begin();
      const_iterator  j ;

      while ( i != end() &&
	     ( ++(j=i) == end() || key_less(i->mapv_key(), j->mapv_key()) )
	     && --j == i )
	{ ++i ; ++count ; }

      return ( i == end() && count == size() );
    }
};

// ---------------------------------------------------------------------

// For the 'no_delete' operator.

template < class Derived_Type >
class Mapv_no_delete
  : public Mapv<Derived_Type,MapvNullPolicy<Derived_Type> > {};

// ---------------------------------------------------------------------

}
#ifdef SIERRA_IA64_OPTIMIZER_FIX
#pragma optimize("", on)
#endif

#undef mapv_assert

#endif // STK_UTIL_DIAG_Mapv_h
