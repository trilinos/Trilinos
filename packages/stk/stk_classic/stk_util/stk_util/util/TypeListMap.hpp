/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


#ifndef STK_UTIL_UTIL_TypeListMap_h
#define STK_UTIL_UTIL_TypeListMap_h

/**
 * @file
 * @author H. Carter Edwards
 * @date   August 2005
 */

#include <stk_util/util/TypeList.hpp>

/**
 * @file
 *
 */

namespace sierra {

/** Map of 'class Tag' to 'Tag::type' value.
 *  The 'class Tag' should be a 'tag'; i.e., it
 *  a) not be virtual,
 *  b) not contain any member data,
 *  c) have a no-op default constructor,
 *  d) have a no-op destructor, and
 *  e) contain a 'typedef <...> type ;' statement.
 */
template< class ListType > class TypeListMap ;

//----------------------------------------------------------------------

template<typename T> class TypeListMapValue ;

template<typename T>
class TypeListMapValue<T const &>
{
private:
  T const * x ;
public:
  typedef TypeListMapValue<T const &> SelfType ;
  typedef T const & const_reference_type ;
  typedef T const & reference_type ;

  TypeListMapValue() : x(0) {}
  TypeListMapValue( SelfType const & v ) : x( v.x ) {}
  explicit TypeListMapValue( T const & v ) : x( & v ) {}

  SelfType & operator = ( SelfType const & v ) { x = v.x ; return *this ; }
  SelfType & operator = ( T const & v )         { x = &v ; return *this ; }

  operator T const & () const { return *x ; }
  T const & get() const { return *x ; }
};

template<typename T>
class TypeListMapValue<T&>
{
private:
  T * x ;
public:
  typedef TypeListMapValue<T&> SelfType ;
  typedef T & const_reference_type ;
  typedef T & reference_type ;

  TypeListMapValue() : x(0) {}
  TypeListMapValue( const SelfType & v ) : x( v.x ) {}
  explicit TypeListMapValue( T & v ) : x( & v ) {}

  SelfType & operator = ( SelfType const & v ) { x = v.x ; return *this ; }
  SelfType & operator = ( T & v )              { x = &v ; return *this ; }

  operator T & () const { return *x ; }
  T & get() const { return *x ; }
};

template<typename T>
class TypeListMapValue<T const>
{
private:
  T x ;
public:
  typedef TypeListMapValue<T> SelfType ;
  typedef T const & const_reference_type ;
  typedef T const & reference_type ;

  TypeListMapValue() {}
  TypeListMapValue( SelfType const & v ) : x( v.x ) {}
  explicit TypeListMapValue( T const & v ) : x( v ) {}

  SelfType & operator = ( SelfType const & v ) { x = v.x ; return *this ; }
  SelfType & operator = ( T const & v )         { x = v ; return *this ; }

  operator T const & () const { return x ; }
  T const & get() const { return x ; }
};

template<typename T>
class TypeListMapValue
{
private:
  T x ;
public:
  typedef TypeListMapValue<T> SelfType ;
  typedef T const & const_reference_type ;
  typedef T & reference_type ;

  TypeListMapValue() {}
  TypeListMapValue( SelfType const & v ) : x( v.x ) {}
  explicit TypeListMapValue( T const & v ) : x( v ) {}

  SelfType & operator = ( SelfType const & v ) { x = v.x ; return *this ; }
  SelfType & operator = ( T const & v )         { x = v ; return *this ; }

  operator T const & () const { return x ; }
  T const & get() const { return x ; }
};

//----------------------------------------------------------------------

template<>
class TypeListMap<TypeListEnd> {};

template<typename Tail>
class TypeListMap< TypeList<TypeListEnd,Tail> > {};

template<class ListType>
class TypeListMap : public TypeListMap<typename ListType::TypeListTail>
{
private:
  template<typename U> friend class TypeListMap ;

  typedef typename ListType::TypeListTail  TailType ;
  typedef typename ListType::TypeListValue TagType ;
  typedef typename TagType::type type ;
  TypeListMapValue<type> m_value ;

public:

  typedef TypeListMap<ListType> SelfType ;

  //----------

  template<class Tag>
    typename TypeListMapValue<typename Tag::type>::const_reference_type
    get() const
    {
      typedef typename TypeListMember<ListType,Tag>::list_type MemberListType ;
      return ((TypeListMap<MemberListType> const &) *this).m_value.get();
    }

  template<class Tag>
    void
    set( typename TypeListMapValue<typename Tag::type>::const_reference_type v )
    {
      typedef typename TypeListMember<ListType,Tag>::list_type MemberListType ;
      ((TypeListMap<MemberListType> &) *this).m_value.operator=( v );
    }

  //----------

  TypeListMap<TailType> & operator <<
    ( typename TypeListMapValue<type>::const_reference_type v )
      { m_value = v ; return *this ; }

  TypeListMap<TailType> const & operator >>
    ( typename TypeListMapValue<type>::reference_type v ) const
      { v = m_value ; return *this ; }

  //----------

  void copy( TypeListMap<TypeListEnd> const & ) {}

  template<class ListB>
  void copy( TypeListMap<TypeList<TypeListEnd,ListB> > const & b ) {}

  template<class ListB>
  void copy( TypeListMap<ListB> const & b )
  {
    typedef typename ListB::TypeListValue TagB ;
    typedef typename ListB::TypeListTail  TailB ;
    this->template set<TagB>( b.template get<TagB>() );
    copy( (TypeListMap<TailB> const &) b );
  }

  //----------

  ~TypeListMap() {}

  TypeListMap() {}

  TypeListMap( const SelfType & m )
    : TypeListMap<TailType>( m ), m_value( m.m_value ) {}

  SelfType & operator = ( const SelfType & m )
    {
      TypeListMap<TailType>::operator=( m );
      m_value = m.m_value ;
      return *this ;
    }
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

}

#endif // STK_UTIL_UTIL_TypeListMap_h
