
#include <stk_mesh/base/DataTraits.hpp>
 
//----------------------------------------------------------------------

namespace stk {
namespace mesh {

template< typename EnumType > class DataTraitsEnum ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
 
namespace {

template< typename T >
class DataTraitsEnum : public DataTraits {
public:
  DataTraitsEnum( const char * name , std::size_t n )
  : DataTraits( typeid(T) , name , sizeof(T) , sizeof(T) )
  {
    is_pod  = true ;
    is_enum = true ;
    enum_info.reserve( n );
  }

  void add_member( const char * n , T v )
  {
    const std::size_t i = enum_info.size();
    enum_info.resize( i + 1 );
    enum_info[i].name.assign( n );
    enum_info[i].value = static_cast<long>( v );
  }

  void construct( void * v , std::size_t n ) const
  {
    const T init = static_cast<T>( enum_info.front().value );
    T * x = reinterpret_cast<T*>(v);
    T * const x_end = x + n ;
    while ( x_end != x ) { *x++ = init ; }
  }

  void destroy( void * v , std::size_t n ) const {}

  void copy( void * vx , const void * vy , std::size_t n ) const
  {
    const T * y = reinterpret_cast<const T*>(vy);
    T * x = reinterpret_cast<T*>(vx);
    T * const x_end = x + n ;
    while ( x_end != x ) { *x++ = *y++ ; }
  }

  void max( void * vx , const void * vy , std::size_t n ) const
  {
    const T * y = reinterpret_cast<const T*>(vy);
    T * x = reinterpret_cast<T*>(vx);
    T * const x_end = x + n ;
    for ( ; x_end != x ; ++x , ++y ) { if ( *x < *y ) { *x = *y ; } }
  }

  void min( void * vx , const void * vy , std::size_t n ) const
  {
    const T * y = reinterpret_cast<const T*>(vy);
    T * x = reinterpret_cast<T*>(vx);
    T * const x_end = x + n ;
    for ( ; x_end != x ; ++x , ++y ) { if ( *x > *y ) { *x = *y ; } }
  }

  void print_one( std::ostream & s , T v ) const
  {
    std::vector<EnumMember>::const_iterator i = enum_info.begin();
    for ( ; i != enum_info.end() && i->value != v ; ++i );
    if ( i != enum_info.end() ) {
      s << i->name ;
    }
    else {
      s << name << "( " << static_cast<long>( v ) << " VALUE_NOT_VALID )" ;
    }
  }

  void print( std::ostream & s , const void * v , std::size_t n ) const
  {
    if ( n ) {
      const T * x = reinterpret_cast<const T*>(v);
      const T * const x_end = x + n ;
      print_one( s , *x++ );
      while ( x_end != x ) { s << " " ; print_one( s , *x++ ); }
    }
  }

  void pack( CommBuffer & buf , const void * v , std::size_t n ) const
  {
    const T * x = reinterpret_cast<const T*>(v);
    buf.pack<T>( x , n );
  }

  void unpack( CommBuffer & buf , void * v , std::size_t n ) const
  {
    T * x = reinterpret_cast<T*>(v);
    buf.unpack<T>( x , n );
  }

  void sum( void * , const void * , std::size_t ) const
  { throw_not_supported( "sum" ); }

  void bit_and( void * , const void * , std::size_t ) const
  { throw_not_supported( "sum" ); }

  void bit_or( void * , const void * , std::size_t ) const
  { throw_not_supported( "sum" ); }

  void bit_xor( void * , const void * , std::size_t ) const
  { throw_not_supported( "sum" ); }
};

}

//----------------------------------------------------------------------

#define DATA_TRAITS_ENUM_1( T , V1 )    \
namespace {     \
class DataTraitsEnum ## T : public DataTraitsEnum<T> {      \
public: \
  DataTraitsEnum ## T () : DataTraitsEnum<T>( # T , 1 ) \
  { add_member( # V1 , V1 ); }  \
};      \
}       \
template<> const DataTraits & data_traits< T >()   \
{ static const DataTraitsEnum ## T traits ; return traits ; }

//----------------------------------------------------------------------

#define DATA_TRAITS_ENUM_2( T , V1 , V2 )    \
namespace {     \
class DataTraitsEnum ## T : public DataTraitsEnum<T> {      \
public: \
  DataTraitsEnum ## T () : DataTraitsEnum<T>( # T , 2 ) \
  { \
    add_member( # V1 , V1 ); \
    add_member( # V2 , V2 ); \
  }  \
};      \
}       \
template<> const DataTraits & data_traits< T >()   \
{ static const DataTraitsEnum ## T traits ; return traits ; }

//----------------------------------------------------------------------

#define DATA_TRAITS_ENUM_3( T , V1 , V2 , V3 )    \
namespace {     \
class DataTraitsEnum ## T : public DataTraitsEnum<T> {      \
public: \
  DataTraitsEnum ## T () : DataTraitsEnum<T>( # T , 3 ) \
  { \
    add_member( # V1 , V1 ); \
    add_member( # V2 , V2 ); \
    add_member( # V3 , V3 ); \
  }  \
};      \
}       \
template<> const DataTraits & data_traits< T >()   \
{ static const DataTraitsEnum ## T traits ; return traits ; }

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

