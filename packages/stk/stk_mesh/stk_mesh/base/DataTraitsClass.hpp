/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/base/DataTraits.hpp>

#include <stk_util/environment/ReportHandler.hpp>

namespace stk {
namespace mesh {
 
//----------------------------------------------------------------------

namespace {

template< typename T >
class DataTraitsClassPOD : public DataTraits {
public:
  DataTraitsClassPOD( const char * name , std::size_t n )
  : DataTraits( typeid(T) , name , sizeof(T) , 1 )
  {
    is_pod = true ;
    is_class = true ;
    class_info.reserve( n );
  }

  void set_stride( const void * first , const void * second )
  {
    stride_of = reinterpret_cast<const unsigned char *>(second) -
                reinterpret_cast<const unsigned char *>(first);
  }

  void add_member( const char * n , const DataTraits & t ,
                   const void * base , const void * member )
  {
    const std::size_t i = class_info.size();
    const std::size_t d = reinterpret_cast<const unsigned char *>(member) -
                          reinterpret_cast<const unsigned char *>(base);
    class_info.resize( i + 1 );
    class_info[i].name.assign( n );
    class_info[i].traits = & t ;
    class_info[i].offset = d ;
    if ( alignment_of < t.alignment_of ) { alignment_of = t.alignment_of ; }
  }

  void construct( void * v , std::size_t n ) const
  {
    T * x = reinterpret_cast<T*>( v ); 
    T * const x_end = x + n ;  
    for ( ; x_end != x ; ++x ) { new(x) T(); }
  }
     
  void destroy( void * v , std::size_t n ) const
  {
    T * x = reinterpret_cast<T*>( v ); 
    T * const x_end = x + n ;  
    for ( ; x_end != x ; ++x ) { x->~T(); }
  }

  void copy( void * vx , const void * vy , std::size_t n ) const
  {
    const T * y = reinterpret_cast<const T*>( vy );
    T * x = reinterpret_cast<T*>( vx ); 
    T * const x_end = x + n ;
    while ( x_end != x ) { *x++ = *y++ ; };
  }

  void pack( CommBuffer & buf , const void * v , std::size_t n ) const
  {  
    const T * x = reinterpret_cast<const T*>( v );
    buf.pack<T>( x , n ); 
  }  
     
  void unpack( CommBuffer & buf , void * v , std::size_t n ) const
  {  
    T * x = reinterpret_cast<T*>( v );
    buf.unpack<T>( x , n ); 
  }  

  void print( std::ostream & s , const void * v , std::size_t n ) const
  { ThrowErrorMsg( "not supported" ); }

  void sum( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }
 
  void max( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }
 
  void min( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }
 
  void bit_and( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }
 
  void bit_or( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }
 
  void bit_xor( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }
};

}

//----------------------------------------------------------------------

#define DATA_TRAITS_POD_CLASS_2( C , M1 , M2 )      \
namespace {     \
class DataTraitsClass ## C : public DataTraitsClassPOD<C> {        \
public: \
  DataTraitsClass ## C () : DataTraitsClassPOD<C>( # C , 2 )       \
  {     \
    C tmp[1] ;     \
    set_stride( tmp , tmp + 1 ); \
    add_member( # M1 , data_traits( tmp->M1 ) , tmp , & tmp->M1 );      \
    add_member( # M2 , data_traits( tmp->M2 ) , tmp , & tmp->M2 );      \
  }     \
};      \
}       \
template<> const DataTraits & data_traits< C >()   \
{ static const DataTraitsClass ## C traits ; return traits ; }

//----------------------------------------------------------------------

#define DATA_TRAITS_POD_CLASS_3( C , M1 , M2 , M3 )      \
namespace {     \
class DataTraitsClass ## C : public DataTraitsClassPOD<C> {        \
public: \
  DataTraitsClass ## C () : DataTraitsClassPOD<C>( # C , 3 )       \
  {     \
    C tmp[1] ;     \
    set_stride( tmp , tmp + 1 ); \
    add_member( # M1 , data_traits( tmp->M1 ) , tmp , & tmp->M1 );      \
    add_member( # M2 , data_traits( tmp->M2 ) , tmp , & tmp->M2 );      \
    add_member( # M3 , data_traits( tmp->M3 ) , tmp , & tmp->M3 );      \
  }     \
};      \
}       \
template<> const DataTraits & data_traits< C >()   \
{ static const DataTraitsClass ## C traits ; return traits ; }

//----------------------------------------------------------------------

#define DATA_TRAITS_POD_CLASS_4( C , M1 , M2 , M3 , M4 )      \
namespace {     \
class DataTraitsClass ## C : public DataTraitsClassPOD<C> {        \
public: \
  DataTraitsClass ## C () : DataTraitsClassPOD<C>( # C , 4 )       \
  {     \
    C tmp[1] ;     \
    set_stride( tmp , tmp + 1 ); \
    add_member( # M1 , data_traits( tmp->M1 ) , tmp , & tmp->M1 );      \
    add_member( # M2 , data_traits( tmp->M2 ) , tmp , & tmp->M2 );      \
    add_member( # M3 , data_traits( tmp->M3 ) , tmp , & tmp->M3 );      \
    add_member( # M4 , data_traits( tmp->M4 ) , tmp , & tmp->M4 );      \
  }     \
};      \
}       \
template<> const DataTraits & data_traits< C >()   \
{ static const DataTraitsClass ## C traits ; return traits ; }

//----------------------------------------------------------------------

#define DATA_TRAITS_POD_CLASS_5( C , M1 , M2 , M3 , M4 , M5 )      \
namespace {     \
class DataTraitsClass ## C : public DataTraitsClassPOD<C> {        \
public: \
  DataTraitsClass ## C () : DataTraitsClassPOD<C>( # C , 5 )       \
  {     \
    C tmp[1] ;     \
    set_stride( tmp , tmp + 1 ); \
    add_member( # M1 , data_traits( tmp->M1 ) , tmp , & tmp->M1 );      \
    add_member( # M2 , data_traits( tmp->M2 ) , tmp , & tmp->M2 );      \
    add_member( # M3 , data_traits( tmp->M3 ) , tmp , & tmp->M3 );      \
    add_member( # M4 , data_traits( tmp->M4 ) , tmp , & tmp->M4 );      \
    add_member( # M5 , data_traits( tmp->M5 ) , tmp , & tmp->M5 );      \
  }     \
};      \
}       \
template<> const DataTraits & data_traits< C >()   \
{ static const DataTraitsClass ## C traits ; return traits ; }

//----------------------------------------------------------------------

#define DATA_TRAITS_POD_CLASS_6( C , M1 , M2 , M3 , M4 , M5 , M6 )      \
namespace {     \
class DataTraitsClass ## C : public DataTraitsClassPOD<C> {        \
public: \
  DataTraitsClass ## C () : DataTraitsClassPOD<C>( # C , 5 )       \
  {     \
    C tmp[1] ;     \
    set_stride( tmp , tmp + 1 ); \
    add_member( # M1 , data_traits( tmp->M1 ) , tmp , & tmp->M1 );      \
    add_member( # M2 , data_traits( tmp->M2 ) , tmp , & tmp->M2 );      \
    add_member( # M3 , data_traits( tmp->M3 ) , tmp , & tmp->M3 );      \
    add_member( # M4 , data_traits( tmp->M4 ) , tmp , & tmp->M4 );      \
    add_member( # M5 , data_traits( tmp->M5 ) , tmp , & tmp->M5 );      \
    add_member( # M6 , data_traits( tmp->M6 ) , tmp , & tmp->M6 );      \
  }     \
};      \
}       \
template<> const DataTraits & data_traits< C >()   \
{ static const DataTraitsClass ## C traits ; return traits ; }

//----------------------------------------------------------------------

}
}

