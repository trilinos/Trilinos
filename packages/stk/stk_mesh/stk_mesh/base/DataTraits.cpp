/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <cstddef>
#include <stdexcept>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_mesh/base/DataTraits.hpp>
#include <stk_mesh/base/DataTraitsEnum.hpp>
#include <stk_mesh/base/DataTraitsClass.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

namespace {

std::size_t stride( std::size_t size , std::size_t align )
{
  if ( align && size % align ) { size += align - size % align ; }
  return size ;
}

}

DataTraits::~DataTraits() {}

DataTraits::DataTraits( const std::type_info & arg_type ,
                        const char * const     arg_name ,
                        std::size_t            arg_size ,
                        std::size_t            arg_align )
  : type_info(         arg_type ),
    size_of(           arg_size ),
    is_void(           false ),
    is_integral(       false ),
    is_floating_point( false ),
    is_pointer(        false ),
    is_enum(           false ),
    is_class(          false ),
    is_pod(            false ),
    is_signed(         false ),
    is_unsigned(       false ),
    alignment_of(      arg_align ),
    stride_of(         stride( arg_size , arg_align ) ),
    remove_pointer(    NULL ),
    name(              arg_name ),
    enum_info(),
    class_info()
{}

DataTraits::DataTraits( const std::type_info & arg_type ,
                        const DataTraits     & arg_traits )
  : type_info(          arg_type ),
    size_of(            sizeof(void*) ),
    is_void(            false ),
    is_integral(        false ),
    is_floating_point(  false ),
    is_pointer(         true ),
    is_enum(            false ),
    is_class(           false ),
    is_pod(             false ),
    is_signed(          false ),
    is_unsigned(        false ),
    alignment_of(       sizeof(void*) ),
    stride_of(          sizeof(void*) ),
    remove_pointer( & arg_traits ),
    name(),
    enum_info(),
    class_info()
{
  name.assign( arg_traits.name ).append("*");
}

//----------------------------------------------------------------------

namespace {

class DataTraitsVoid : public DataTraits {
public:

  DataTraitsVoid()
    : DataTraits( typeid(void) , "void" , 0 , 0 )
    { is_void = true ; }

  void construct( void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  void destroy( void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  void pack( CommBuffer & , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  void unpack( CommBuffer & , void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  void print( std::ostream & , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  void copy( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  void sum( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  void max( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  void min( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  virtual void bit_and( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  virtual void bit_or( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  virtual void bit_xor( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }
};

}

template<> const DataTraits & data_traits<void>()
{ static const DataTraitsVoid traits ; return traits ; }

//----------------------------------------------------------------------

namespace {

template< typename A , typename B >
struct IsSameType { enum { value = false }; };

template< typename A >
struct IsSameType<A,A> { enum { value = true }; };


template< typename T >
class DataTraitsCommon : public DataTraits {
public:

  explicit DataTraitsCommon( const char * arg_name )
    : DataTraits( typeid(T) , arg_name , sizeof(T) , sizeof(T) )
  {
    is_pod            = true ;

    is_integral       = IsSameType<T,char>::value ||
                        IsSameType<T,unsigned char>::value ||
                        IsSameType<T,short>::value ||
                        IsSameType<T,unsigned short>::value ||
                        IsSameType<T,int>::value ||
                        IsSameType<T,unsigned int>::value ||
                        IsSameType<T,long>::value ||
                        IsSameType<T,unsigned long>::value ;

    is_signed         = IsSameType<T,char>::value ||
                        IsSameType<T,short>::value ||
                        IsSameType<T,int>::value ||
                        IsSameType<T,long>::value ;

    is_unsigned       = IsSameType<T,unsigned char>::value ||
                        IsSameType<T,unsigned short>::value ||
                        IsSameType<T,unsigned int>::value ||
                        IsSameType<T,unsigned long>::value ;

    is_floating_point = IsSameType<T,double>::value ||
                        IsSameType<T,float>::value ;
  }

  void construct( void * v , std::size_t n ) const
  {
    T * x = reinterpret_cast<T*>( v );
    T * const x_end = x + n ;
    while ( x_end != x ) { *x++ = 0 ; }
  }

  void destroy( void * v , std::size_t n ) const {}

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

  void copy( void * vx , const void * vy , std::size_t n ) const
  {
    const T * y = reinterpret_cast<const T*>( vy );
    T * x = reinterpret_cast<T*>( vx );
    T * const x_end = x + n ;
    while ( x_end != x ) { *x++ = *y++ ; };
  }

  void sum( void * vx , const void * vy , std::size_t n ) const
  {
    const T * y = reinterpret_cast<const T*>( vy );
    T * x = reinterpret_cast<T*>( vx );
    T * const x_end = x + n ;
    while ( x_end != x ) { *x++ += *y++ ; };
  }

  virtual void print( std::ostream & s , const void * v , std::size_t n ) const
  {
    if ( n ) {
      const T * x = reinterpret_cast<const T*>( v );
      const T * const x_end = x + n ;
      s << *x++ ;
      while ( x_end != x ) { s << " " << *x++ ; }
    }
  }

  virtual void max( void * vx , const void * vy , std::size_t n ) const
  { ThrowErrorMsg( "not supported" ); }

  virtual void min( void * vx , const void * vy , std::size_t n ) const
  { ThrowErrorMsg( "not supported" ); }

  virtual void bit_and( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  virtual void bit_or( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  virtual void bit_xor( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }
};

template< typename T >
class DataTraitsNumeric : public DataTraitsCommon<T> {
public:

  explicit DataTraitsNumeric( const char * arg_name )
    : DataTraitsCommon<T>( arg_name )  {}

  virtual void max( void * vx , const void * vy , std::size_t n ) const
  {
    const T * y = reinterpret_cast<const T*>( vy );
    T * x = reinterpret_cast<T*>( vx );
    T * const x_end = x + n ;
    for ( ; x_end != x ; ++x , ++y ) { if ( *x < *y ) { *x = *y ; } }
  }

  virtual void min( void * vx , const void * vy , std::size_t n ) const
  {
    const T * y = reinterpret_cast<const T*>( vy );
    T * x = reinterpret_cast<T*>( vx );
    T * const x_end = x + n ;
    for ( ; x_end != x ; ++x , ++y ) { if ( *x > *y ) { *x = *y ; } }
  }
};

template< typename T >
class DataTraitsComplex : public DataTraitsCommon<T> {
public:

  explicit DataTraitsComplex( const char * arg_name )
    : DataTraitsCommon<T>( arg_name ) {}
};

template< typename T >
class DataTraitsIntegral : public DataTraitsNumeric<T> {
public:
  DataTraitsIntegral( const char * name ) : DataTraitsNumeric<T>( name ) {}

  virtual void bit_and( void * vx , const void * vy , std::size_t n ) const
  {
    const T * y = reinterpret_cast<const T*>( vy );
    T * x = reinterpret_cast<T*>( vx );
    T * const x_end = x + n ;
    while ( x_end != x ) { *x++ &= *y++ ; }
  }

  virtual void bit_or( void * vx , const void * vy , std::size_t n ) const
  {
    const T * y = reinterpret_cast<const T*>( vy );
    T * x = reinterpret_cast<T*>( vx );
    T * const x_end = x + n ;
    while ( x_end != x ) { *x++ |= *y++ ; }
  }

  virtual void bit_xor( void * vx , const void * vy , std::size_t n ) const
  {
    const T * y = reinterpret_cast<const T*>( vy );
    T * x = reinterpret_cast<T*>( vx );
    T * const x_end = x + n ;
    while ( x_end != x ) { *x++ ^= *y++ ; }
  }
};

class DataTraitsChar : public DataTraitsIntegral<char> {
public:
  DataTraitsChar() : DataTraitsIntegral<char>( "char" ) {}

  virtual void print( std::ostream & s , const void * v , std::size_t n ) const
  {
    if ( n ) {
      const char * x = reinterpret_cast<const char*>( v );
      const char * const x_end = x + n ;
      s << int(*x++) ;
      while ( x_end != x ) { s << " " << int(*x++) ; }
    }
  }
};

class DataTraitsUnsignedChar : public DataTraitsIntegral<unsigned char> {
public:
  DataTraitsUnsignedChar()
    : DataTraitsIntegral<unsigned char>( "unsigned char" ) {}

  virtual void print( std::ostream & s , const void * v , std::size_t n ) const
  {
    if ( n ) {
      const unsigned char * x = reinterpret_cast<const unsigned char*>( v );
      const unsigned char * const x_end = x + n ;
      s << unsigned(*x++) ;
      while ( x_end != x ) { s << " " << unsigned(*x++) ; }
    }
  }
};

}

#define DATA_TRAITS_NUMERIC( T )        \
template<>      \
const DataTraits & data_traits<T>()     \
{ static const DataTraitsNumeric<T> traits( #T ); return traits ; }

#define DATA_TRAITS_COMPLEX( T )        \
template<>      \
const DataTraits & data_traits<T>()     \
{ static const DataTraitsComplex<T> traits( #T ); return traits ; }

#define DATA_TRAITS_INTEGRAL( T )        \
template<>      \
const DataTraits & data_traits<T>()     \
{ static const DataTraitsIntegral<T> traits( #T ); return traits ; }

template<>
const DataTraits & data_traits<char>()
{ static const DataTraitsChar traits ; return traits ; }

template<>
const DataTraits & data_traits<unsigned char>()
{ static const DataTraitsUnsignedChar traits ; return traits ; }

DATA_TRAITS_INTEGRAL( short )
DATA_TRAITS_INTEGRAL( unsigned short )
DATA_TRAITS_INTEGRAL( int )
DATA_TRAITS_INTEGRAL( unsigned int )
DATA_TRAITS_INTEGRAL( long )
DATA_TRAITS_INTEGRAL( unsigned long )
DATA_TRAITS_NUMERIC( float )
DATA_TRAITS_NUMERIC( double )
DATA_TRAITS_COMPLEX( std::complex<float> ) // TODO: Probably not right
DATA_TRAITS_COMPLEX( std::complex<double> ) // TODO: Probably not right

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

template< typename T >
class DataTraitsPointerToFundamental : public DataTraits {
public:

  DataTraitsPointerToFundamental()
    : DataTraits( typeid(T*) , data_traits<T>() ) {}

  void construct( void * v , std::size_t n ) const
  {
    void ** x = reinterpret_cast<void**>(v);
    void ** const x_end = x + n ;
    while ( x_end != x ) { *x++ = NULL ; }
  }

  void destroy( void * v , std::size_t n ) const
  {
    void ** x = reinterpret_cast<void**>(v);
    void ** const x_end = x + n ;
    while ( x_end != x ) { *x++ = NULL ; }
  }

  void copy( void * vx , const void * vy , std::size_t n ) const
  {
    void * const * y = reinterpret_cast<void* const *>(vy);
    void ** x = reinterpret_cast<void**>(vx);
    void ** const x_end = x + n ;
    while ( x_end != x ) { *x++ = *y++ ; }
  }

  void pack( CommBuffer & , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  void unpack( CommBuffer & , void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  void print( std::ostream & , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  void sum( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  void max( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  void min( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  virtual void bit_and( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  virtual void bit_or( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }

  virtual void bit_xor( void * , const void * , std::size_t ) const
  { ThrowErrorMsg( "not supported" ); }
};

}

#define DATA_TRAITS_POINTER( T )        \
template<> const DataTraits & data_traits<T*>()  \
{ static const DataTraitsPointerToFundamental<T> traits ; return traits ; }

DATA_TRAITS_POINTER( char )
DATA_TRAITS_POINTER( unsigned char )
DATA_TRAITS_POINTER( short )
DATA_TRAITS_POINTER( unsigned short )
DATA_TRAITS_POINTER( int )
DATA_TRAITS_POINTER( unsigned int )
DATA_TRAITS_POINTER( long )
DATA_TRAITS_POINTER( unsigned long )
DATA_TRAITS_POINTER( float )
DATA_TRAITS_POINTER( double )
DATA_TRAITS_POINTER( void )
DATA_TRAITS_POINTER( std::complex<float> )
DATA_TRAITS_POINTER( std::complex<double> )

//----------------------------------------------------------------------

}
}



