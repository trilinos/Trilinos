
#include <stddef.h>
#include <mpi.h>
#include <stdexcept>
#include <stk_mesh/base/DataTraits.hpp>
#include <stk_mesh/base/DataTraitsEnum.hpp>
#include <stk_mesh/base/DataTraitsClass.hpp>

#include <unit_tests/stk_utest_macros.hpp>

namespace stk {
namespace mesh {

class UnitTestDataTraits {
public:
  UnitTestDataTraits() {}

  void testVoid( bool );
  void testFundemental( bool );
  void testEnum( bool );
  void testClass( bool );
  void testPointer( bool );
};

}//namespace mesh
}//namespace stk

namespace {
STKUNIT_UNIT_TEST(TestDataTraits, testUnit)
{
  int mpi_rank = 0;
  int mpi_size = 0;
  
  STKUNIT_ASSERT_EQUAL(MPI_SUCCESS, MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank) );
  STKUNIT_ASSERT_EQUAL(MPI_SUCCESS, MPI_Comm_size(MPI_COMM_WORLD, &mpi_size) );
  STKUNIT_ASSERT(mpi_rank < mpi_size);

  bool testComm = mpi_size <= 1 ;

  stk::mesh::UnitTestDataTraits udt;
  if ( 0 == mpi_rank ) {
    udt.testVoid( testComm );
    udt.testFundemental( testComm );
    udt.testEnum( testComm );
    udt.testClass( testComm );
    udt.testPointer( testComm );
  }
}

}//namespace <anonymous>

namespace stk {
namespace mesh {
//----------------------------------------------------------------------

void UnitTestDataTraits::testVoid( bool testComm )
{
  const DataTraits & traits = data_traits<void>();

  STKUNIT_ASSERT(       traits.type_info        == typeid(void) );
  STKUNIT_ASSERT_EQUAL( traits.size_of           , size_t(0) );
  STKUNIT_ASSERT_EQUAL( traits.alignment_of      , size_t(0) );
  STKUNIT_ASSERT_EQUAL( traits.stride_of         , size_t(0) );
  STKUNIT_ASSERT_EQUAL( traits.is_void           , true );
  STKUNIT_ASSERT_EQUAL( traits.is_integral       , false );
  STKUNIT_ASSERT_EQUAL( traits.is_floating_point , false );
  STKUNIT_ASSERT_EQUAL( traits.is_array          , false );
  STKUNIT_ASSERT_EQUAL( traits.is_pointer        , false );
  STKUNIT_ASSERT_EQUAL( traits.is_enum           , false );
  STKUNIT_ASSERT_EQUAL( traits.is_class          , false );
  STKUNIT_ASSERT_EQUAL( traits.is_pod            , false );
  STKUNIT_ASSERT_EQUAL( traits.is_signed         , false );
  STKUNIT_ASSERT_EQUAL( traits.is_unsigned       , false );

  STKUNIT_ASSERT( ! traits.remove_pointer   );
  STKUNIT_ASSERT( traits.enum_info.empty()  );
  STKUNIT_ASSERT( traits.class_info.empty() );

  STKUNIT_ASSERT_THROW( traits.construct( NULL , 0 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.destroy( NULL , 0 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.copy( NULL , NULL , 0 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.sum( NULL , NULL , 0 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.max( NULL , NULL , 0 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.min( NULL , NULL , 0 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.bit_and( NULL, NULL, 0 ), std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.bit_or( NULL , NULL, 0 ), std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.bit_xor( NULL, NULL, 0 ), std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.print( std::cout, NULL, 0 ), std::runtime_error );

  if ( testComm ) {
    CommAll all( MPI_COMM_WORLD );
    STKUNIT_ASSERT_THROW( traits.pack( all.send_buffer(0) , NULL , 0 ), std::runtime_error );
    all.allocate_buffers( 0 );
    all.communicate();
    STKUNIT_ASSERT_THROW( traits.unpack( all.recv_buffer(0) , NULL , 0 ), std::runtime_error );
  }
}



//----------------------------------------------------------------------

template< typename T , bool is_integral , bool is_signed >
void test_fundemental_type( bool testComm )
{
  const DataTraits & traits = data_traits<T>();
  STKUNIT_ASSERT(       traits.type_info        == typeid(T) );
  STKUNIT_ASSERT_EQUAL( traits.size_of           , sizeof(T) );
  STKUNIT_ASSERT_EQUAL( traits.alignment_of      , sizeof(T) );
  STKUNIT_ASSERT_EQUAL( traits.stride_of         , sizeof(T) );
  STKUNIT_ASSERT_EQUAL( traits.is_void           , false );
  STKUNIT_ASSERT_EQUAL( traits.is_integral       , is_integral );
  STKUNIT_ASSERT_EQUAL( traits.is_floating_point , ! is_integral );
  STKUNIT_ASSERT_EQUAL( traits.is_array          , false );
  STKUNIT_ASSERT_EQUAL( traits.is_pointer        , false );
  STKUNIT_ASSERT_EQUAL( traits.is_enum           , false );
  STKUNIT_ASSERT_EQUAL( traits.is_class          , false );
  STKUNIT_ASSERT_EQUAL( traits.is_pod            , true );
  STKUNIT_ASSERT_EQUAL( traits.is_signed         , is_signed );
  STKUNIT_ASSERT_EQUAL( traits.is_unsigned       , is_integral && ! is_signed );

  STKUNIT_ASSERT( ! traits.remove_pointer   );
  STKUNIT_ASSERT( traits.enum_info.empty()  );
  STKUNIT_ASSERT( traits.class_info.empty() );


  const T a[3] = { T(1) , T(2) , T(4) };
  T b[3] ;

  traits.construct( b , 3 );
  STKUNIT_ASSERT_EQUAL( T(0) , b[0] );
  STKUNIT_ASSERT_EQUAL( T(0) , b[1] );
  STKUNIT_ASSERT_EQUAL( T(0) , b[2] );

  traits.copy( b , a , 3 );
  STKUNIT_ASSERT_EQUAL( T(1) , b[0] );
  STKUNIT_ASSERT_EQUAL( T(2) , b[1] );
  STKUNIT_ASSERT_EQUAL( T(4) , b[2] );

  traits.sum( b , a , 3 );
  STKUNIT_ASSERT_EQUAL( T(2) , b[0] );
  STKUNIT_ASSERT_EQUAL( T(4) , b[1] );
  STKUNIT_ASSERT_EQUAL( T(8) , b[2] );

  traits.min( b , a , 3 );
  STKUNIT_ASSERT_EQUAL( T(1) , b[0] );
  STKUNIT_ASSERT_EQUAL( T(2) , b[1] );
  STKUNIT_ASSERT_EQUAL( T(4) , b[2] );

  traits.sum( b , a , 3 );
  traits.max( b , a , 3 );
  STKUNIT_ASSERT_EQUAL( T(2) , b[0] );
  STKUNIT_ASSERT_EQUAL( T(4) , b[1] );
  STKUNIT_ASSERT_EQUAL( T(8) , b[2] );

  if ( is_integral ) {
    traits.bit_or( b , a , 3 );
    STKUNIT_ASSERT_EQUAL( T(3) , b[0] );
    STKUNIT_ASSERT_EQUAL( T(6) , b[1] );
    STKUNIT_ASSERT_EQUAL( T(12) , b[2] );

    traits.bit_and( b , a , 3 );
    STKUNIT_ASSERT_EQUAL( T(1) , b[0] );
    STKUNIT_ASSERT_EQUAL( T(2) , b[1] );
    STKUNIT_ASSERT_EQUAL( T(4) , b[2] );

    traits.bit_xor( b , a , 3 );
    STKUNIT_ASSERT_EQUAL( T(0) , b[0] );
    STKUNIT_ASSERT_EQUAL( T(0) , b[1] );
    STKUNIT_ASSERT_EQUAL( T(0) , b[2] );
  }
  else {
    STKUNIT_ASSERT_THROW( traits.bit_or(  b , a , 3 ), std::runtime_error );
    STKUNIT_ASSERT_THROW( traits.bit_and( b , a , 3 ), std::runtime_error );
    STKUNIT_ASSERT_THROW( traits.bit_xor( b , a , 3 ), std::runtime_error );
  }

  if ( testComm ) {
    traits.construct( b , 3 );
    CommAll all( MPI_COMM_WORLD );
    traits.pack( all.send_buffer(0) , a , 3 );
    all.allocate_buffers( 0 );
    traits.pack( all.send_buffer(0) , a , 3 );
    all.communicate();
    traits.unpack( all.recv_buffer(0) , b , 3 );
    STKUNIT_ASSERT_EQUAL( T(1) , b[0] );
    STKUNIT_ASSERT_EQUAL( T(2) , b[1] );
    STKUNIT_ASSERT_EQUAL( T(4) , b[2] );
  }

  std::cout << traits.name << " " ;
  traits.print( std::cout , a , 3 );
  std::cout << std::endl ;

  traits.destroy( b , 3 );
}


void UnitTestDataTraits::testFundemental( bool testComm )
{
  test_fundemental_type<char,           true, true >( testComm );
  test_fundemental_type<unsigned char,  true, false>( testComm );
  test_fundemental_type<short,          true, true >( testComm );
  test_fundemental_type<unsigned short, true, false>( testComm );
  test_fundemental_type<int,            true, true >( testComm );
  test_fundemental_type<unsigned int,   true, false>( testComm );
  test_fundemental_type<long,           true, true >( testComm );
  test_fundemental_type<unsigned long,  true, false>( testComm );
  test_fundemental_type<float,          false,false>( testComm );
  test_fundemental_type<double,         false,false>( testComm );
}

//----------------------------------------------------------------------

namespace {

template< typename T >
void test_fundemental_pointer( bool testComm )
{
  const DataTraits & traits = data_traits<T*>();
  STKUNIT_ASSERT(       traits.type_info        == typeid(T*) );
  STKUNIT_ASSERT_EQUAL( traits.size_of           , sizeof(T*) );
  STKUNIT_ASSERT_EQUAL( traits.alignment_of      , sizeof(T*) );
  STKUNIT_ASSERT_EQUAL( traits.stride_of         , sizeof(T*) );
  STKUNIT_ASSERT_EQUAL( traits.is_void           , false );
  STKUNIT_ASSERT_EQUAL( traits.is_integral       , false );
  STKUNIT_ASSERT_EQUAL( traits.is_floating_point , false );
  STKUNIT_ASSERT_EQUAL( traits.is_array          , false );
  STKUNIT_ASSERT_EQUAL( traits.is_pointer        , true );
  STKUNIT_ASSERT_EQUAL( traits.is_enum           , false );
  STKUNIT_ASSERT_EQUAL( traits.is_class          , false );
  STKUNIT_ASSERT_EQUAL( traits.is_pod            , false );
  STKUNIT_ASSERT_EQUAL( traits.is_signed         , false );
  STKUNIT_ASSERT_EQUAL( traits.is_unsigned       , false );

  STKUNIT_ASSERT( traits.remove_pointer == & data_traits<T>()  );
  STKUNIT_ASSERT( traits.enum_info.empty()  );
  STKUNIT_ASSERT( traits.class_info.empty() );


  T val[3] = { T(1) , T(2) , T(4) };
  T * const a[3] = { val , val + 1 , val + 2 };
  T * b[3] ;

  traits.construct( b , 3 );
  STKUNIT_ASSERT_EQUAL( reinterpret_cast<T*>(NULL) , b[0] );
  STKUNIT_ASSERT_EQUAL( reinterpret_cast<T*>(NULL) , b[1] );
  STKUNIT_ASSERT_EQUAL( reinterpret_cast<T*>(NULL) , b[2] );

  traits.copy( b , a , 3 );
  STKUNIT_ASSERT_EQUAL( val + 0 , b[0] );
  STKUNIT_ASSERT_EQUAL( val + 1 , b[1] );
  STKUNIT_ASSERT_EQUAL( val + 2 , b[2] );

  traits.destroy( b , 3 );
  STKUNIT_ASSERT_EQUAL( reinterpret_cast<T*>(NULL) , b[0] );
  STKUNIT_ASSERT_EQUAL( reinterpret_cast<T*>(NULL) , b[1] );
  STKUNIT_ASSERT_EQUAL( reinterpret_cast<T*>(NULL) , b[2] );

  STKUNIT_ASSERT_THROW( traits.sum( b , a , 3 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.max( b , a , 3 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.min( b , a , 3 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.bit_and( b, a, 3 ), std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.bit_or(  b, a, 3 ), std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.bit_xor( b, a, 3 ), std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.print( std::cout, NULL, 0 ), std::runtime_error );

  if ( testComm ) {
    CommAll all( MPI_COMM_WORLD );
    STKUNIT_ASSERT_THROW( traits.pack( all.send_buffer(0) , a , 3 ), std::runtime_error );
    all.allocate_buffers( 0 );
    all.communicate();
    STKUNIT_ASSERT_THROW( traits.unpack( all.recv_buffer(0) , b , 3 ), std::runtime_error );
  }
}

}


void UnitTestDataTraits::testPointer( bool testComm )
{
  test_fundemental_pointer<char          >( testComm );
  test_fundemental_pointer<unsigned char >( testComm );
  test_fundemental_pointer<short         >( testComm );
  test_fundemental_pointer<unsigned short>( testComm );
  test_fundemental_pointer<int           >( testComm );
  test_fundemental_pointer<unsigned int  >( testComm );
  test_fundemental_pointer<long          >( testComm );
  test_fundemental_pointer<unsigned long >( testComm );
  test_fundemental_pointer<float         >( testComm );
  test_fundemental_pointer<double        >( testComm );
}

//----------------------------------------------------------------------

enum EType { val_a = 'a' , val_b = 'b' , val_c = 'c' };

DATA_TRAITS_ENUM_3( EType , val_a , val_b , val_c )

void UnitTestDataTraits::testEnum( bool testComm )
{
  typedef EType T ;
  const DataTraits & traits = data_traits<T>();

  STKUNIT_ASSERT(       traits.type_info        == typeid(T) );
  STKUNIT_ASSERT_EQUAL( traits.size_of           , sizeof(T) );
  STKUNIT_ASSERT_EQUAL( traits.alignment_of      , sizeof(T) );
  STKUNIT_ASSERT_EQUAL( traits.stride_of         , sizeof(T) );
  STKUNIT_ASSERT_EQUAL( traits.is_integral       , false );
  STKUNIT_ASSERT_EQUAL( traits.is_floating_point , false );
  STKUNIT_ASSERT_EQUAL( traits.is_array          , false );
  STKUNIT_ASSERT_EQUAL( traits.is_pointer        , false );
  STKUNIT_ASSERT_EQUAL( traits.is_enum           , true );
  STKUNIT_ASSERT_EQUAL( traits.is_class          , false );
  STKUNIT_ASSERT_EQUAL( traits.is_pod            , true );
  STKUNIT_ASSERT_EQUAL( traits.is_signed         , false );
  STKUNIT_ASSERT_EQUAL( traits.is_unsigned       , false );

  STKUNIT_ASSERT( ! traits.remove_pointer   );
  STKUNIT_ASSERT( traits.class_info.empty() );

  STKUNIT_ASSERT_EQUAL( traits.enum_info.size()   , size_t(3) );
  STKUNIT_ASSERT_EQUAL( traits.enum_info[0].name  , std::string("val_a") );
  STKUNIT_ASSERT_EQUAL( traits.enum_info[1].name  , std::string("val_b") );
  STKUNIT_ASSERT_EQUAL( traits.enum_info[2].name  , std::string("val_c") );
  STKUNIT_ASSERT_EQUAL( traits.enum_info[0].value , long(val_a) );
  STKUNIT_ASSERT_EQUAL( traits.enum_info[1].value , long(val_b) );
  STKUNIT_ASSERT_EQUAL( traits.enum_info[2].value , long(val_c) );

  EType a[3] = { val_a , val_b , val_c };
  EType b[3] ;

  traits.construct( b , 3 );
  STKUNIT_ASSERT_EQUAL( val_a , b[0] );
  STKUNIT_ASSERT_EQUAL( val_a , b[1] );
  STKUNIT_ASSERT_EQUAL( val_a , b[2] );

  traits.copy( b , a , 3 );
  STKUNIT_ASSERT_EQUAL( a[0] , b[0] );
  STKUNIT_ASSERT_EQUAL( a[1] , b[1] );
  STKUNIT_ASSERT_EQUAL( a[2] , b[2] );

  b[0] = val_b ; b[1] = val_b ; b[2] = val_b ;

  traits.min( b , a , 3 );
  STKUNIT_ASSERT_EQUAL( val_a , b[0] );
  STKUNIT_ASSERT_EQUAL( val_b , b[1] );
  STKUNIT_ASSERT_EQUAL( val_b , b[2] );

  b[0] = val_b ; b[1] = val_b ; b[2] = val_b ;

  traits.max( b , a , 3 );
  STKUNIT_ASSERT_EQUAL( val_b , b[0] );
  STKUNIT_ASSERT_EQUAL( val_b , b[1] );
  STKUNIT_ASSERT_EQUAL( val_c , b[2] );

  STKUNIT_ASSERT_THROW( traits.sum( b , a , 3 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.bit_and( b, a, 3 ), std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.bit_or(  b, a, 3 ), std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.bit_xor( b, a, 3 ), std::runtime_error );

  if ( testComm ) {
    traits.construct( b , 3 );
    CommAll all( MPI_COMM_WORLD );
    traits.pack( all.send_buffer(0) , a , 3 );
    all.allocate_buffers( 0 );
    traits.pack( all.send_buffer(0) , a , 3 );
    all.communicate();
    traits.unpack( all.recv_buffer(0) , b , 3 );
    STKUNIT_ASSERT_EQUAL( a[0] , b[0] );
    STKUNIT_ASSERT_EQUAL( a[1] , b[1] );
    STKUNIT_ASSERT_EQUAL( a[2] , b[2] );
  }

  b[2] = static_cast<EType>( 'd' );
  std::cout << traits.name << " " ;
  traits.print( std::cout , b , 3 );
  std::cout << std::endl ;

  traits.destroy( b , 3 );
}

//----------------------------------------------------------------------

struct Vec3 { double x , y , z ; };

DATA_TRAITS_POD_CLASS_3( Vec3 , x , y , z )

void UnitTestDataTraits::testClass( bool testComm )
{
  typedef Vec3 T ;
  const DataTraits & traits = data_traits<T>();

  STKUNIT_ASSERT(       traits.type_info        == typeid(T) );
  STKUNIT_ASSERT_EQUAL( traits.size_of           , sizeof(T) );
  STKUNIT_ASSERT_EQUAL( traits.alignment_of      , sizeof(double) );
  STKUNIT_ASSERT_EQUAL( traits.stride_of         , sizeof(T) );
  STKUNIT_ASSERT_EQUAL( traits.is_integral       , false );
  STKUNIT_ASSERT_EQUAL( traits.is_floating_point , false );
  STKUNIT_ASSERT_EQUAL( traits.is_array          , false );
  STKUNIT_ASSERT_EQUAL( traits.is_pointer        , false );
  STKUNIT_ASSERT_EQUAL( traits.is_enum           , false );
  STKUNIT_ASSERT_EQUAL( traits.is_class          , true );
  STKUNIT_ASSERT_EQUAL( traits.is_pod            , true );
  STKUNIT_ASSERT_EQUAL( traits.is_signed         , false );
  STKUNIT_ASSERT_EQUAL( traits.is_unsigned       , false );

  STKUNIT_ASSERT( ! traits.remove_pointer   );
  STKUNIT_ASSERT( traits.enum_info.empty() );

  const Vec3 a = { 1.0 , 2.0 , 3.0 };
  const size_t dx = reinterpret_cast<const unsigned char *>( & a.x ) -
                    reinterpret_cast<const unsigned char *>( & a );
  const size_t dy = reinterpret_cast<const unsigned char *>( & a.y ) -
                    reinterpret_cast<const unsigned char *>( & a );
  const size_t dz = reinterpret_cast<const unsigned char *>( & a.z ) -
                    reinterpret_cast<const unsigned char *>( & a );

  STKUNIT_ASSERT_EQUAL( traits.class_info.size() , size_t(3) );
  STKUNIT_ASSERT_EQUAL( traits.class_info[0].name, std::string("x") );
  STKUNIT_ASSERT_EQUAL( traits.class_info[1].name, std::string("y") );
  STKUNIT_ASSERT_EQUAL( traits.class_info[2].name, std::string("z") );
  STKUNIT_ASSERT_EQUAL( traits.class_info[0].traits, & data_traits<double>() );
  STKUNIT_ASSERT_EQUAL( traits.class_info[1].traits, & data_traits<double>() );
  STKUNIT_ASSERT_EQUAL( traits.class_info[2].traits, & data_traits<double>() );
  STKUNIT_ASSERT_EQUAL( traits.class_info[0].offset, dx );
  STKUNIT_ASSERT_EQUAL( traits.class_info[1].offset, dy );
  STKUNIT_ASSERT_EQUAL( traits.class_info[2].offset, dz );

  Vec3 b ;
  traits.construct( & b , 1 );
  traits.copy( & b , & a , 1 );
  STKUNIT_ASSERT_EQUAL( b.x , a.x );
  STKUNIT_ASSERT_EQUAL( b.y , a.y );
  STKUNIT_ASSERT_EQUAL( b.z , a.z );

  STKUNIT_ASSERT_THROW( traits.sum( NULL , NULL , 0 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.max( NULL , NULL , 0 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.min( NULL , NULL , 0 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.bit_and( NULL, NULL, 0 ), std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.bit_or( NULL , NULL, 0 ), std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.bit_xor( NULL, NULL, 0 ), std::runtime_error );
  STKUNIT_ASSERT_THROW( traits.print( std::cout, NULL, 0 ), std::runtime_error );

  if ( testComm ) {
    traits.construct( & b , 1 );
    CommAll all( MPI_COMM_WORLD );
    traits.pack( all.send_buffer(0) , & a , 1 );
    all.allocate_buffers( 0 );
    traits.pack( all.send_buffer(0) , & a , 1 );
    all.communicate();
    traits.unpack( all.recv_buffer(0) , & b , 1 );
    STKUNIT_ASSERT_EQUAL( a.x , b.x );
    STKUNIT_ASSERT_EQUAL( a.y , b.y );
    STKUNIT_ASSERT_EQUAL( a.z , b.z );
  }

  traits.destroy( & b , 1 );

}

}
}

