/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stddef.h>
#include <stdexcept>
#include <sstream>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/DataTraits.hpp>
#include <stk_mesh/base/DataTraitsEnum.hpp>
#include <stk_mesh/base/DataTraitsClass.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

using stk_classic::mesh::DataTraits;
using stk_classic::mesh::data_traits;
using stk_classic::CommAll;

//----------------------------------------------------------------------

namespace {

STKUNIT_UNIT_TEST(TestDataTraits, testVoid)
{
  // Test the DataTrait for void

  stk_classic::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

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
  STKUNIT_ASSERT_THROW( traits.print( std::cout, NULL, 0), std::runtime_error);

  CommAll comm( pm );
  STKUNIT_ASSERT_THROW( traits.pack( comm.send_buffer(0) , NULL , 0 ), std::runtime_error );
  comm.allocate_buffers( 0 );
  comm.communicate();
  STKUNIT_ASSERT_THROW( traits.unpack( comm.recv_buffer(0) , NULL , 0 ), std::runtime_error );
}

//----------------------------------------------------------------------

template< typename T , bool is_integral , bool is_signed >
void test_fundamental_type()
{
  // Test DataTrait for fundamental type T

  stk_classic::ParallelMachine pm = MPI_COMM_WORLD;
  int p_rank = stk_classic::parallel_machine_rank( pm );
  int p_size = stk_classic::parallel_machine_size( pm );
  MPI_Barrier( pm );

  // Test data trait properties of type T
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
  STKUNIT_ASSERT_EQUAL( traits.is_unsigned       , is_integral && ! is_signed);

  STKUNIT_ASSERT( ! traits.remove_pointer   );
  STKUNIT_ASSERT( traits.enum_info.empty()  );
  STKUNIT_ASSERT( traits.class_info.empty() );

  const unsigned array_size = 3;
  const T a[array_size] = { T(1) , T(2) , T(4) };
  T b[array_size] ;

  // Test data trait basic operations on type T
  traits.construct( b , array_size );
  STKUNIT_ASSERT_EQUAL( T(0) , b[0] );
  STKUNIT_ASSERT_EQUAL( T(0) , b[1] );
  STKUNIT_ASSERT_EQUAL( T(0) , b[2] );

  traits.copy( b , a , array_size );
  STKUNIT_ASSERT_EQUAL( T(1) , b[0] );
  STKUNIT_ASSERT_EQUAL( T(2) , b[1] );
  STKUNIT_ASSERT_EQUAL( T(4) , b[2] );

  traits.sum( b , a , array_size );
  STKUNIT_ASSERT_EQUAL( T(2) , b[0] );
  STKUNIT_ASSERT_EQUAL( T(4) , b[1] );
  STKUNIT_ASSERT_EQUAL( T(8) , b[2] );

  traits.min( b , a , array_size );
  STKUNIT_ASSERT_EQUAL( T(1) , b[0] );
  STKUNIT_ASSERT_EQUAL( T(2) , b[1] );
  STKUNIT_ASSERT_EQUAL( T(4) , b[2] );

  traits.sum( b , a , array_size );
  traits.max( b , a , array_size );
  STKUNIT_ASSERT_EQUAL( T(2) , b[0] );
  STKUNIT_ASSERT_EQUAL( T(4) , b[1] );
  STKUNIT_ASSERT_EQUAL( T(8) , b[2] );

  if ( is_integral ) {
    // Test integral-specific operations
    traits.bit_or( b , a , array_size );
    STKUNIT_ASSERT_EQUAL( T(3) , b[0] );
    STKUNIT_ASSERT_EQUAL( T(6) , b[1] );
    STKUNIT_ASSERT_EQUAL( T(12) , b[2] );

    traits.bit_and( b , a , array_size );
    STKUNIT_ASSERT_EQUAL( T(1) , b[0] );
    STKUNIT_ASSERT_EQUAL( T(2) , b[1] );
    STKUNIT_ASSERT_EQUAL( T(4) , b[2] );

    traits.bit_xor( b , a , array_size );
    STKUNIT_ASSERT_EQUAL( T(0) , b[0] );
    STKUNIT_ASSERT_EQUAL( T(0) , b[1] );
    STKUNIT_ASSERT_EQUAL( T(0) , b[2] );
  }
  else {
    // Test unsupported operations
    STKUNIT_ASSERT_THROW(traits.bit_or (b, a, array_size), std::runtime_error);
    STKUNIT_ASSERT_THROW(traits.bit_and(b, a, array_size), std::runtime_error);
    STKUNIT_ASSERT_THROW(traits.bit_xor(b, a, array_size), std::runtime_error);
  }

  // Test data trait pack/unpack (communication) of type T
  traits.construct( b , array_size );
  CommAll comm( pm );
  traits.pack( comm.send_buffer(0) , a , array_size );
  comm.allocate_buffers( 0 );
  traits.pack( comm.send_buffer(0) , a , array_size );
  comm.communicate();
  if (p_rank == 0) {
    for (int proc_id = 0; proc_id < p_size; ++proc_id) {
      traits.unpack( comm.recv_buffer(proc_id) , b , array_size );
      STKUNIT_ASSERT_EQUAL( T(1) , b[0] );
      STKUNIT_ASSERT_EQUAL( T(2) , b[1] );
      STKUNIT_ASSERT_EQUAL( T(4) , b[2] );
    }
  }

  // Test data trait print of type T
  std::ostringstream oss;
  oss << traits.name << " " ;
  traits.print( oss , a , array_size );
  oss << std::endl ;

  // Test data trait destruction (no-op in this case)
  traits.destroy( b , array_size );
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_bool)
{
  test_fundamental_type<char, true, true>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_unsignedchar)
{
  test_fundamental_type<unsigned char, true, false>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_short)
{
  test_fundamental_type<short, true, true>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_unsignedshort)
{
  test_fundamental_type<unsigned short, true, false>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_int)
{
  test_fundamental_type<int, true, true>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_unsignedint)
{
  test_fundamental_type<unsigned int, true, false>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_long)
{
  test_fundamental_type<long, true, true>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_unsignedlong)
{
  test_fundamental_type<unsigned long, true, false>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_float)
{
  test_fundamental_type<float, false, false>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_double)
{
  test_fundamental_type<double, false, false>();
}

//----------------------------------------------------------------------

template< typename T >
void test_fundamental_pointer()
{
  // Test DataTrait for fundamenter pointer type T*

  stk_classic::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  // Test data trait properties of type T*
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

  const unsigned array_size = 3;
  T val[array_size] = { T(1) , T(2) , T(4) };
  T * const a[array_size] = { val , val + 1 , val + 2 };
  T * b[array_size] ;

  // Test data trait basic operations on type T*
  traits.construct( b , array_size );
  STKUNIT_ASSERT_EQUAL( static_cast<T*>(NULL) , b[0] );
  STKUNIT_ASSERT_EQUAL( static_cast<T*>(NULL) , b[1] );
  STKUNIT_ASSERT_EQUAL( static_cast<T*>(NULL) , b[2] );

  traits.copy( b , a , array_size );
  STKUNIT_ASSERT_EQUAL( val + 0 , b[0] );
  STKUNIT_ASSERT_EQUAL( val + 1 , b[1] );
  STKUNIT_ASSERT_EQUAL( val + 2 , b[2] );

  traits.destroy( b , array_size );
  STKUNIT_ASSERT_EQUAL( static_cast<T*>(NULL) , b[0] );
  STKUNIT_ASSERT_EQUAL( static_cast<T*>(NULL) , b[1] );
  STKUNIT_ASSERT_EQUAL( static_cast<T*>(NULL) , b[2] );

  // Test unsupported operations
  STKUNIT_ASSERT_THROW(traits.sum    (b, a, array_size),   std::runtime_error);
  STKUNIT_ASSERT_THROW(traits.max    (b, a, array_size),   std::runtime_error);
  STKUNIT_ASSERT_THROW(traits.min    (b, a, array_size),   std::runtime_error);
  STKUNIT_ASSERT_THROW(traits.bit_and(b, a, array_size),   std::runtime_error);
  STKUNIT_ASSERT_THROW(traits.bit_or( b, a, array_size),   std::runtime_error);
  STKUNIT_ASSERT_THROW(traits.bit_xor(b, a, array_size),   std::runtime_error);
  STKUNIT_ASSERT_THROW(traits.print  (std::cout, NULL, 0), std::runtime_error);

  CommAll comm( pm );
  STKUNIT_ASSERT_THROW( traits.pack( comm.send_buffer(0) , a , array_size ), std::runtime_error );
  comm.allocate_buffers( 0 );
  comm.communicate();
  STKUNIT_ASSERT_THROW( traits.unpack( comm.recv_buffer(0) , b , array_size ), std::runtime_error );
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_char_ptr)
{
  test_fundamental_pointer<char>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_unsignedchar_ptr)
{
  test_fundamental_pointer<unsigned char>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_short_ptr)
{
  test_fundamental_pointer<short>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_unsignedshort_ptr)
{
  test_fundamental_pointer<unsigned short>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_int_ptr)
{
  test_fundamental_pointer<int>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_unsignedint_ptr)
{
  test_fundamental_pointer<unsigned int>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_long_ptr)
{
  test_fundamental_pointer<long>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_unsignedlong_ptr)
{
  test_fundamental_pointer<unsigned long>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_float_ptr)
{
  test_fundamental_pointer<float>();
}

STKUNIT_UNIT_TEST(TestDataTraits, testFundamental_double_ptr)
{
  test_fundamental_pointer<double>();
}

}
//----------------------------------------------------------------------

#ifndef __PGI

enum EType { val_a = 'a' , val_b = 'b' , val_c = 'c' };

namespace stk_classic {
namespace mesh {

// This enum will only work within the stk_classic::mesh namespace
DATA_TRAITS_ENUM_3( EType , val_a , val_b , val_c )

}
}

namespace {

STKUNIT_UNIT_TEST(TestDataTraits, testEnum)
{
  // Test interaction of DataTraits with enums

  stk_classic::ParallelMachine pm = MPI_COMM_WORLD;
  int p_rank = stk_classic::parallel_machine_rank( pm );
  int p_size = stk_classic::parallel_machine_size( pm );
  MPI_Barrier( pm );

  typedef EType T ;
  const DataTraits & traits = data_traits<T>();

  // Test data trait properties of enum type
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
  STKUNIT_ASSERT_EQUAL( (traits.enum_info[0].name == "val_a"), true );
  STKUNIT_ASSERT_EQUAL( (traits.enum_info[1].name  == "val_b"), true );
  STKUNIT_ASSERT_EQUAL( (traits.enum_info[2].name  == "val_c"), true );
  STKUNIT_ASSERT_EQUAL( traits.enum_info[0].value , long(val_a) );
  STKUNIT_ASSERT_EQUAL( traits.enum_info[1].value , long(val_b) );
  STKUNIT_ASSERT_EQUAL( traits.enum_info[2].value , long(val_c) );

  const unsigned array_size = 3;
  EType a[array_size] = { val_a , val_b , val_c };
  EType b[array_size] ;

  // Test data trait basic operations on enum type
  traits.construct( b , array_size );
  STKUNIT_ASSERT_EQUAL( val_a , b[0] );
  STKUNIT_ASSERT_EQUAL( val_a , b[1] );
  STKUNIT_ASSERT_EQUAL( val_a , b[2] );

  traits.copy( b , a , array_size );
  STKUNIT_ASSERT_EQUAL( a[0] , b[0] );
  STKUNIT_ASSERT_EQUAL( a[1] , b[1] );
  STKUNIT_ASSERT_EQUAL( a[2] , b[2] );

  b[0] = val_b ; b[1] = val_b ; b[2] = val_b ;

  traits.min( b , a , array_size );
  STKUNIT_ASSERT_EQUAL( val_a , b[0] );
  STKUNIT_ASSERT_EQUAL( val_b , b[1] );
  STKUNIT_ASSERT_EQUAL( val_b , b[2] );

  b[0] = val_b ; b[1] = val_b ; b[2] = val_b ;

  traits.max( b , a , array_size );
  STKUNIT_ASSERT_EQUAL( val_b , b[0] );
  STKUNIT_ASSERT_EQUAL( val_b , b[1] );
  STKUNIT_ASSERT_EQUAL( val_c , b[2] );

  // Test unsupported operations
  STKUNIT_ASSERT_THROW(traits.sum    (b, a, array_size), std::runtime_error);
  STKUNIT_ASSERT_THROW(traits.bit_and(b, a, array_size), std::runtime_error);
  STKUNIT_ASSERT_THROW(traits.bit_or (b, a, array_size), std::runtime_error);
  STKUNIT_ASSERT_THROW(traits.bit_xor(b, a, array_size), std::runtime_error);

  // Test pack/unpack (communication) of enum type
  traits.construct( b , array_size );
  CommAll comm( pm );
  traits.pack( comm.send_buffer(0) , a , array_size );
  comm.allocate_buffers( 0 );
  traits.pack( comm.send_buffer(0) , a , array_size );
  comm.communicate();
  if (p_rank == 0) {
    for (int proc_id = 0; proc_id < p_size; ++proc_id) {
      traits.unpack( comm.recv_buffer(proc_id) , b , array_size );
      STKUNIT_ASSERT_EQUAL( a[0] , b[0] );
      STKUNIT_ASSERT_EQUAL( a[1] , b[1] );
      STKUNIT_ASSERT_EQUAL( a[2] , b[2] );
    }
  }

  // Test printing of enum type
  std::ostringstream oss;
  b[2] = static_cast<EType>( 'd' );
  oss << traits.name << " " ;
  traits.print( oss , b , array_size );
  oss << std::endl ;

  // Test destruction of enum type (no-op in this case)
  traits.destroy( b , array_size );
}

}

//----------------------------------------------------------------------

struct Vec3 { double x , y , z ; };

namespace stk_classic {
namespace mesh {

// This enum will only work within the stk_classic::mesh namespace
DATA_TRAITS_POD_CLASS_3( Vec3 , x , y , z )

}
}

namespace {

STKUNIT_UNIT_TEST(TestDataTraits, testClass)
{
  // Test interaction of DataTraits with classes

  stk_classic::ParallelMachine pm = MPI_COMM_WORLD;
  int p_rank = stk_classic::parallel_machine_rank( pm );
  int p_size = stk_classic::parallel_machine_size( pm );
  MPI_Barrier( pm );

  typedef Vec3 T ;
  const DataTraits & traits = data_traits<T>();

  // Test data trait properties of class type
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
  STKUNIT_ASSERT_EQUAL( (traits.class_info[0].name == "x"), true );
  STKUNIT_ASSERT_EQUAL( (traits.class_info[1].name == "y"), true );
  STKUNIT_ASSERT_EQUAL( (traits.class_info[2].name == "z"), true );
  STKUNIT_ASSERT_EQUAL( traits.class_info[0].traits, & data_traits<double>() );
  STKUNIT_ASSERT_EQUAL( traits.class_info[1].traits, & data_traits<double>() );
  STKUNIT_ASSERT_EQUAL( traits.class_info[2].traits, & data_traits<double>() );
  STKUNIT_ASSERT_EQUAL( traits.class_info[0].offset, dx );
  STKUNIT_ASSERT_EQUAL( traits.class_info[1].offset, dy );
  STKUNIT_ASSERT_EQUAL( traits.class_info[2].offset, dz );

  // Test data trait basic operations on class type
  const unsigned array_size = 1;
  Vec3 b ;
  traits.construct( & b , array_size );
  traits.copy( & b , & a , array_size );
  STKUNIT_ASSERT_EQUAL( b.x , a.x );
  STKUNIT_ASSERT_EQUAL( b.y , a.y );
  STKUNIT_ASSERT_EQUAL( b.z , a.z );

  // Test unsupport operations on class type
  STKUNIT_ASSERT_THROW(traits.sum    (NULL, NULL, 0 ),     std::runtime_error);
  STKUNIT_ASSERT_THROW(traits.max    (NULL, NULL, 0 ),     std::runtime_error);
  STKUNIT_ASSERT_THROW(traits.min    (NULL, NULL, 0 ),     std::runtime_error);
  STKUNIT_ASSERT_THROW(traits.bit_and(NULL, NULL, 0 ),     std::runtime_error);
  STKUNIT_ASSERT_THROW(traits.bit_or (NULL, NULL, 0 ),     std::runtime_error);
  STKUNIT_ASSERT_THROW(traits.bit_xor(NULL, NULL, 0 ),     std::runtime_error);
  STKUNIT_ASSERT_THROW(traits.print  (std::cout, NULL, 0), std::runtime_error);

  // Test data trait pack/unpack (communication) of class type
  traits.construct( & b , array_size );
  CommAll comm( pm );
  traits.pack( comm.send_buffer(0) , & a , array_size );
  comm.allocate_buffers( 0 );
  traits.pack( comm.send_buffer(0) , & a , array_size );
  comm.communicate();
  if (p_rank == 0) {
    for (int proc_id = 0; proc_id < p_size; ++proc_id) {
      traits.unpack( comm.recv_buffer(proc_id) , & b , array_size );
      STKUNIT_ASSERT_EQUAL( a.x , b.x );
      STKUNIT_ASSERT_EQUAL( a.y , b.y );
      STKUNIT_ASSERT_EQUAL( a.z , b.z );
    }
  }

  traits.destroy( & b , array_size );
}

}

#endif
