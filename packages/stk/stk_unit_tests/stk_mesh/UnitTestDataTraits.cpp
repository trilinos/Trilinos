// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

#include <stddef.h>                     // for NULL, size_t
#include <iostream>                     // for ostringstream, etc
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/DataTraits.hpp>  // for DataTraits, etc
#include <stk_mesh/base/DataTraitsClass.hpp>  // for data_traits, etc
#include <stk_mesh/base/DataTraitsEnum.hpp>  // for data_traits, etc
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <gtest/gtest.h>
#include <string>                       // for operator==, string, etc
#include <typeinfo>                     // for type_info
#include <vector>                       // for vector
#include "stk_util/parallel/ParallelComm.hpp"  // for CommAll




using stk::mesh::DataTraits;
using stk::mesh::data_traits;
using stk::CommAll;

//----------------------------------------------------------------------

namespace {

TEST(TestDataTraits, testVoid)
{
  // Test the DataTrait for void

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  const DataTraits & traits = data_traits<void>();

  ASSERT_TRUE(       traits.type_info        == typeid(void) );
  ASSERT_EQ( traits.size_of           , size_t(0) );
  ASSERT_EQ( traits.alignment_of      , size_t(0) );
  ASSERT_EQ( traits.stride_of         , size_t(0) );
  ASSERT_EQ( traits.is_void           , true );
  ASSERT_EQ( traits.is_integral       , false );
  ASSERT_EQ( traits.is_floating_point , false );
  ASSERT_EQ( traits.is_array          , false );
  ASSERT_EQ( traits.is_pointer        , false );
  ASSERT_EQ( traits.is_enum           , false );
  ASSERT_EQ( traits.is_class          , false );
  ASSERT_EQ( traits.is_pod            , false );
  ASSERT_EQ( traits.is_signed         , false );
  ASSERT_EQ( traits.is_unsigned       , false );

  ASSERT_TRUE( ! traits.remove_pointer   );
  ASSERT_TRUE( traits.enum_info.empty()  );
  ASSERT_TRUE( traits.class_info.empty() );

  ASSERT_THROW( traits.construct( NULL , 0 ) , std::runtime_error );
  ASSERT_THROW( traits.destroy( NULL , 0 ) , std::runtime_error );
  ASSERT_THROW( traits.copy( NULL , NULL , 0 ) , std::runtime_error );
  ASSERT_THROW( traits.sum( NULL , NULL , 0 ) , std::runtime_error );
  ASSERT_THROW( traits.max( NULL , NULL , 0 ) , std::runtime_error );
  ASSERT_THROW( traits.min( NULL , NULL , 0 ) , std::runtime_error );
  ASSERT_THROW( traits.bit_and( NULL, NULL, 0 ), std::runtime_error );
  ASSERT_THROW( traits.bit_or( NULL , NULL, 0 ), std::runtime_error );
  ASSERT_THROW( traits.bit_xor( NULL, NULL, 0 ), std::runtime_error );
  ASSERT_THROW( traits.print( std::cout, NULL, 0), std::runtime_error);

  CommAll comm( pm );
  ASSERT_THROW( traits.pack( comm.send_buffer(0) , NULL , 0 ), std::runtime_error );
  comm.allocate_buffers( 0 );
  comm.communicate();
  ASSERT_THROW( traits.unpack( comm.recv_buffer(0) , NULL , 0 ), std::runtime_error );
}

//----------------------------------------------------------------------

template< typename T , bool is_integral , bool is_signed >
void test_fundamental_type()
{
  // Test DataTrait for fundamental type T

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_rank = stk::parallel_machine_rank( pm );
  int p_size = stk::parallel_machine_size( pm );
  MPI_Barrier( pm );

  // Test data trait properties of type T
  const DataTraits & traits = data_traits<T>();
  ASSERT_TRUE(       traits.type_info        == typeid(T) );
  ASSERT_EQ( traits.size_of           , sizeof(T) );
  ASSERT_EQ( traits.alignment_of      , sizeof(T) );
  ASSERT_EQ( traits.stride_of         , sizeof(T) );
  ASSERT_EQ( traits.is_void           , false );
  ASSERT_EQ( traits.is_integral       , is_integral );
  ASSERT_EQ( traits.is_floating_point , ! is_integral );
  ASSERT_EQ( traits.is_array          , false );
  ASSERT_EQ( traits.is_pointer        , false );
  ASSERT_EQ( traits.is_enum           , false );
  ASSERT_EQ( traits.is_class          , false );
  ASSERT_EQ( traits.is_pod            , true );
  ASSERT_EQ( traits.is_signed         , is_signed );
  ASSERT_EQ( traits.is_unsigned       , is_integral && ! is_signed);

  ASSERT_TRUE( ! traits.remove_pointer   );
  ASSERT_TRUE( traits.enum_info.empty()  );
  ASSERT_TRUE( traits.class_info.empty() );

  const unsigned array_size = 3;
  const T a[array_size] = { T(1) , T(2) , T(4) };
  T b[array_size] ;

  // Test data trait basic operations on type T
  traits.construct( b , array_size );
  ASSERT_EQ( T(0) , b[0] );
  ASSERT_EQ( T(0) , b[1] );
  ASSERT_EQ( T(0) , b[2] );

  traits.copy( b , a , array_size );
  ASSERT_EQ( T(1) , b[0] );
  ASSERT_EQ( T(2) , b[1] );
  ASSERT_EQ( T(4) , b[2] );

  traits.sum( b , a , array_size );
  ASSERT_EQ( T(2) , b[0] );
  ASSERT_EQ( T(4) , b[1] );
  ASSERT_EQ( T(8) , b[2] );

  traits.min( b , a , array_size );
  ASSERT_EQ( T(1) , b[0] );
  ASSERT_EQ( T(2) , b[1] );
  ASSERT_EQ( T(4) , b[2] );

  traits.sum( b , a , array_size );
  traits.max( b , a , array_size );
  ASSERT_EQ( T(2) , b[0] );
  ASSERT_EQ( T(4) , b[1] );
  ASSERT_EQ( T(8) , b[2] );

  if ( is_integral ) {
    // Test integral-specific operations
    traits.bit_or( b , a , array_size );
    ASSERT_EQ( T(3) , b[0] );
    ASSERT_EQ( T(6) , b[1] );
    ASSERT_EQ( T(12) , b[2] );

    traits.bit_and( b , a , array_size );
    ASSERT_EQ( T(1) , b[0] );
    ASSERT_EQ( T(2) , b[1] );
    ASSERT_EQ( T(4) , b[2] );

    traits.bit_xor( b , a , array_size );
    ASSERT_EQ( T(0) , b[0] );
    ASSERT_EQ( T(0) , b[1] );
    ASSERT_EQ( T(0) , b[2] );
  }
  else {
    // Test unsupported operations
    ASSERT_THROW(traits.bit_or (b, a, array_size), std::runtime_error);
    ASSERT_THROW(traits.bit_and(b, a, array_size), std::runtime_error);
    ASSERT_THROW(traits.bit_xor(b, a, array_size), std::runtime_error);
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
      ASSERT_EQ( T(1) , b[0] );
      ASSERT_EQ( T(2) , b[1] );
      ASSERT_EQ( T(4) , b[2] );
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

TEST(TestDataTraits, testFundamental_bool)
{
  test_fundamental_type<char, true, true>();
}

TEST(TestDataTraits, testFundamental_unsignedchar)
{
  test_fundamental_type<unsigned char, true, false>();
}

TEST(TestDataTraits, testFundamental_short)
{
  test_fundamental_type<short, true, true>();
}

TEST(TestDataTraits, testFundamental_unsignedshort)
{
  test_fundamental_type<unsigned short, true, false>();
}

TEST(TestDataTraits, testFundamental_int)
{
  test_fundamental_type<int, true, true>();
}

TEST(TestDataTraits, testFundamental_unsignedint)
{
  test_fundamental_type<unsigned int, true, false>();
}

TEST(TestDataTraits, testFundamental_long)
{
  test_fundamental_type<long, true, true>();
}

TEST(TestDataTraits, testFundamental_unsignedlong)
{
  test_fundamental_type<unsigned long, true, false>();
}

TEST(TestDataTraits, testFundamental_float)
{
  test_fundamental_type<float, false, false>();
}

TEST(TestDataTraits, testFundamental_double)
{
  test_fundamental_type<double, false, false>();
}

//----------------------------------------------------------------------

template< typename T >
void test_fundamental_pointer()
{
  // Test DataTrait for fundamenter pointer type T*

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  // Test data trait properties of type T*
  const DataTraits & traits = data_traits<T*>();
  ASSERT_TRUE(       traits.type_info        == typeid(T*) );
  ASSERT_EQ( traits.size_of           , sizeof(T*) );
  ASSERT_EQ( traits.alignment_of      , sizeof(T*) );
  ASSERT_EQ( traits.stride_of         , sizeof(T*) );
  ASSERT_EQ( traits.is_void           , false );
  ASSERT_EQ( traits.is_integral       , false );
  ASSERT_EQ( traits.is_floating_point , false );
  ASSERT_EQ( traits.is_array          , false );
  ASSERT_EQ( traits.is_pointer        , true );
  ASSERT_EQ( traits.is_enum           , false );
  ASSERT_EQ( traits.is_class          , false );
  ASSERT_EQ( traits.is_pod            , false );
  ASSERT_EQ( traits.is_signed         , false );
  ASSERT_EQ( traits.is_unsigned       , false );

  ASSERT_TRUE( traits.remove_pointer == & data_traits<T>()  );
  ASSERT_TRUE( traits.enum_info.empty()  );
  ASSERT_TRUE( traits.class_info.empty() );

  const unsigned array_size = 3;
  T val[array_size] = { T(1) , T(2) , T(4) };
  T * const a[array_size] = { val , val + 1 , val + 2 };
  T * b[array_size] ;

  // Test data trait basic operations on type T*
  traits.construct( b , array_size );
  ASSERT_EQ( static_cast<T*>(NULL) , b[0] );
  ASSERT_EQ( static_cast<T*>(NULL) , b[1] );
  ASSERT_EQ( static_cast<T*>(NULL) , b[2] );

  traits.copy( b , a , array_size );
  ASSERT_EQ( val + 0 , b[0] );
  ASSERT_EQ( val + 1 , b[1] );
  ASSERT_EQ( val + 2 , b[2] );

  traits.destroy( b , array_size );
  ASSERT_EQ( static_cast<T*>(NULL) , b[0] );
  ASSERT_EQ( static_cast<T*>(NULL) , b[1] );
  ASSERT_EQ( static_cast<T*>(NULL) , b[2] );

  // Test unsupported operations
  ASSERT_THROW(traits.sum    (b, a, array_size),   std::runtime_error);
  ASSERT_THROW(traits.max    (b, a, array_size),   std::runtime_error);
  ASSERT_THROW(traits.min    (b, a, array_size),   std::runtime_error);
  ASSERT_THROW(traits.bit_and(b, a, array_size),   std::runtime_error);
  ASSERT_THROW(traits.bit_or( b, a, array_size),   std::runtime_error);
  ASSERT_THROW(traits.bit_xor(b, a, array_size),   std::runtime_error);
  ASSERT_THROW(traits.print  (std::cout, NULL, 0), std::runtime_error);

  CommAll comm( pm );
  ASSERT_THROW( traits.pack( comm.send_buffer(0) , a , array_size ), std::runtime_error );
  comm.allocate_buffers( 0 );
  comm.communicate();
  ASSERT_THROW( traits.unpack( comm.recv_buffer(0) , b , array_size ), std::runtime_error );
}

TEST(TestDataTraits, testFundamental_char_ptr)
{
  test_fundamental_pointer<char>();
}

TEST(TestDataTraits, testFundamental_unsignedchar_ptr)
{
  test_fundamental_pointer<unsigned char>();
}

TEST(TestDataTraits, testFundamental_short_ptr)
{
  test_fundamental_pointer<short>();
}

TEST(TestDataTraits, testFundamental_unsignedshort_ptr)
{
  test_fundamental_pointer<unsigned short>();
}

TEST(TestDataTraits, testFundamental_int_ptr)
{
  test_fundamental_pointer<int>();
}

TEST(TestDataTraits, testFundamental_unsignedint_ptr)
{
  test_fundamental_pointer<unsigned int>();
}

TEST(TestDataTraits, testFundamental_long_ptr)
{
  test_fundamental_pointer<long>();
}

TEST(TestDataTraits, testFundamental_unsignedlong_ptr)
{
  test_fundamental_pointer<unsigned long>();
}

TEST(TestDataTraits, testFundamental_float_ptr)
{
  test_fundamental_pointer<float>();
}

TEST(TestDataTraits, testFundamental_double_ptr)
{
  test_fundamental_pointer<double>();
}

}
//----------------------------------------------------------------------

#ifndef __PGI

enum EType { val_a = 'a' , val_b = 'b' , val_c = 'c' };

namespace stk {
namespace mesh {

// This enum will only work within the stk::mesh namespace
DATA_TRAITS_ENUM_3( EType , val_a , val_b , val_c )

}
}

namespace {

TEST(TestDataTraits, testEnum)
{
  // Test interaction of DataTraits with enums

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_rank = stk::parallel_machine_rank( pm );
  int p_size = stk::parallel_machine_size( pm );
  MPI_Barrier( pm );

  typedef EType T ;
  const DataTraits & traits = data_traits<T>();

  // Test data trait properties of enum type
  ASSERT_TRUE(       traits.type_info        == typeid(T) );
  ASSERT_EQ( traits.size_of           , sizeof(T) );
  ASSERT_EQ( traits.alignment_of      , sizeof(T) );
  ASSERT_EQ( traits.stride_of         , sizeof(T) );
  ASSERT_EQ( traits.is_integral       , false );
  ASSERT_EQ( traits.is_floating_point , false );
  ASSERT_EQ( traits.is_array          , false );
  ASSERT_EQ( traits.is_pointer        , false );
  ASSERT_EQ( traits.is_enum           , true );
  ASSERT_EQ( traits.is_class          , false );
  ASSERT_EQ( traits.is_pod            , true );
  ASSERT_EQ( traits.is_signed         , false );
  ASSERT_EQ( traits.is_unsigned       , false );

  ASSERT_TRUE( ! traits.remove_pointer   );
  ASSERT_TRUE( traits.class_info.empty() );

  ASSERT_EQ( traits.enum_info.size()   , size_t(3) );
  ASSERT_EQ( (traits.enum_info[0].name == "val_a"), true );
  ASSERT_EQ( (traits.enum_info[1].name  == "val_b"), true );
  ASSERT_EQ( (traits.enum_info[2].name  == "val_c"), true );
  ASSERT_EQ( traits.enum_info[0].value , long(val_a) );
  ASSERT_EQ( traits.enum_info[1].value , long(val_b) );
  ASSERT_EQ( traits.enum_info[2].value , long(val_c) );

  const unsigned array_size = 3;
  EType a[array_size] = { val_a , val_b , val_c };
  EType b[array_size] ;

  // Test data trait basic operations on enum type
  traits.construct( b , array_size );
  ASSERT_EQ( val_a , b[0] );
  ASSERT_EQ( val_a , b[1] );
  ASSERT_EQ( val_a , b[2] );

  traits.copy( b , a , array_size );
  ASSERT_EQ( a[0] , b[0] );
  ASSERT_EQ( a[1] , b[1] );
  ASSERT_EQ( a[2] , b[2] );

  b[0] = val_b ; b[1] = val_b ; b[2] = val_b ;

  traits.min( b , a , array_size );
  ASSERT_EQ( val_a , b[0] );
  ASSERT_EQ( val_b , b[1] );
  ASSERT_EQ( val_b , b[2] );

  b[0] = val_b ; b[1] = val_b ; b[2] = val_b ;

  traits.max( b , a , array_size );
  ASSERT_EQ( val_b , b[0] );
  ASSERT_EQ( val_b , b[1] );
  ASSERT_EQ( val_c , b[2] );

  // Test unsupported operations
  ASSERT_THROW(traits.sum    (b, a, array_size), std::runtime_error);
  ASSERT_THROW(traits.bit_and(b, a, array_size), std::runtime_error);
  ASSERT_THROW(traits.bit_or (b, a, array_size), std::runtime_error);
  ASSERT_THROW(traits.bit_xor(b, a, array_size), std::runtime_error);

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
      ASSERT_EQ( a[0] , b[0] );
      ASSERT_EQ( a[1] , b[1] );
      ASSERT_EQ( a[2] , b[2] );
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

namespace stk {
namespace mesh {

// This enum will only work within the stk::mesh namespace
DATA_TRAITS_POD_CLASS_3( Vec3 , x , y , z )

}
}

namespace {

TEST(TestDataTraits, testClass)
{
  // Test interaction of DataTraits with classes

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_rank = stk::parallel_machine_rank( pm );
  int p_size = stk::parallel_machine_size( pm );
  MPI_Barrier( pm );

  typedef Vec3 T ;
  const DataTraits & traits = data_traits<T>();

  // Test data trait properties of class type
  ASSERT_TRUE(       traits.type_info        == typeid(T) );
  ASSERT_EQ( traits.size_of           , sizeof(T) );
  ASSERT_EQ( traits.alignment_of      , sizeof(double) );
  ASSERT_EQ( traits.stride_of         , sizeof(T) );
  ASSERT_EQ( traits.is_integral       , false );
  ASSERT_EQ( traits.is_floating_point , false );
  ASSERT_EQ( traits.is_array          , false );
  ASSERT_EQ( traits.is_pointer        , false );
  ASSERT_EQ( traits.is_enum           , false );
  ASSERT_EQ( traits.is_class          , true );
  ASSERT_EQ( traits.is_pod            , true );
  ASSERT_EQ( traits.is_signed         , false );
  ASSERT_EQ( traits.is_unsigned       , false );

  ASSERT_TRUE( ! traits.remove_pointer   );
  ASSERT_TRUE( traits.enum_info.empty() );

  const Vec3 a = { 1.0 , 2.0 , 3.0 };
  const size_t dx = reinterpret_cast<const unsigned char *>( & a.x ) -
                    reinterpret_cast<const unsigned char *>( & a );
  const size_t dy = reinterpret_cast<const unsigned char *>( & a.y ) -
                    reinterpret_cast<const unsigned char *>( & a );
  const size_t dz = reinterpret_cast<const unsigned char *>( & a.z ) -
                    reinterpret_cast<const unsigned char *>( & a );

  ASSERT_EQ( traits.class_info.size() , size_t(3) );
  ASSERT_EQ( (traits.class_info[0].name == "x"), true );
  ASSERT_EQ( (traits.class_info[1].name == "y"), true );
  ASSERT_EQ( (traits.class_info[2].name == "z"), true );
  ASSERT_EQ( traits.class_info[0].traits, & data_traits<double>() );
  ASSERT_EQ( traits.class_info[1].traits, & data_traits<double>() );
  ASSERT_EQ( traits.class_info[2].traits, & data_traits<double>() );
  ASSERT_EQ( traits.class_info[0].offset, dx );
  ASSERT_EQ( traits.class_info[1].offset, dy );
  ASSERT_EQ( traits.class_info[2].offset, dz );

  // Test data trait basic operations on class type
  const unsigned array_size = 1;
  Vec3 b ;
  traits.construct( & b , array_size );
  traits.copy( & b , & a , array_size );
  ASSERT_EQ( b.x , a.x );
  ASSERT_EQ( b.y , a.y );
  ASSERT_EQ( b.z , a.z );

  // Test unsupport operations on class type
  ASSERT_THROW(traits.sum    (NULL, NULL, 0 ),     std::runtime_error);
  ASSERT_THROW(traits.max    (NULL, NULL, 0 ),     std::runtime_error);
  ASSERT_THROW(traits.min    (NULL, NULL, 0 ),     std::runtime_error);
  ASSERT_THROW(traits.bit_and(NULL, NULL, 0 ),     std::runtime_error);
  ASSERT_THROW(traits.bit_or (NULL, NULL, 0 ),     std::runtime_error);
  ASSERT_THROW(traits.bit_xor(NULL, NULL, 0 ),     std::runtime_error);
  ASSERT_THROW(traits.print  (std::cout, NULL, 0), std::runtime_error);

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
      ASSERT_EQ( a.x , b.x );
      ASSERT_EQ( a.y , b.y );
      ASSERT_EQ( a.z , b.z );
    }
  }

  traits.destroy( & b , array_size );
}

}

#endif
