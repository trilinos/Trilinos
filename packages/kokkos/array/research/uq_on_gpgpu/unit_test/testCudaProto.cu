/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/


#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>

#include <KokkosArray_Cuda.hpp>
#include <KokkosArray_Host.hpp>

#include <impl/KokkosArray_Timer.hpp>

#include <bbcsr.hpp>
#include <bsbcsr.hpp>

typedef KokkosArray::BigBlockCRSGraph<KokkosArray::Cuda>  graph_type ;
typedef KokkosArray::View< float[] , KokkosArray::Cuda >  matrix_type ;
typedef KokkosArray::View< double[] , KokkosArray::Cuda > vector_type ;

template< class VectorView >
void print_vector( const VectorView v )
{
  for ( int i = 0 ; i < v.length() ; ) {
    std::cout << std::setw(4) << i << " :" ;
    for ( int j = 0 ; j < 4 && i < v.length() ; ++i , ++j ) {
      std::cout << " " << std::setw(16) << v(i) ;
    }
    std::cout << std::endl ;
  }
}

struct PerformanceData {
  double seconds ;
  size_t nonzeros ;
  size_t flops ;
};


namespace diagonal_test { // Unit test with diagonal block graph

__global__
void initialize( const graph_type  graph ,
                 const matrix_type matrix ,
                 const vector_type input )
{
  const int tid = threadIdx.x + blockDim.x * blockIdx.x ;
  const int vector_begin = blockIdx.x * graph.block_size ;
  const int matrix_begin = blockIdx.x * graph.block_stride * graph.block_size ;

  if ( tid == 0 ) {
    for ( int i = 0 ; i <= gridDim.x ; ++i ) {
      graph.block_column_offset(i) = i ;
    }
    for ( int i = 0 ; i < gridDim.x ; ++i ) {
      graph.block_column_index(i) = i ;
    }
  }

  input( vector_begin + threadIdx.x ) = threadIdx.x ;

  if ( threadIdx.x % 2 ) {
    for ( int i = 0 ; i < graph.block_size ; ++i ) {

      matrix( matrix_begin + threadIdx.x + graph.block_stride * i ) = i * ( blockIdx.x + 1 );
    }
  }
}

void run( int block_count , int block_size )
{
  graph_type  graph ;
  matrix_type matrix ;
  vector_type input , output ;

  graph.block_system_size = block_count ;
  graph.block_size   = block_size ;
  graph.block_stride = block_size ; // No alignment consideration
  graph.block_maximum_columns = 1 ;

  graph.block_column_offset = graph_type::vector_type( block_count + 1 );
  graph.block_column_index  = graph_type::vector_type( block_count );

  const int vector_size = block_count * block_size ;
  const int matrix_size = block_count * block_size * block_size ;

  matrix = matrix_type( matrix_size );
  input  = vector_type( vector_size );
  output = vector_type( vector_size );

  dim3 grid_dim( block_count , 1 , 1 );
  dim3 block_dim( block_size , 1 , 1 );

  std::cout << "diagonal_test::run( "
            << block_count << " , "
            << block_size << " )" << std::endl ;

  initialize<<< grid_dim , block_dim >>>( graph , matrix , input );

  KokkosArray::multiply( graph , matrix , input , output );

  vector_type::HostMirror h_output = KokkosArray::mirror_create( output );
  vector_type::HostMirror h_input  = KokkosArray::mirror_create( input );
  matrix_type::HostMirror h_matrix = KokkosArray::mirror_create( matrix );

  KokkosArray::mirror_update( h_input , input );
  KokkosArray::mirror_update( h_output , output );
  KokkosArray::mirror_update( h_matrix , matrix );

  int sum = 0 ;
  for ( int i = 0 ; i < block_size ; ++i ) {
    sum += i * i ;
  }

  std::cout << "sum = " << sum << std::endl ;

//  std::cout << "matrix = " << std::endl ; print_vector( h_matrix );
//  std::cout << "input  = " << std::endl ; print_vector( h_input );
  std::cout << "output = " << std::endl ; print_vector( h_output );
}

}

namespace bbcsr_test {

// Unit test with checkerboard block graph
// starting at the diagonal and going right

__global__
void initialize( const graph_type  graph ,
                 const matrix_type matrix ,
                 const vector_type input )
{
  const int vector_begin = blockIdx.x * graph.block_size ;

  input( vector_begin + threadIdx.x ) = threadIdx.x + blockIdx.x * graph.block_size ;

  const int block_begin = graph.block_column_offset( blockIdx.x );
  const int block_end   = graph.block_column_offset( blockIdx.x + 1 );

  if ( threadIdx.x % 2 ) {
    for ( int i = block_begin ; i < block_end ; ++i ) {

      const int matrix_begin = i * graph.block_stride * graph.block_size ;

      for ( int k = 0 ; k < graph.block_size ; ++k ) {
        matrix( matrix_begin + threadIdx.x + graph.block_stride * k ) = blockIdx.x + 1 ;
      }
    }
  }
}

void run( const int block_count ,
          const int block_size ,
          const int number_iterations ,
          const bool check_answer ,
          PerformanceData & perf )
{
  graph_type  graph ;
  matrix_type matrix ;
  vector_type input , output ;

  int count = 0 ;
  int max_column_count = 0 ;
  for ( int i = 0 ; i < block_count ; ++i ) {
    int column_count = 0 ;
    for ( int j = i ; j < block_count ; j += 2 ) {
      ++column_count ;
    }
    if ( max_column_count < column_count ) max_column_count = column_count ;
    count += column_count ;
  }

  // Floating point multiply-add for each matrix entry
  perf.nonzeros = ((size_t) count ) *
                  ((size_t) block_size ) *
                  ((size_t) block_size );
  perf.flops = 2 * perf.nonzeros * ((size_t) number_iterations );

  graph.block_system_size = block_count ;
  graph.block_size   = block_size ;
  graph.block_stride = block_size ; // No alignment consideration
  graph.block_maximum_columns = max_column_count ;

  const int vector_size = graph.block_system_size * graph.block_size ;
  const int matrix_size = count * graph.block_stride * graph.block_size ;

  graph.block_column_offset = graph_type::vector_type( block_count + 1 );
  graph.block_column_index  = graph_type::vector_type( count );

  graph_type::vector_type::HostMirror h_column_offset =
    KokkosArray::mirror_create( graph.block_column_offset );
  graph_type::vector_type::HostMirror h_column_index  =
    KokkosArray::mirror_create( graph.block_column_index );

  count = 0 ;
  h_column_offset(0) = count ;
  for ( int i = 0 ; i < block_count ; ++i ) {
    for ( int j = i ; j < block_count ; j += 2 ) {
      h_column_index( count ) = j ;
      ++count ;
    }
    h_column_offset(i+1) = count ;
  }

  KokkosArray::mirror_update( graph.block_column_offset , h_column_offset );
  KokkosArray::mirror_update( graph.block_column_index  , h_column_index );

  matrix = matrix_type( matrix_size );
  input  = vector_type( vector_size );
  output = vector_type( vector_size );

  dim3 grid_dim( block_count , 1 , 1 );
  dim3 block_dim( block_size , 1 , 1 );

  initialize<<< grid_dim , block_dim >>>( graph , matrix , input );

  KokkosArray::Impl::Timer wall_clock ;

  for ( int i = 0 ; i < number_iterations ; ++i ) {
    KokkosArray::multiply( graph , matrix , input , output );
  }

  KokkosArray::Cuda::fence();

  perf.seconds = wall_clock.seconds();

  vector_type::HostMirror h_output = KokkosArray::mirror_create( output );
  vector_type::HostMirror h_input  = KokkosArray::mirror_create( input );
  matrix_type::HostMirror h_matrix = KokkosArray::mirror_create( matrix );

  KokkosArray::mirror_update( h_input , input );
  KokkosArray::mirror_update( h_output , output );
  KokkosArray::mirror_update( h_matrix , matrix );

//  std::cout << "matrix = " << std::endl ; print_vector( h_matrix );
//  std::cout << "input  = " << std::endl ; print_vector( h_input );
//  std::cout << "output = " << std::endl ; print_vector( h_output );

  if ( check_answer ) {

    for ( int i = 0 ; i < block_count ; ++i ) {
      const int output_begin = i * block_size ;
      const int matrix_coeff = i + 1 ;
      int sum = 0 ;
      for ( int j = i ; j < block_count ; j += 2 ) {
        int vector_begin = j * block_size ;
        for ( int k = 0 ; k < block_size ; ++k ) {
          sum += matrix_coeff * ( vector_begin + k );
        }
      }

      for ( int j = output_begin ; j < output_begin + block_size ; ++j ) {
        const int expected_value = j % 2 ? sum : 0 ;
        const int actual_value   = h_output(j);
        if ( expected_value != actual_value ) {
          std::cout << "check_test error at output(" << j << ") = "
                    << actual_value << " != " << expected_value 
                    << std::endl ;
          throw std::runtime_error(std::string("FAILED"));
        }
      }
    }
  }
}

}

//----------------------------------------------------------------------------

namespace bsbcsr_test {

__global__
void initialize(
  const KokkosArray::BigSymmetricBlockCSRGraph< KokkosArray::Cuda > graph ,
  const matrix_type matrix ,
  const vector_type input )
{
  const KokkosArray::DenseSymmetricMatrixDiagonalStorageMap
    storage_map( graph.block_length , graph.block_length );

  const int vector_begin = blockIdx.x * graph.block_length ;

  input( vector_begin + threadIdx.x ) =
    threadIdx.x + blockIdx.x * graph.block_length ;

  const int block_begin = graph.block_column_offset( blockIdx.x );
  const int block_end   = graph.block_column_offset( blockIdx.x + 1 );

  for ( int i = block_begin ; i < block_end ; ++i ) {

    const int matrix_begin = i * graph.diag_stride * graph.diag_count ;

    for ( int k = threadIdx.x ; k < graph.block_length ; k += 2 ) {
      const int offset = storage_map.offset( threadIdx.x , k );
      matrix( matrix_begin + offset ) = blockIdx.x + 1 ;
    }
  }
}

void run( const int block_count ,
          const int block_size ,
          const int number_iterations ,
          const bool check_answer ,
          PerformanceData & perf )
{
  KokkosArray::DenseSymmetricMatrixDiagonalStorageMap
    storage_map( block_size , block_size );

#if 0
  std::cout << "Location" << std::endl ;
  for ( int i = 0 ; i < block_size ; ++i ) {
    for ( int j = 0 ; j < block_size ; ) {
      std::cout << "M(" << i << "," << j << ") =" ;
      for ( int k = 0 ; k < 10 && j < block_size ; ++k , ++j ) {
        const int offset = storage_map.offset( i , j );
        std::cout << " " << offset ;
      }
      std::cout << std::endl ;
    }
  }
#endif

  KokkosArray::BigSymmetricBlockCSRGraph< KokkosArray::Cuda > graph ;
  matrix_type matrix ;
  vector_type input , output ;

  int count = 0 ;
  int max_column_count = 0 ;
  for ( int i = 0 ; i < block_count ; ++i ) {
    int column_count = 0 ;
    for ( int j = i ; j < block_count ; j += 2 ) {
      ++column_count ;
    }
    if ( max_column_count < column_count ) max_column_count = column_count ;
    count += column_count ;
  }

  // Floating point multiply-add for each matrix entry
  perf.nonzeros = ((size_t) count ) *
                  ((size_t) block_size ) *
                  ((size_t) block_size );
  perf.flops = 2 * perf.nonzeros * ((size_t) number_iterations );

  graph.block_system_size = block_count ;
  graph.block_length      = storage_map.m_length ;
  graph.diag_stride       = storage_map.m_stride ;
  graph.diag_count        = storage_map.m_diag_count ;

  const int vector_size = graph.block_system_size * graph.block_length ;
  const int matrix_size = count * graph.diag_stride * graph.diag_count ;

  graph.block_column_offset = graph_type::vector_type( block_count + 1 );
  graph.block_column_index  = graph_type::vector_type( count );

  graph_type::vector_type::HostMirror h_column_offset =
    KokkosArray::mirror_create( graph.block_column_offset );
  graph_type::vector_type::HostMirror h_column_index  =
    KokkosArray::mirror_create( graph.block_column_index );

  count = 0 ;
  h_column_offset(0) = count ;
  for ( int i = 0 ; i < block_count ; ++i ) {
    for ( int j = i ; j < block_count ; j += 2 ) {
      h_column_index( count ) = j ;
      ++count ;
    }
    h_column_offset(i+1) = count ;
  }

  KokkosArray::mirror_update( graph.block_column_offset , h_column_offset );
  KokkosArray::mirror_update( graph.block_column_index  , h_column_index );

  matrix = matrix_type( matrix_size );
  input  = vector_type( vector_size );
  output = vector_type( vector_size );

  dim3 grid_dim( block_count , 1 , 1 );
  dim3 block_dim( block_size , 1 , 1 );

  initialize<<< grid_dim , block_dim >>>( graph , matrix , input );

  KokkosArray::Impl::Timer wall_clock ;

  for ( int i = 0 ; i < number_iterations ; ++i ) {
    KokkosArray::multiply( graph , matrix , input , output );
  }

  KokkosArray::Cuda::fence();

  perf.seconds = wall_clock.seconds();

  vector_type::HostMirror h_output = KokkosArray::mirror_create( output );
  vector_type::HostMirror h_input  = KokkosArray::mirror_create( input );
  matrix_type::HostMirror h_matrix = KokkosArray::mirror_create( matrix );

  KokkosArray::mirror_update( h_input , input );
  KokkosArray::mirror_update( h_output , output );
  KokkosArray::mirror_update( h_matrix , matrix );

//  std::cout << "matrix = " << std::endl ; print_vector( h_matrix );
//  std::cout << "input  = " << std::endl ; print_vector( h_input );
//  std::cout << "output = " << std::endl ; print_vector( h_output );

  if ( check_answer ) {

    for ( int i = 0 ; i < block_count ; ++i ) {
      const int block_begin = h_column_offset(i);
      const int block_end   = h_column_offset(i+1);

      for ( int iy = 0 ; iy < block_size ; ++iy ) {
        const int global_y = iy + i * block_size ;
        double y = 0 ;

        for ( int j = block_begin ; j < block_end ; ++j ) {
          int vector_begin = h_column_index( j ) * block_size ;
          int matrix_begin = j * graph.diag_stride * graph.diag_count ;

          for ( int ix = 0 ; ix < block_size ; ++ix ) {
            const int offset = storage_map.offset( iy , ix );
            y += h_matrix( matrix_begin + offset ) *
                 h_input(  vector_begin + ix );
          }
        }

        const int expected_value = y ;
        const int actual_value   = h_output( global_y );

        if ( expected_value != actual_value ) {
          std::cout << "check_test error at output(" << global_y << ") = "
                    << actual_value << " != " << expected_value 
                    << std::endl ;
          throw std::runtime_error(std::string("FAILED"));
        }
      }
    }
  }
}

}

//----------------------------------------------------------------------------


int mainCuda()
{
  // diagonal_test::run( 1 , 32 );
  // diagonal_test::run( 1 , 48 );
  // diagonal_test::run( 10 , 32 );

  PerformanceData perf ;

  bbcsr_test::run( 8 , 32 , 1 , true , perf );
  bsbcsr_test::run( 8 , 32 , 1 , true , perf );

  std::cout << "\"Test\" , \"BlockSize\" , \"Nonzeros\" , \"Flops\" , \"Time\" , \"GFlop/sec\"" << std::endl ;

  const int block_size = 64 ;

  std::cout << std::endl ;

  for ( int i = 20 ; i <= 300 ; i += 20 ) {
    bbcsr_test::run( i , block_size , 1000 , false , perf );

    std::cout << "\"bbcsr_test(" << i << "," << block_size << ",1000) , "
              << block_size << " , " << perf.nonzeros << " , "
              << perf.flops << " , " << perf.seconds << " , "
              << ( (double) perf.flops / ( perf.seconds * 1e9 ) )
              << std::endl ;
  }

  std::cout << std::endl ;

  for ( int i = 20 ; i <= 300 ; i += 20 ) {
    bsbcsr_test::run( i , block_size , 1000 , false , perf );

    std::cout << "\"bsbcsr_test(" << i << "," << block_size << ",1000) , "
              << block_size << " , " << perf.nonzeros << " , "
              << perf.flops << " , " << perf.seconds << " , "
              << ( (double) perf.flops / ( perf.seconds * 1e9 ) )
              << std::endl ;
  }

  return 0 ;
}

