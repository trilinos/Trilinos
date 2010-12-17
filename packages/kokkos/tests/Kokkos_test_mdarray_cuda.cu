
#include <iostream>

#include <Kokkos_CudaMDArrayView.hpp>

#include <Kokkos_test_hex_grad.hpp>

typedef Kokkos::CudaMap array_map_type ;

//------------------------------------------------------------------------
/* TODO: Make part of 'run_kernel' implementation */

inline __device__
int work_stride() { return blockDim.x * gridDim.x ; }

inline __device__
int work_begin() { return threadIdx.x + blockDim.x * blockIdx.x ; }

__global__
void global_hex_simple_fill( Kokkos::MDArrayView<float,array_map_type> coord )
{
  const int stride = work_stride();
  const int nP = coord.dimension( coord.rank() - 1 );

  for ( int iP = work_begin() ; iP < nP ; iP += stride ) {
    kernel_hex_simple_fill( coord , iP );
  }
}

__global__
void global_hex_grad( Kokkos::MDArrayView<float,array_map_type> coord ,
                      Kokkos::MDArrayView<float,array_map_type> grad )
{
  const int stride = work_stride();
  const int nP = coord.dimension( coord.rank() - 1 );

  for ( int iP = work_begin() ; iP < nP ; iP += stride ) {
    kernel_hex_grad( coord , grad , iP );
  }
}

#define KOKKOS_RUN_1( function , arg0 ) \
{ \
  dim3 dimGrid( 256 , 1 , 1 ); \
  dim3 dimBlock( 256 , 1 , 1 ); \
  function <<< dimGrid , dimBlock >>> ( arg0 ); \
}

#define KOKKOS_RUN_2( function , arg0 , arg1 ) \
{ \
  dim3 dimGrid( 256 , 1 , 1 ); \
  dim3 dimBlock( 256 , 1 , 1 ); \
  function <<< dimGrid , dimBlock >>> ( arg0 , arg1 ); \
}


//------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  const int parallel_work_length = 1000000 ;

  array_map_type map( parallel_work_length );
  array_map_type map_2( parallel_work_length * 2 );

  Kokkos::MDArrayView< float , array_map_type > coord , grad ;

  coord = map.create_mdarray<float>( 3 , 8 );
  grad  = map.create_mdarray<float>( 3 , 8 );

  KOKKOS_RUN_1( global_hex_simple_fill , coord );
  KOKKOS_RUN_2( global_hex_grad , coord , grad );

  return 0 ;
}

