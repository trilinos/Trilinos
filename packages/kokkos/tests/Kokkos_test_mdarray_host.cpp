
#include <iostream>
#include <Kokkos_HostMDArrayView.hpp>

#include <Kokkos_test_hex_grad.hpp>

typedef Kokkos::HostMap array_map_type ;

//------------------------------------------------------------------------
/* TODO: Make part of 'run_kernel' implementation */

inline
int work_stride() { return 1 ; }

inline
int work_begin() { return 0 ; }

void global_hex_simple_fill( Kokkos::MDArrayView<float,array_map_type> coord )
{
  const int stride = work_stride();
  const int nP = coord.dimension( coord.rank() - 1 );

  for ( int iP = work_begin() ; iP < nP ; iP += stride ) {
    kernel_hex_simple_fill( coord , iP );
  }
}

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
  function( arg0 ); \
}

#define KOKKOS_RUN_2( function , arg0 , arg1 ) \
{ \
  function( arg0 , arg1 ); \
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

  {
    Kokkos::MDArrayView< float , array_map_type > tmp1 = coord ;
    Kokkos::MDArrayView< float , array_map_type > tmp2 = tmp1 ;
  }

  KOKKOS_RUN_1( global_hex_simple_fill , coord );
  KOKKOS_RUN_2( global_hex_grad , coord , grad );

  return 0 ;
}

