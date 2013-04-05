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
#include <KokkosArray_View.hpp>

namespace Explicit {

template< typename CoordScalarType ,
          typename GradScalarType >
KOKKOSARRAY_INLINE_FUNCTION
void grad( const CoordScalarType * const x ,
           const CoordScalarType * const z ,
                 GradScalarType  * const grad_y )
{
  const GradScalarType R42 = x[3] - x[1] ;
  const GradScalarType R52 = x[4] - x[1] ;
  const GradScalarType R54 = x[4] - x[3] ;

  const GradScalarType R63 = x[5] - x[2] ;
  const GradScalarType R83 = x[7] - x[2] ;
  const GradScalarType R86 = x[7] - x[5] ;

  const GradScalarType R31 = x[2] - x[0] ;
  const GradScalarType R61 = x[5] - x[0] ;
  const GradScalarType R74 = x[6] - x[3] ;

  const GradScalarType R72 = x[6] - x[1] ;
  const GradScalarType R75 = x[6] - x[4] ;
  const GradScalarType R81 = x[7] - x[0] ;

  const GradScalarType t1 = R63 + R54 ;
  const GradScalarType t2 = R61 + R74 ;
  const GradScalarType t3 = R72 + R81 ;

  const GradScalarType t4 = R86 + R42 ;
  const GradScalarType t5 = R83 + R52 ;
  const GradScalarType t6 = R75 + R31 ;

  //  Calculate Y gradient from X and Z data

  grad_y[0] = (z[1] *  t1) - (z[2] * R42) - (z[3] *  t5) + (z[4] *  t4) + (z[5] * R52) - (z[7] * R54);
  grad_y[1] = (z[2] *  t2) + (z[3] * R31) - (z[0] *  t1) - (z[5] *  t6) + (z[6] * R63) - (z[4] * R61);
  grad_y[2] = (z[3] *  t3) + (z[0] * R42) - (z[1] *  t2) - (z[6] *  t4) + (z[7] * R74) - (z[5] * R72);
  grad_y[3] = (z[0] *  t5) - (z[1] * R31) - (z[2] *  t3) + (z[7] *  t6) + (z[4] * R81) - (z[6] * R83);
  grad_y[4] = (z[5] *  t3) + (z[6] * R86) - (z[7] *  t2) - (z[0] *  t4) - (z[3] * R81) + (z[1] * R61);
  grad_y[5] = (z[6] *  t5) - (z[4] *  t3) - (z[7] * R75) + (z[1] *  t6) - (z[0] * R52) + (z[2] * R72);
  grad_y[6] = (z[7] *  t1) - (z[5] *  t5) - (z[4] * R86) + (z[2] *  t4) - (z[1] * R63) + (z[3] * R83);
  grad_y[7] = (z[4] *  t2) - (z[6] *  t1) + (z[5] * R75) - (z[3] *  t6) - (z[2] * R74) + (z[0] * R54);
}

template< typename CoordScalarType ,
          typename GradScalarType >
KOKKOSARRAY_INLINE_FUNCTION
void grad( const CoordScalarType x[] ,
           const CoordScalarType y[] ,
           const CoordScalarType z[] ,
                 GradScalarType grad_x[] ,
                 GradScalarType grad_y[] ,
                 GradScalarType grad_z[] )
{
  grad( x , z , grad_y );
  grad( z , y , grad_x );
  grad( y , x , grad_z );
}

template< typename CoordScalarType ,
          typename GradScalarType ,
          class    Device >
struct TestHexGrad
{
  typedef Device device_type ;

  typedef KokkosArray::View< CoordScalarType*[3][8] , device_type > coord_type ;
  typedef KokkosArray::View< GradScalarType *[3][8] , device_type > grad_type ;

  enum { NodeCount = 8 };

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const unsigned i ) const
  {

if ( 0 && i == 0 ) {
  printf("sizeof(CoordScalarType) = %d , sizeof(GradScalarType) = %d\n",
         (int) sizeof(CoordScalarType), (int) sizeof(GradScalarType) );
}

    CoordScalarType x[ NodeCount ];
    CoordScalarType y[ NodeCount ];
    CoordScalarType z[ NodeCount ];

    GradScalarType gx[ NodeCount ];
    GradScalarType gy[ NodeCount ];
    GradScalarType gz[ NodeCount ];

    for ( unsigned j = 0 ; j < NodeCount ; ++j ) {
      x[j] = coord(i,0,j);
      y[j] = coord(i,1,j);
      z[j] = coord(i,2,j);
    }

    grad( x, y, z, gx, gy, gz );

    for ( unsigned j = 0 ; j < NodeCount ; ++j ) {
      gradient(i,0,j) = gx[j] ;
      gradient(i,1,j) = gy[j] ;
      gradient(i,2,j) = gz[j] ;
    }
  }

  struct InitCoord {
    typedef Device device_type ;

    const coord_type coords ;

    InitCoord( const coord_type & c ) : coords( c ) {}

    KOKKOSARRAY_INLINE_FUNCTION
    void operator()( const unsigned ielem ) const
    {
      typedef Device device_type ;

      coords(ielem,0,0) = 0.;
      coords(ielem,1,0) = 0.;
      coords(ielem,2,0) = 0.;

      coords(ielem,0,1) = 1.;
      coords(ielem,1,1) = 0.;
      coords(ielem,2,1) = 0.;

      coords(ielem,0,2) = 1.;
      coords(ielem,1,2) = 1.;
      coords(ielem,2,2) = 0.;

      coords(ielem,0,3) = 0.;
      coords(ielem,1,3) = 1.;
      coords(ielem,2,3) = 0.;


      coords(ielem,0,4) = 0.;
      coords(ielem,1,4) = 0.;
      coords(ielem,2,4) = 1.;

      coords(ielem,0,5) = 1.;
      coords(ielem,1,5) = 0.;
      coords(ielem,2,5) = 1.;

      coords(ielem,0,6) = 1.;
      coords(ielem,1,6) = 1.;
      coords(ielem,2,6) = 1.;

      coords(ielem,0,7) = 0.;
      coords(ielem,1,7) = 1.;
      coords(ielem,2,7) = 1.;
    }
  };

  const coord_type coord ;
  const grad_type  gradient ;

  explicit
  TestHexGrad( size_t n )
  : coord("coord",n)
  , gradient("grad",n)
  {
    std::cout << "TestHexGrad::coord("
              << " " << coord.dimension_0() 
              << " " << coord.dimension_1()
              << " " << coord.dimension_2()
              << " " << coord.dimension_3()
              << " )" << std::endl ;
    KokkosArray::parallel_for(n,InitCoord(coord));
  }

  void apply() const
  {
    KokkosArray::parallel_for(coord.dimension_0(),*this);
  }
};

template< class Device >
void test( const std::string & label ,
           const size_t elem_count ,
           const size_t iter_count )
{
  KokkosArray::Impl::Timer timer ;

  double seconds_scalar ;
  double seconds_multi ;
  double seconds_array1 ;
  double seconds_array4 ;
  double seconds_array16 ;

  { // Loop 16 times:
    Explicit::TestHexGrad<double,float,Device> test_scalar( elem_count );

    timer.reset();

    for ( size_t i = 0 ; i < iter_count * 16 ; ++i ) {
      test_scalar.apply();
    }

    Device::fence();

    seconds_scalar = timer.seconds() / ( 16 * iter_count * elem_count );
  }

  { // 16 x elements
    Explicit::TestHexGrad<double,float,Device> test_multiple( elem_count * 16 );

    timer.reset();

    for ( size_t i = 0 ; i < iter_count ; ++i ) {
      test_multiple.apply();
    }

    Device::fence();

    seconds_multi = timer.seconds() / ( 16 * iter_count * elem_count );
  }

  { // 16 x elements with Array<1>
    typedef KokkosArray::Array<double,1> coord_scalar_type ;
    typedef KokkosArray::Array<float,1>  grad_scalar_type ;

    Explicit::TestHexGrad<coord_scalar_type,grad_scalar_type,Device>
      test_array( elem_count * 16 );

    timer.reset();

    for ( size_t i = 0 ; i < iter_count ; ++i ) {
      test_array.apply();
    }

    Device::fence();

    seconds_array1 = timer.seconds() / ( 16 * iter_count * elem_count );
  }

  { // 4 x elements with Array<4>
    typedef KokkosArray::Array<double,4> coord_scalar_type ;
    typedef KokkosArray::Array<float,4>  grad_scalar_type ;

    Explicit::TestHexGrad<coord_scalar_type,grad_scalar_type,Device>
      test_array( elem_count * 4 );

    timer.reset();

    for ( size_t i = 0 ; i < iter_count ; ++i ) {
      test_array.apply();
    }

    Device::fence();

    seconds_array4 = timer.seconds() / ( 16 * iter_count * elem_count );
  }

  { // 1 x elements with Array<16>
    typedef KokkosArray::Array<double,16> coord_scalar_type ;
    typedef KokkosArray::Array<float,16>  grad_scalar_type ;

    Explicit::TestHexGrad<coord_scalar_type,grad_scalar_type,Device> test_array( elem_count );

    timer.reset();

    for ( size_t i = 0 ; i < iter_count ; ++i ) {
      test_array.apply();
    }

    Device::fence();

    seconds_array16 = timer.seconds() / ( 16 * iter_count * elem_count );
  }

  std::cout << label
            << " scalar( " << seconds_scalar
            << " ) multi( " << seconds_multi << " )"
            << " ) array1( " << seconds_array1 << " )"
            << " ) array4( " << seconds_array4 << " )"
            << " ) array16( " << seconds_array16 << " )"
            << std::endl ;
}

} // namespace Explicit







