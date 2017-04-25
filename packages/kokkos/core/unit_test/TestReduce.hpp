/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
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

#include <stdexcept>
#include <sstream>
#include <iostream>
#include <limits>

#include <Kokkos_Core.hpp>

namespace Test {

template< typename ScalarType, class DeviceType >
class ReduceFunctor
{
public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;

  struct value_type {
    ScalarType value[3];
  };

  const size_type nwork;

  ReduceFunctor( const size_type & arg_nwork )
    : nwork( arg_nwork ) {}

  ReduceFunctor( const ReduceFunctor & rhs )
    : nwork( rhs.nwork ) {}

/*
  KOKKOS_INLINE_FUNCTION
  void init( value_type & dst ) const
  {
    dst.value[0] = 0;
    dst.value[1] = 0;
    dst.value[2] = 0;
  }
*/

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & dst,
             const volatile value_type & src ) const
  {
    dst.value[0] += src.value[0];
    dst.value[1] += src.value[1];
    dst.value[2] += src.value[2];
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type iwork, value_type & dst ) const
  {
    dst.value[0] += 1;
    dst.value[1] += iwork + 1;
    dst.value[2] += nwork - iwork;
  }
};

template< class DeviceType >
class ReduceFunctorFinal : public ReduceFunctor< long, DeviceType > {
public:
  typedef typename ReduceFunctor< long, DeviceType >::value_type value_type;

  ReduceFunctorFinal( const size_t n )
    : ReduceFunctor< long, DeviceType >( n ) {}

  KOKKOS_INLINE_FUNCTION
  void final( value_type & dst ) const
  {
    dst.value[0] = -dst.value[0];
    dst.value[1] = -dst.value[1];
    dst.value[2] = -dst.value[2];
  }
};

template< typename ScalarType, class DeviceType >
class RuntimeReduceFunctor
{
public:
  // Required for functor:
  typedef DeviceType  execution_space;
  typedef ScalarType  value_type[];
  const unsigned      value_count;

  // Unit test details:

  typedef typename execution_space::size_type size_type;

  const size_type     nwork;

  RuntimeReduceFunctor( const size_type arg_nwork,
                        const size_type arg_count )
    : value_count( arg_count )
    , nwork( arg_nwork ) {}

  KOKKOS_INLINE_FUNCTION
  void init( ScalarType dst[] ) const
  {
    for ( unsigned i = 0; i < value_count; ++i ) dst[i] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile ScalarType dst[],
             const volatile ScalarType src[] ) const
  {
    for ( unsigned i = 0; i < value_count; ++i ) dst[i] += src[i];
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type iwork, ScalarType dst[] ) const
  {
    const size_type tmp[3] = { 1, iwork + 1, nwork - iwork };

    for ( size_type i = 0; i < value_count; ++i ) {
      dst[i] += tmp[ i % 3 ];
    }
  }
};

template< typename ScalarType, class DeviceType >
class RuntimeReduceMinMax
{
public:
  // Required for functor:
  typedef DeviceType  execution_space;
  typedef ScalarType  value_type[];
  const unsigned      value_count;

  // Unit test details:

  typedef typename execution_space::size_type size_type;

  const size_type     nwork;
  const ScalarType    amin;
  const ScalarType    amax;

  RuntimeReduceMinMax( const size_type arg_nwork,
                       const size_type arg_count )
    : value_count( arg_count )
    , nwork( arg_nwork )
    , amin( std::numeric_limits< ScalarType >::min() )
    , amax( std::numeric_limits< ScalarType >::max() )
    {}

  KOKKOS_INLINE_FUNCTION
  void init( ScalarType dst[] ) const
  {
    for ( unsigned i = 0; i < value_count; ++i ) {
      dst[i] = i % 2 ? amax : amin;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile ScalarType dst[],
             const volatile ScalarType src[] ) const
  {
    for ( unsigned i = 0; i < value_count; ++i ) {
      dst[i] = i % 2 ? ( dst[i] < src[i] ? dst[i] : src[i] )  // min
                     : ( dst[i] > src[i] ? dst[i] : src[i] ); // max
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type iwork, ScalarType dst[] ) const
  {
    const ScalarType tmp[2] = { ScalarType( iwork + 1 )
                              , ScalarType( nwork - iwork ) };

    for ( size_type i = 0; i < value_count; ++i ) {
      dst[i] = i % 2 ? ( dst[i] < tmp[i % 2] ? dst[i] : tmp[i % 2] )
                     : ( dst[i] > tmp[i % 2] ? dst[i] : tmp[i % 2] );
    }
  }
};

template< class DeviceType >
class RuntimeReduceFunctorFinal : public RuntimeReduceFunctor< long, DeviceType > {
public:
  typedef RuntimeReduceFunctor< long, DeviceType > base_type;
  typedef typename base_type::value_type value_type;
  typedef long scalar_type;

  RuntimeReduceFunctorFinal( const size_t theNwork, const size_t count )
    : base_type( theNwork, count ) {}

  KOKKOS_INLINE_FUNCTION
  void final( value_type dst ) const
  {
    for ( unsigned i = 0; i < base_type::value_count; ++i ) {
      dst[i] = -dst[i];
    }
  }
};

} // namespace Test

namespace {

template< typename ScalarType, class DeviceType >
class TestReduce
{
public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;

  TestReduce( const size_type & nwork )
  {
    run_test( nwork );
    run_test_final( nwork );
  }

  void run_test( const size_type & nwork )
  {
    typedef Test::ReduceFunctor< ScalarType, execution_space > functor_type;
    typedef typename functor_type::value_type value_type;

    enum { Count = 3 };
    enum { Repeat = 100 };

    value_type result[ Repeat ];

    const unsigned long nw   = nwork;
    const unsigned long nsum = nw % 2 ? nw * ( ( nw + 1 ) / 2 )
                                      : ( nw / 2 ) * ( nw + 1 );

    for ( unsigned i = 0; i < Repeat; ++i ) {
      Kokkos::parallel_reduce( nwork, functor_type( nwork ), result[i] );
    }

    for ( unsigned i = 0; i < Repeat; ++i ) {
      for ( unsigned j = 0; j < Count; ++j ) {
        const unsigned long correct = 0 == j % 3 ? nw : nsum;
        ASSERT_EQ( (ScalarType) correct, result[i].value[j] );
      }
    }
  }

  void run_test_final( const size_type & nwork )
  {
    typedef Test::ReduceFunctorFinal< execution_space > functor_type;
    typedef typename functor_type::value_type value_type;

    enum { Count = 3 };
    enum { Repeat = 100 };

    value_type result[ Repeat ];

    const unsigned long nw   = nwork;
    const unsigned long nsum = nw % 2 ? nw * ( ( nw + 1 ) / 2 )
                                      : ( nw / 2 ) * ( nw + 1 );

    for ( unsigned i = 0; i < Repeat; ++i ) {
      if ( i % 2 == 0 ) {
        Kokkos::parallel_reduce( nwork, functor_type( nwork ), result[i] );
      }
      else {
        Kokkos::parallel_reduce( "Reduce", nwork, functor_type( nwork ), result[i] );
      }
    }

    for ( unsigned i = 0; i < Repeat; ++i ) {
      for ( unsigned j = 0; j < Count; ++j ) {
        const unsigned long correct = 0 == j % 3 ? nw : nsum;
        ASSERT_EQ( (ScalarType) correct, -result[i].value[j] );
      }
    }
  }
};

template< typename ScalarType, class DeviceType >
class TestReduceDynamic
{
public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;

  TestReduceDynamic( const size_type nwork )
  {
    run_test_dynamic( nwork );
    run_test_dynamic_minmax( nwork );
    run_test_dynamic_final( nwork );
  }

  void run_test_dynamic( const size_type nwork )
  {
    typedef Test::RuntimeReduceFunctor< ScalarType, execution_space > functor_type;

    enum { Count = 3 };
    enum { Repeat = 100 };

    ScalarType result[ Repeat ][ Count ];

    const unsigned long nw   = nwork;
    const unsigned long nsum = nw % 2 ? nw * ( ( nw + 1 ) / 2 )
                                      : ( nw / 2 ) * ( nw + 1 );

    for ( unsigned i = 0; i < Repeat; ++i ) {
      if ( i % 2 == 0 ) {
        Kokkos::parallel_reduce( nwork, functor_type( nwork, Count ), result[i] );
      }
      else {
        Kokkos::parallel_reduce( "Reduce", nwork, functor_type( nwork, Count ), result[i] );
      }
    }

    for ( unsigned i = 0; i < Repeat; ++i ) {
      for ( unsigned j = 0; j < Count; ++j ) {
        const unsigned long correct = 0 == j % 3 ? nw : nsum;
        ASSERT_EQ( (ScalarType) correct, result[i][j] );
      }
    }
  }

  void run_test_dynamic_minmax( const size_type nwork )
  {
    typedef Test::RuntimeReduceMinMax< ScalarType, execution_space > functor_type;

    enum { Count = 2 };
    enum { Repeat = 100 };

    ScalarType result[ Repeat ][ Count ];

    for ( unsigned i = 0; i < Repeat; ++i ) {
      if ( i % 2 == 0 ) {
        Kokkos::parallel_reduce( nwork, functor_type( nwork, Count ), result[i] );
      }
      else {
        Kokkos::parallel_reduce( "Reduce", nwork, functor_type( nwork, Count ), result[i] );
      }
    }

    for ( unsigned i = 0; i < Repeat; ++i ) {
      for ( unsigned j = 0; j < Count; ++j ) {
        if ( nwork == 0 )
        {
          ScalarType amin( std::numeric_limits< ScalarType >::min() );
          ScalarType amax( std::numeric_limits< ScalarType >::max() );
          const ScalarType correct = ( j % 2 ) ? amax : amin;
          ASSERT_EQ( (ScalarType) correct, result[i][j] );
        }
        else {
          const unsigned long correct = j % 2 ? 1 : nwork;
          ASSERT_EQ( (ScalarType) correct, result[i][j] );
        }
      }
    }
  }

  void run_test_dynamic_final( const size_type nwork )
  {
    typedef Test::RuntimeReduceFunctorFinal< execution_space > functor_type;

    enum { Count = 3 };
    enum { Repeat = 100 };

    typename functor_type::scalar_type result[ Repeat ][ Count ];

    const unsigned long nw   = nwork;
    const unsigned long nsum = nw % 2 ? nw * ( ( nw + 1 ) / 2 )
                                      : ( nw / 2 ) * ( nw + 1 );

    for ( unsigned i = 0; i < Repeat; ++i ) {
      if ( i % 2 == 0 ) {
        Kokkos::parallel_reduce( nwork, functor_type( nwork, Count ), result[i] );
      }
      else {
        Kokkos::parallel_reduce( "TestKernelReduce", nwork, functor_type( nwork, Count ), result[i] );
      }

    }

    for ( unsigned i = 0; i < Repeat; ++i ) {
      for ( unsigned j = 0; j < Count; ++j ) {
        const unsigned long correct = 0 == j % 3 ? nw : nsum;
        ASSERT_EQ( (ScalarType) correct, -result[i][j] );
      }
    }
  }
};

template< typename ScalarType, class DeviceType >
class TestReduceDynamicView
{
public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;

  TestReduceDynamicView( const size_type nwork )
  {
    run_test_dynamic_view( nwork );
  }

  void run_test_dynamic_view( const size_type nwork )
  {
    typedef Test::RuntimeReduceFunctor< ScalarType, execution_space > functor_type;

    typedef Kokkos::View< ScalarType*, DeviceType > result_type;
    typedef typename result_type::HostMirror result_host_type;

    const unsigned CountLimit = 23;

    const unsigned long nw   = nwork;
    const unsigned long nsum = nw % 2 ? nw * ( ( nw + 1 ) / 2 )
                                      : ( nw / 2 ) * ( nw + 1 );

    for ( unsigned count = 0; count < CountLimit; ++count ) {

      result_type result( "result", count );
      result_host_type host_result = Kokkos::create_mirror( result );

      // Test result to host pointer:

      std::string str( "TestKernelReduce" );
      if ( count % 2 == 0 ) {
        Kokkos::parallel_reduce( nw, functor_type( nw, count ), host_result.ptr_on_device() );
      }
      else {
        Kokkos::parallel_reduce( str, nw, functor_type( nw, count ), host_result.ptr_on_device() );
      }

      for ( unsigned j = 0; j < count; ++j ) {
        const unsigned long correct = 0 == j % 3 ? nw : nsum;
        ASSERT_EQ( host_result( j ), (ScalarType) correct );
        host_result( j ) = 0;
      }
    }
  }
};

} // namespace

// Computes y^T*A*x
// ( modified from kokkos-tutorials/GTC2016/Exercises/ThreeLevelPar )

#if ( ! defined( KOKKOS_ENABLE_CUDA ) ) || defined( KOKKOS_ENABLE_CUDA_LAMBDA )

template< typename ScalarType, class DeviceType >
class TestTripleNestedReduce
{
public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;

  TestTripleNestedReduce( const size_type & nrows, const size_type & ncols
                        , const size_type & team_size, const size_type & vector_length )
  {
    run_test( nrows, ncols, team_size, vector_length );
  }

  void run_test( const size_type & nrows, const size_type & ncols
               , const size_type & team_size, const size_type & vector_length )
  {
    //typedef Kokkos::LayoutLeft Layout;
    typedef Kokkos::LayoutRight Layout;

    typedef Kokkos::View< ScalarType*, DeviceType >            ViewVector;
    typedef Kokkos::View< ScalarType**, Layout, DeviceType >   ViewMatrix;

    ViewVector y( "y", nrows );
    ViewVector x( "x", ncols );
    ViewMatrix A( "A", nrows, ncols );

    typedef Kokkos::RangePolicy<DeviceType> range_policy;

    // Initialize y vector.
    Kokkos::parallel_for( range_policy( 0, nrows ), KOKKOS_LAMBDA ( const int i ) { y( i ) = 1; } );

    // Initialize x vector.
    Kokkos::parallel_for( range_policy( 0, ncols ), KOKKOS_LAMBDA ( const int i ) { x( i ) = 1; } );

    typedef Kokkos::TeamPolicy< DeviceType >                        team_policy;
    typedef typename Kokkos::TeamPolicy< DeviceType >::member_type  member_type;

    // Initialize A matrix, note 2D indexing computation.
    Kokkos::parallel_for( team_policy( nrows, Kokkos::AUTO ), KOKKOS_LAMBDA ( const member_type & teamMember ) {
      const int j = teamMember.league_rank();
      Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember, ncols ), [&] ( const int i ) {
        A( j, i ) = 1;
      } );
    } );

    // Three level parallelism kernel to force caching of vector x.
    ScalarType result = 0.0;
    int chunk_size = 128;
    Kokkos::parallel_reduce( team_policy( nrows / chunk_size, team_size, vector_length ),
                             KOKKOS_LAMBDA ( const member_type & teamMember, double & update ) {
      const int row_start = teamMember.league_rank() * chunk_size;
      const int row_end   = row_start + chunk_size;
      Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember, row_start, row_end ), [&] ( const int i ) {
        ScalarType sum_i = 0.0;
        Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( teamMember, ncols ), [&] ( const int j, ScalarType &innerUpdate ) {
          innerUpdate += A( i, j ) * x( j );
        }, sum_i );
        Kokkos::single( Kokkos::PerThread( teamMember ), [&] () {
          update += y( i ) * sum_i;
        } );
      } );
    }, result );

    const ScalarType solution = (ScalarType) nrows * (ScalarType) ncols;
    ASSERT_EQ( solution, result );
  }
};

#else // #if ( ! defined( KOKKOS_ENABLE_CUDA ) ) || defined( KOKKOS_ENABLE_CUDA_LAMBDA )

template< typename ScalarType, class DeviceType >
class TestTripleNestedReduce
{
public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;

  TestTripleNestedReduce( const size_type &, const size_type
                        , const size_type &, const size_type )
  {}
};

#endif

//--------------------------------------------------------------------------

namespace Test {

namespace ReduceCombinatorical {

template< class Scalar, class Space = Kokkos::HostSpace >
struct AddPlus {
public:
  // Required.
  typedef AddPlus reducer_type;
  typedef Scalar value_type;

  typedef Kokkos::View< value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  result_view_type result;

public:
  AddPlus( value_type & result_ ) : result( &result_ ) {}

  // Required.
  KOKKOS_INLINE_FUNCTION
  void join( value_type & dest, const value_type & src ) const {
    dest += src + 1;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & dest, const volatile value_type & src ) const {
    dest += src + 1;
  }

  // Optional.
  KOKKOS_INLINE_FUNCTION
  void init( value_type & val )  const {
    val = value_type();
  }

  result_view_type result_view() const {
    return result;
  }
};

template< int ISTEAM >
struct FunctorScalar;

template<>
struct FunctorScalar< 0 > {
  Kokkos::View< double > result;

  FunctorScalar( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }
};

template<>
struct FunctorScalar< 1 > {
  typedef Kokkos::TeamPolicy<>::member_type team_type;

  Kokkos::View< double > result;

  FunctorScalar( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_type & team, double & update ) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }
};

template< int ISTEAM >
struct FunctorScalarInit;

template<>
struct FunctorScalarInit< 0 > {
  Kokkos::View< double > result;

  FunctorScalarInit( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }

  KOKKOS_INLINE_FUNCTION
  void init( double & update ) const {
    update = 0.0;
  }
};

template<>
struct FunctorScalarInit< 1 > {
  typedef Kokkos::TeamPolicy<>::member_type team_type;

  Kokkos::View< double > result;

  FunctorScalarInit( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_type & team, double & update ) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void init( double & update ) const {
    update = 0.0;
  }
};

template< int ISTEAM >
struct FunctorScalarFinal;

template<>
struct FunctorScalarFinal< 0 > {
  Kokkos::View<double> result;

  FunctorScalarFinal( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }

  KOKKOS_INLINE_FUNCTION
  void final( double & update ) const {
    result() = update;
  }
};

template<>
struct FunctorScalarFinal< 1 > {
  typedef Kokkos::TeamPolicy<>::member_type team_type;

  Kokkos::View< double > result;

  FunctorScalarFinal( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_type & team, double & update ) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void final( double & update ) const {
    result() = update;
  }
};

template< int ISTEAM >
struct FunctorScalarJoin;

template<>
struct FunctorScalarJoin< 0 > {
  Kokkos::View<double> result;

  FunctorScalarJoin( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }
};

template<>
struct FunctorScalarJoin< 1 > {
  typedef Kokkos::TeamPolicy<>::member_type team_type;

  Kokkos::View< double > result;

  FunctorScalarJoin( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_type & team, double & update ) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }
};

template< int ISTEAM >
struct FunctorScalarJoinFinal;

template<>
struct FunctorScalarJoinFinal< 0 > {
  Kokkos::View< double > result;

  FunctorScalarJoinFinal( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void final( double & update ) const {
    result() = update;
  }
};

template<>
struct FunctorScalarJoinFinal< 1 > {
  typedef Kokkos::TeamPolicy<>::member_type team_type;

  Kokkos::View< double > result;

  FunctorScalarJoinFinal( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_type & team, double & update ) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void final( double & update ) const {
    result() = update;
  }
};

template< int ISTEAM >
struct FunctorScalarJoinInit;

template<>
struct FunctorScalarJoinInit< 0 > {
  Kokkos::View< double > result;

  FunctorScalarJoinInit( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void init( double & update ) const {
    update = 0.0;
  }
};

template<>
struct FunctorScalarJoinInit< 1 > {
  typedef Kokkos::TeamPolicy<>::member_type team_type;

  Kokkos::View< double > result;

  FunctorScalarJoinInit( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_type & team, double & update ) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void init( double & update ) const {
    update = 0.0;
  }
};

template< int ISTEAM >
struct FunctorScalarJoinFinalInit;

template<>
struct FunctorScalarJoinFinalInit< 0 > {
  Kokkos::View<double> result;

  FunctorScalarJoinFinalInit( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void final( double & update ) const {
    result() = update;
  }

  KOKKOS_INLINE_FUNCTION
  void init( double & update ) const {
    update = 0.0;
  }
};

template<>
struct FunctorScalarJoinFinalInit< 1 > {
  typedef Kokkos::TeamPolicy<>::member_type team_type;

  Kokkos::View< double > result;

  FunctorScalarJoinFinalInit( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_type & team, double & update ) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void final( double & update ) const {
    result() = update;
  }

  KOKKOS_INLINE_FUNCTION
  void init( double & update ) const {
    update = 0.0;
  }
};

struct Functor1 {
  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }
};

struct Functor2 {
  typedef double value_type[];

  const unsigned value_count;

  Functor2( unsigned n ) : value_count( n ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned & i, double update[] ) const {
    for ( unsigned j = 0; j < value_count; j++ ) {
      update[j] += i;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init( double dst[] ) const
  {
    for ( unsigned i = 0; i < value_count; ++i ) dst[i] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double dst[],
             const volatile double src[] ) const
  {
    for ( unsigned i = 0; i < value_count; ++i ) dst[i] += src[i];
  }
};

} // namespace ReduceCombinatorical

} // namespace Test

namespace Test {

template< class ExecSpace = Kokkos::DefaultExecutionSpace >
struct TestReduceCombinatoricalInstantiation {
  template< class ... Args >
  static void CallParallelReduce( Args... args ) {
    Kokkos::parallel_reduce( args... );
  }

  template< class ... Args >
  static void AddReturnArgument( Args... args ) {
    Kokkos::View< double, Kokkos::HostSpace > result_view( "ResultView" );
    double expected_result = 1000.0 * 999.0 / 2.0;

    double value = 0;
    Kokkos::parallel_reduce( args..., value );
    ASSERT_EQ( expected_result, value );

    result_view() = 0;
    CallParallelReduce( args..., result_view );
    ASSERT_EQ( expected_result, result_view() );

    value = 0;
    CallParallelReduce( args..., Kokkos::View< double, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >( &value ) );
    ASSERT_EQ( expected_result, value );

    result_view() = 0;
    const Kokkos::View< double, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_const_um = result_view;
    CallParallelReduce( args..., result_view_const_um );
    ASSERT_EQ( expected_result, result_view_const_um() );

    value = 0;
    CallParallelReduce( args..., Test::ReduceCombinatorical::AddPlus< double >( value ) );
    if ( ( Kokkos::DefaultExecutionSpace::concurrency() > 1 ) && ( ExecSpace::concurrency() > 1 ) ) {
      ASSERT_TRUE( expected_result < value );
    }
    else if ( ( Kokkos::DefaultExecutionSpace::concurrency() > 1 ) || ( ExecSpace::concurrency() > 1 ) ) {
      ASSERT_TRUE( expected_result <= value );
    }
    else {
      ASSERT_EQ( expected_result, value );
    }

    value = 0;
    Test::ReduceCombinatorical::AddPlus< double > add( value );
    CallParallelReduce( args..., add );
    if ( ( Kokkos::DefaultExecutionSpace::concurrency() > 1 ) && ( ExecSpace::concurrency() > 1 ) ) {
      ASSERT_TRUE( expected_result < value );
    }
    else if ( ( Kokkos::DefaultExecutionSpace::concurrency() > 1 ) || ( ExecSpace::concurrency() > 1 ) ) {
      ASSERT_TRUE( expected_result <= value );
    }
    else {
      ASSERT_EQ( expected_result, value );
    }
  }

  template< class ... Args >
  static void AddLambdaRange( void*, Args... args ) {
    AddReturnArgument( args..., KOKKOS_LAMBDA ( const int & i, double & lsum ) {
      lsum += i;
    });
  }

  template< class ... Args >
  static void AddLambdaTeam( void*, Args... args ) {
    AddReturnArgument( args..., KOKKOS_LAMBDA ( const Kokkos::TeamPolicy<>::member_type & team, double & update ) {
      update += 1.0 / team.team_size() * team.league_rank();
    });
  }

  template< class ... Args >
  static void AddLambdaRange( Kokkos::InvalidType, Args... args ) {}

  template< class ... Args >
  static void AddLambdaTeam( Kokkos::InvalidType, Args... args ) {}

  template< int ISTEAM, class ... Args >
  static void AddFunctor( Args... args ) {
    Kokkos::View< double > result_view( "FunctorView" );
    auto h_r = Kokkos::create_mirror_view( result_view );
    Test::ReduceCombinatorical::FunctorScalar< ISTEAM > functor( result_view );
    double expected_result = 1000.0 * 999.0 / 2.0;

    AddReturnArgument( args..., functor );
    AddReturnArgument( args..., Test::ReduceCombinatorical::FunctorScalar< ISTEAM >( result_view ) );
    AddReturnArgument( args..., Test::ReduceCombinatorical::FunctorScalarInit< ISTEAM >( result_view ) );
    AddReturnArgument( args..., Test::ReduceCombinatorical::FunctorScalarJoin< ISTEAM >( result_view ) );
    AddReturnArgument( args..., Test::ReduceCombinatorical::FunctorScalarJoinInit< ISTEAM >( result_view ) );

    h_r() = 0;
    Kokkos::deep_copy( result_view, h_r );
    CallParallelReduce( args..., Test::ReduceCombinatorical::FunctorScalarFinal< ISTEAM >( result_view ) );
    Kokkos::deep_copy( h_r, result_view );
    ASSERT_EQ( expected_result, h_r() );

    h_r() = 0;
    Kokkos::deep_copy( result_view, h_r );
    CallParallelReduce( args..., Test::ReduceCombinatorical::FunctorScalarJoinFinal< ISTEAM >( result_view ) );
    Kokkos::deep_copy( h_r, result_view );
    ASSERT_EQ( expected_result, h_r() );

    h_r() = 0;
    Kokkos::deep_copy( result_view, h_r );
    CallParallelReduce( args..., Test::ReduceCombinatorical::FunctorScalarJoinFinalInit< ISTEAM >( result_view ) );
    Kokkos::deep_copy( h_r, result_view );
    ASSERT_EQ( expected_result, h_r() );
  }

  template< class ... Args >
  static void AddFunctorLambdaRange( Args... args ) {
    AddFunctor< 0, Args... >( args... );
#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
    AddLambdaRange( typename std::conditional< std::is_same<ExecSpace, Kokkos::DefaultExecutionSpace>::value, void*, Kokkos::InvalidType >::type(), args... );
#endif
  }

  template< class ... Args >
  static void AddFunctorLambdaTeam( Args... args ) {
    AddFunctor< 1, Args... >( args... );
#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
    AddLambdaTeam( typename std::conditional< std::is_same<ExecSpace, Kokkos::DefaultExecutionSpace>::value, void*, Kokkos::InvalidType >::type(), args... );
#endif
  }

  template< class ... Args >
  static void AddPolicy( Args... args ) {
    int N = 1000;
    Kokkos::RangePolicy< ExecSpace > policy( 0, N );

    AddFunctorLambdaRange( args..., 1000 );
    AddFunctorLambdaRange( args..., N );
    AddFunctorLambdaRange( args..., policy );
    AddFunctorLambdaRange( args..., Kokkos::RangePolicy< ExecSpace >( 0, N ) );
    AddFunctorLambdaRange( args..., Kokkos::RangePolicy< ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >( 0, N ) );
    AddFunctorLambdaRange( args..., Kokkos::RangePolicy< ExecSpace, Kokkos::Schedule<Kokkos::Static> >( 0, N ).set_chunk_size( 10 ) );
    AddFunctorLambdaRange( args..., Kokkos::RangePolicy< ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >( 0, N ).set_chunk_size( 10 ) );

    AddFunctorLambdaTeam( args..., Kokkos::TeamPolicy< ExecSpace >( N, Kokkos::AUTO ) );
    AddFunctorLambdaTeam( args..., Kokkos::TeamPolicy< ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >( N, Kokkos::AUTO ) );
    AddFunctorLambdaTeam( args..., Kokkos::TeamPolicy< ExecSpace, Kokkos::Schedule<Kokkos::Static> >( N, Kokkos::AUTO ).set_chunk_size( 10 ) );
    AddFunctorLambdaTeam( args..., Kokkos::TeamPolicy< ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >( N, Kokkos::AUTO ).set_chunk_size( 10 ) );
  }

  static void execute_a() {
    AddPolicy();
  }

  static void execute_b() {
    std::string s( "Std::String" );
    AddPolicy( s.c_str() );
    AddPolicy( "Char Constant" );
  }

  static void execute_c() {
    std::string s( "Std::String" );
    AddPolicy( s );
  }
};

template< class Scalar, class ExecSpace = Kokkos::DefaultExecutionSpace >
struct TestReducers {
  struct SumFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      value += values( i );
    }
  };

  struct ProdFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      value *= values( i );
    }
  };

  struct MinFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      if ( values( i ) < value ) value = values( i );
    }
  };

  struct MaxFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      if ( values( i ) > value ) value = values( i );
    }
  };

  struct MinLocFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, typename Kokkos::Experimental::MinLoc< Scalar, int >::value_type & value ) const {
      if ( values( i ) < value.val ) {
        value.val = values( i );
        value.loc = i;
      }
    }
  };

  struct MaxLocFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, typename Kokkos::Experimental::MaxLoc< Scalar, int >::value_type & value ) const {
      if ( values( i ) > value.val ) {
        value.val = values( i );
        value.loc = i;
      }
    }
  };

  struct MinMaxLocFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, typename Kokkos::Experimental::MinMaxLoc< Scalar, int >::value_type & value ) const {
      if ( values( i ) > value.max_val ) {
        value.max_val = values( i );
        value.max_loc = i;
      }

      if ( values( i ) < value.min_val ) {
        value.min_val = values( i );
        value.min_loc = i;
      }
    }
  };

  struct BAndFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      value = value & values( i );
    }
  };

  struct BOrFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      value = value | values( i );
    }
  };

  struct BXorFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      value = value ^ values( i );
    }
  };

  struct LAndFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      value = value && values( i );
    }
  };

  struct LOrFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      value = value || values( i );
    }
  };

  struct LXorFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      value = value ? ( !values( i ) ) : values( i );
    }
  };

  static void test_sum( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_sum = 0;

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 100 );
      reference_sum += h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    SumFunctor f;
    f.values = values;
    Scalar init = 0;

    {
      Scalar sum_scalar = init;
      Kokkos::Experimental::Sum< Scalar > reducer_scalar( sum_scalar );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );

      ASSERT_EQ( sum_scalar, reference_sum );

      Scalar sum_scalar_view = reducer_scalar.result_view()();
      ASSERT_EQ( sum_scalar_view, reference_sum );
    }

    {
      Scalar sum_scalar_init = init;
      Kokkos::Experimental::Sum< Scalar > reducer_scalar_init( sum_scalar_init, init );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar_init );

      ASSERT_EQ( sum_scalar_init, reference_sum );

      Scalar sum_scalar_init_view = reducer_scalar_init.result_view()();
      ASSERT_EQ( sum_scalar_init_view, reference_sum );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace> sum_view( "View" );
      sum_view() = init;
      Kokkos::Experimental::Sum< Scalar > reducer_view( sum_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );

      Scalar sum_view_scalar = sum_view();
      ASSERT_EQ( sum_view_scalar, reference_sum );

      Scalar sum_view_view = reducer_view.result_view()();
      ASSERT_EQ( sum_view_view, reference_sum );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > sum_view_init( "View" );
      sum_view_init() = init;
      Kokkos::Experimental::Sum< Scalar > reducer_view_init( sum_view_init, init );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view_init );

      Scalar sum_view_init_scalar = sum_view_init();
      ASSERT_EQ( sum_view_init_scalar, reference_sum );

      Scalar sum_view_init_view = reducer_view_init.result_view()();
      ASSERT_EQ( sum_view_init_view, reference_sum );
    }
  }

  static void test_prod( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_prod = 1;

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 4 + 1 );
      reference_prod *= h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    ProdFunctor f;
    f.values = values;
    Scalar init = 1;

    if ( std::is_arithmetic< Scalar >::value )
    {
      Scalar prod_scalar = init;
      Kokkos::Experimental::Prod< Scalar > reducer_scalar( prod_scalar );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );

      ASSERT_EQ( prod_scalar, reference_prod );

      Scalar prod_scalar_view = reducer_scalar.result_view()();
      ASSERT_EQ( prod_scalar_view, reference_prod );
    }

    {
      Scalar prod_scalar_init = init;
      Kokkos::Experimental::Prod< Scalar > reducer_scalar_init( prod_scalar_init, init );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar_init );

      ASSERT_EQ( prod_scalar_init, reference_prod );

      Scalar prod_scalar_init_view = reducer_scalar_init.result_view()();
      ASSERT_EQ( prod_scalar_init_view, reference_prod );
    }

    if ( std::is_arithmetic< Scalar >::value )
    {
      Kokkos::View< Scalar, Kokkos::HostSpace > prod_view( "View" );
      prod_view() = init;
      Kokkos::Experimental::Prod< Scalar > reducer_view( prod_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );

      Scalar prod_view_scalar = prod_view();
      ASSERT_EQ( prod_view_scalar, reference_prod );

      Scalar prod_view_view = reducer_view.result_view()();
      ASSERT_EQ( prod_view_view, reference_prod );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > prod_view_init( "View" );
      prod_view_init() = init;
      Kokkos::Experimental::Prod< Scalar > reducer_view_init( prod_view_init, init );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view_init );

      Scalar prod_view_init_scalar = prod_view_init();
      ASSERT_EQ( prod_view_init_scalar, reference_prod );

      Scalar prod_view_init_view = reducer_view_init.result_view()();
      ASSERT_EQ( prod_view_init_view, reference_prod );
    }
  }

  static void test_min( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_min = std::numeric_limits< Scalar >::max();

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 100000 );

      if ( h_values( i ) < reference_min ) reference_min = h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    MinFunctor f;
    f.values = values;
    Scalar init = std::numeric_limits< Scalar >::max();

    {
      Scalar min_scalar = init;
      Kokkos::Experimental::Min< Scalar > reducer_scalar( min_scalar );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );

      ASSERT_EQ( min_scalar, reference_min );

      Scalar min_scalar_view = reducer_scalar.result_view()();
      ASSERT_EQ( min_scalar_view, reference_min );
    }

    {
      Scalar min_scalar_init = init;
      Kokkos::Experimental::Min< Scalar > reducer_scalar_init( min_scalar_init, init );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar_init );

      ASSERT_EQ( min_scalar_init, reference_min );

      Scalar min_scalar_init_view = reducer_scalar_init.result_view()();
      ASSERT_EQ( min_scalar_init_view, reference_min );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > min_view( "View" );
      min_view() = init;
      Kokkos::Experimental::Min< Scalar > reducer_view( min_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );

      Scalar min_view_scalar = min_view();
      ASSERT_EQ( min_view_scalar, reference_min );

      Scalar min_view_view = reducer_view.result_view()();
      ASSERT_EQ( min_view_view, reference_min );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > min_view_init( "View" );
      min_view_init() = init;
      Kokkos::Experimental::Min< Scalar > reducer_view_init( min_view_init, init );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view_init );

      Scalar min_view_init_scalar = min_view_init();
      ASSERT_EQ( min_view_init_scalar, reference_min );

      Scalar min_view_init_view = reducer_view_init.result_view()();
      ASSERT_EQ( min_view_init_view, reference_min );
    }
  }

  static void test_max( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_max = std::numeric_limits< Scalar >::min();

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 100000 + 1 );

      if ( h_values( i ) > reference_max ) reference_max = h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    MaxFunctor f;
    f.values = values;
    Scalar init = std::numeric_limits< Scalar >::min();

    {
      Scalar max_scalar = init;
      Kokkos::Experimental::Max< Scalar > reducer_scalar( max_scalar );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );

      ASSERT_EQ( max_scalar, reference_max );

      Scalar max_scalar_view = reducer_scalar.result_view()();
      ASSERT_EQ( max_scalar_view, reference_max );
    }

    {
      Scalar max_scalar_init = init;
      Kokkos::Experimental::Max< Scalar > reducer_scalar_init( max_scalar_init, init );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar_init );

      ASSERT_EQ( max_scalar_init, reference_max );

      Scalar max_scalar_init_view = reducer_scalar_init.result_view()();
      ASSERT_EQ( max_scalar_init_view, reference_max );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > max_view( "View" );
      max_view() = init;
      Kokkos::Experimental::Max< Scalar > reducer_view( max_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );

      Scalar max_view_scalar = max_view();
      ASSERT_EQ( max_view_scalar, reference_max );

      Scalar max_view_view = reducer_view.result_view()();
      ASSERT_EQ( max_view_view, reference_max );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > max_view_init( "View" );
      max_view_init() = init;
      Kokkos::Experimental::Max< Scalar > reducer_view_init( max_view_init, init );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view_init );

      Scalar max_view_init_scalar = max_view_init();
      ASSERT_EQ( max_view_init_scalar, reference_max );

      Scalar max_view_init_view = reducer_view_init.result_view()();
      ASSERT_EQ( max_view_init_view, reference_max );
    }
  }

  static void test_minloc( int N ) {
    typedef typename Kokkos::Experimental::MinLoc< Scalar, int >::value_type value_type;

    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_min = std::numeric_limits< Scalar >::max();
    int reference_loc = -1;

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 100000 );

      if ( h_values( i ) < reference_min ) {
        reference_min = h_values( i );
        reference_loc = i;
      }
      else if ( h_values( i ) == reference_min ) {
        // Make min unique.
        h_values( i ) += std::numeric_limits< Scalar >::epsilon();
      }
    }
    Kokkos::deep_copy( values, h_values );

    MinLocFunctor f;
    f.values = values;
    Scalar init = std::numeric_limits< Scalar >::max();

    {
      value_type min_scalar;
      Kokkos::Experimental::MinLoc< Scalar, int > reducer_scalar( min_scalar );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );

      ASSERT_EQ( min_scalar.val, reference_min );
      ASSERT_EQ( min_scalar.loc, reference_loc );

      value_type min_scalar_view = reducer_scalar.result_view()();
      ASSERT_EQ( min_scalar_view.val, reference_min );
      ASSERT_EQ( min_scalar_view.loc, reference_loc );
    }

    {
      value_type min_scalar_init;
      Kokkos::Experimental::MinLoc< Scalar, int > reducer_scalar_init( min_scalar_init, init );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar_init );

      ASSERT_EQ( min_scalar_init.val, reference_min );
      ASSERT_EQ( min_scalar_init.loc, reference_loc );

      value_type min_scalar_init_view = reducer_scalar_init.result_view()();
      ASSERT_EQ( min_scalar_init_view.val, reference_min );
      ASSERT_EQ( min_scalar_init_view.loc, reference_loc );
    }

    {
      Kokkos::View< value_type, Kokkos::HostSpace > min_view( "View" );
      Kokkos::Experimental::MinLoc< Scalar, int > reducer_view( min_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );

      value_type min_view_scalar = min_view();
      ASSERT_EQ( min_view_scalar.val, reference_min );
      ASSERT_EQ( min_view_scalar.loc, reference_loc );

      value_type min_view_view = reducer_view.result_view()();
      ASSERT_EQ( min_view_view.val, reference_min );
      ASSERT_EQ( min_view_view.loc, reference_loc );
    }

    {
      Kokkos::View< value_type, Kokkos::HostSpace > min_view_init( "View" );
      Kokkos::Experimental::MinLoc< Scalar, int > reducer_view_init( min_view_init, init );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view_init );

      value_type min_view_init_scalar = min_view_init();
      ASSERT_EQ( min_view_init_scalar.val, reference_min );
      ASSERT_EQ( min_view_init_scalar.loc, reference_loc );

      value_type min_view_init_view = reducer_view_init.result_view()();
      ASSERT_EQ( min_view_init_view.val, reference_min );
      ASSERT_EQ( min_view_init_view.loc, reference_loc );
    }
  }

  static void test_maxloc( int N ) {
    typedef typename Kokkos::Experimental::MaxLoc< Scalar, int >::value_type value_type;

    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_max = std::numeric_limits< Scalar >::min();
    int reference_loc = -1;

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 100000 );

      if ( h_values( i ) > reference_max ) {
        reference_max = h_values( i );
        reference_loc = i;
      }
      else if ( h_values( i ) == reference_max ) {
        // Make max unique.
        h_values( i ) -= std::numeric_limits< Scalar >::epsilon();
      }
    }
    Kokkos::deep_copy( values, h_values );

    MaxLocFunctor f;
    f.values = values;
    Scalar init = std::numeric_limits< Scalar >::min();

    {
      value_type max_scalar;
      Kokkos::Experimental::MaxLoc< Scalar, int > reducer_scalar( max_scalar );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );

      ASSERT_EQ( max_scalar.val, reference_max );
      ASSERT_EQ( max_scalar.loc, reference_loc );

      value_type max_scalar_view = reducer_scalar.result_view()();
      ASSERT_EQ( max_scalar_view.val, reference_max );
      ASSERT_EQ( max_scalar_view.loc, reference_loc );
    }

    {
      value_type max_scalar_init;
      Kokkos::Experimental::MaxLoc< Scalar, int > reducer_scalar_init( max_scalar_init, init );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar_init );

      ASSERT_EQ( max_scalar_init.val, reference_max );
      ASSERT_EQ( max_scalar_init.loc, reference_loc );

      value_type max_scalar_init_view = reducer_scalar_init.result_view()();
      ASSERT_EQ( max_scalar_init_view.val, reference_max );
      ASSERT_EQ( max_scalar_init_view.loc, reference_loc );
    }

    {
      Kokkos::View< value_type, Kokkos::HostSpace > max_view( "View" );
      Kokkos::Experimental::MaxLoc< Scalar, int > reducer_view( max_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );

      value_type max_view_scalar = max_view();
      ASSERT_EQ( max_view_scalar.val, reference_max );
      ASSERT_EQ( max_view_scalar.loc, reference_loc );

      value_type max_view_view = reducer_view.result_view()();
      ASSERT_EQ( max_view_view.val, reference_max );
      ASSERT_EQ( max_view_view.loc, reference_loc );
    }

    {
      Kokkos::View< value_type, Kokkos::HostSpace > max_view_init( "View" );
      Kokkos::Experimental::MaxLoc< Scalar, int > reducer_view_init( max_view_init, init );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view_init );

      value_type max_view_init_scalar = max_view_init();
      ASSERT_EQ( max_view_init_scalar.val, reference_max );
      ASSERT_EQ( max_view_init_scalar.loc, reference_loc );

      value_type max_view_init_view = reducer_view_init.result_view()();
      ASSERT_EQ( max_view_init_view.val, reference_max );
      ASSERT_EQ( max_view_init_view.loc, reference_loc );
    }
  }

  static void test_minmaxloc( int N ) {
     typedef typename Kokkos::Experimental::MinMaxLoc< Scalar, int >::value_type value_type;

     Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
     auto h_values = Kokkos::create_mirror_view( values );
     Scalar reference_max = std::numeric_limits< Scalar >::min();
     Scalar reference_min = std::numeric_limits< Scalar >::max();
     int reference_minloc = -1;
     int reference_maxloc = -1;

     for ( int i = 0; i < N; i++ ) {
       h_values( i ) = (Scalar) ( rand() % 100000 );
     }

     for ( int i = 0; i < N; i++ ) {
       if ( h_values( i ) > reference_max ) {
         reference_max = h_values( i );
         reference_maxloc = i;
       }
       else if ( h_values( i ) == reference_max ) {
         // Make max unique.
         h_values( i ) -= std::numeric_limits< Scalar >::epsilon();
       }
     }

     for ( int i = 0; i < N; i++ ) {
       if ( h_values( i ) < reference_min ) {
         reference_min = h_values( i );
         reference_minloc = i;
       }
       else if ( h_values( i ) == reference_min ) {
         // Make min unique.
         h_values( i ) += std::numeric_limits< Scalar >::epsilon();
       }
     }

     Kokkos::deep_copy( values, h_values );

     MinMaxLocFunctor f;
     f.values = values;
     Scalar init_min = std::numeric_limits< Scalar >::max();
     Scalar init_max = std::numeric_limits< Scalar >::min();

     {
       value_type minmax_scalar;
       Kokkos::Experimental::MinMaxLoc< Scalar, int > reducer_scalar( minmax_scalar );
       Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );

       ASSERT_EQ( minmax_scalar.min_val, reference_min );

       for ( int i = 0; i < N; i++ ) {
         if ( ( i == minmax_scalar.min_loc ) && ( h_values( i ) == reference_min ) ) {
           reference_minloc = i;
         }
       }

       ASSERT_EQ( minmax_scalar.min_loc, reference_minloc );
       ASSERT_EQ( minmax_scalar.max_val, reference_max );

       for ( int i = 0; i < N; i++ ) {
         if ( ( i == minmax_scalar.max_loc ) && ( h_values( i ) == reference_max ) ) {
           reference_maxloc = i;
         }
       }

       ASSERT_EQ( minmax_scalar.max_loc, reference_maxloc );

       value_type minmax_scalar_view = reducer_scalar.result_view()();
       ASSERT_EQ( minmax_scalar_view.min_val, reference_min );
       ASSERT_EQ( minmax_scalar_view.min_loc, reference_minloc );
       ASSERT_EQ( minmax_scalar_view.max_val, reference_max );
       ASSERT_EQ( minmax_scalar_view.max_loc, reference_maxloc );
     }

     {
       value_type minmax_scalar_init;
       Kokkos::Experimental::MinMaxLoc< Scalar, int > reducer_scalar_init( minmax_scalar_init, init_min, init_max );
       Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar_init );

       ASSERT_EQ( minmax_scalar_init.min_val, reference_min );
       ASSERT_EQ( minmax_scalar_init.min_loc, reference_minloc );
       ASSERT_EQ( minmax_scalar_init.max_val, reference_max );
       ASSERT_EQ( minmax_scalar_init.max_loc, reference_maxloc );

       value_type minmax_scalar_init_view = reducer_scalar_init.result_view()();
       ASSERT_EQ( minmax_scalar_init_view.min_val, reference_min );
       ASSERT_EQ( minmax_scalar_init_view.min_loc, reference_minloc );
       ASSERT_EQ( minmax_scalar_init_view.max_val, reference_max );
       ASSERT_EQ( minmax_scalar_init_view.max_loc, reference_maxloc );
     }

     {
       Kokkos::View< value_type, Kokkos::HostSpace > minmax_view( "View" );
       Kokkos::Experimental::MinMaxLoc< Scalar, int > reducer_view( minmax_view );
       Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );

       value_type minmax_view_scalar = minmax_view();
       ASSERT_EQ( minmax_view_scalar.min_val, reference_min );
       ASSERT_EQ( minmax_view_scalar.min_loc, reference_minloc );
       ASSERT_EQ( minmax_view_scalar.max_val, reference_max );
       ASSERT_EQ( minmax_view_scalar.max_loc, reference_maxloc );

       value_type minmax_view_view = reducer_view.result_view()();
       ASSERT_EQ( minmax_view_view.min_val, reference_min );
       ASSERT_EQ( minmax_view_view.min_loc, reference_minloc );
       ASSERT_EQ( minmax_view_view.max_val, reference_max );
       ASSERT_EQ( minmax_view_view.max_loc, reference_maxloc );
     }

     {
       Kokkos::View< value_type, Kokkos::HostSpace > minmax_view_init( "View" );
       Kokkos::Experimental::MinMaxLoc< Scalar, int > reducer_view_init( minmax_view_init, init_min, init_max );
       Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view_init );

       value_type minmax_view_init_scalar = minmax_view_init();
       ASSERT_EQ( minmax_view_init_scalar.min_val, reference_min );
       ASSERT_EQ( minmax_view_init_scalar.min_loc, reference_minloc );
       ASSERT_EQ( minmax_view_init_scalar.max_val, reference_max );
       ASSERT_EQ( minmax_view_init_scalar.max_loc, reference_maxloc );

       value_type minmax_view_init_view = reducer_view_init.result_view()();
       ASSERT_EQ( minmax_view_init_view.min_val, reference_min );
       ASSERT_EQ( minmax_view_init_view.min_loc, reference_minloc );
       ASSERT_EQ( minmax_view_init_view.max_val, reference_max );
       ASSERT_EQ( minmax_view_init_view.max_loc, reference_maxloc );
     }
   }

  static void test_BAnd( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_band = Scalar() | ( ~Scalar() );

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 100000 + 1 );
      reference_band = reference_band & h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    BAndFunctor f;
    f.values = values;
    Scalar init = Scalar() | ( ~Scalar() );

    {
      Scalar band_scalar = init;
      Kokkos::Experimental::BAnd< Scalar > reducer_scalar( band_scalar );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );

      ASSERT_EQ( band_scalar, reference_band );
      Scalar band_scalar_view = reducer_scalar.result_view()();

      ASSERT_EQ( band_scalar_view, reference_band );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > band_view( "View" );
      band_view() = init;
      Kokkos::Experimental::BAnd< Scalar > reducer_view( band_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );

      Scalar band_view_scalar = band_view();
      ASSERT_EQ( band_view_scalar, reference_band );

      Scalar band_view_view = reducer_view.result_view()();
      ASSERT_EQ( band_view_view, reference_band );
    }
  }

  static void test_BOr( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_bor = Scalar() & ( ~Scalar() );

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( ( rand() % 100000 + 1 ) * 2 );
      reference_bor = reference_bor | h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    BOrFunctor f;
    f.values = values;
    Scalar init = Scalar() & ( ~Scalar() );

    {
      Scalar bor_scalar = init;
      Kokkos::Experimental::BOr< Scalar > reducer_scalar( bor_scalar );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );

      ASSERT_EQ( bor_scalar, reference_bor );

      Scalar bor_scalar_view = reducer_scalar.result_view()();
      ASSERT_EQ( bor_scalar_view, reference_bor );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > bor_view( "View" );
      bor_view() = init;
      Kokkos::Experimental::BOr< Scalar > reducer_view( bor_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );

      Scalar bor_view_scalar = bor_view();
      ASSERT_EQ( bor_view_scalar, reference_bor );

      Scalar bor_view_view = reducer_view.result_view()();
      ASSERT_EQ( bor_view_view, reference_bor );
    }
  }

  static void test_BXor( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_bxor = Scalar() & ( ~Scalar() );

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( ( rand() % 100000 + 1 ) * 2 );
      reference_bxor = reference_bxor ^ h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    BXorFunctor f;
    f.values = values;
    Scalar init = Scalar() & ( ~Scalar() );

    {
      Scalar bxor_scalar = init;
      Kokkos::Experimental::BXor< Scalar > reducer_scalar( bxor_scalar );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );

      ASSERT_EQ( bxor_scalar, reference_bxor );

      Scalar bxor_scalar_view = reducer_scalar.result_view()();
      ASSERT_EQ( bxor_scalar_view, reference_bxor );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > bxor_view( "View" );
      bxor_view() = init;
      Kokkos::Experimental::BXor< Scalar > reducer_view( bxor_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );

      Scalar bxor_view_scalar = bxor_view();
      ASSERT_EQ( bxor_view_scalar, reference_bxor );

      Scalar bxor_view_view = reducer_view.result_view()();
      ASSERT_EQ( bxor_view_view, reference_bxor );
    }
  }

  static void test_LAnd( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_land = 1;

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 2 );
      reference_land = reference_land && h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    LAndFunctor f;
    f.values = values;
    Scalar init = 1;

    {
      Scalar land_scalar = init;
      Kokkos::Experimental::LAnd< Scalar > reducer_scalar( land_scalar );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );

      ASSERT_EQ( land_scalar, reference_land );

      Scalar land_scalar_view = reducer_scalar.result_view()();
      ASSERT_EQ( land_scalar_view, reference_land );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > land_view( "View" );
      land_view() = init;
      Kokkos::Experimental::LAnd< Scalar > reducer_view( land_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );

      Scalar land_view_scalar = land_view();
      ASSERT_EQ( land_view_scalar, reference_land );

      Scalar land_view_view = reducer_view.result_view()();
      ASSERT_EQ( land_view_view, reference_land );
    }
  }

  static void test_LOr( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_lor = 0;

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 2 );
      reference_lor = reference_lor || h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    LOrFunctor f;
    f.values = values;
    Scalar init = 0;

    {
      Scalar lor_scalar = init;
      Kokkos::Experimental::LOr< Scalar > reducer_scalar( lor_scalar );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );

      ASSERT_EQ( lor_scalar, reference_lor );

      Scalar lor_scalar_view = reducer_scalar.result_view()();
      ASSERT_EQ( lor_scalar_view, reference_lor );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > lor_view( "View" );
      lor_view() = init;
      Kokkos::Experimental::LOr< Scalar > reducer_view( lor_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );

      Scalar lor_view_scalar = lor_view();
      ASSERT_EQ( lor_view_scalar, reference_lor );

      Scalar lor_view_view = reducer_view.result_view()();
      ASSERT_EQ( lor_view_view, reference_lor );
    }
  }

  static void test_LXor( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_lxor = 0;

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 2 );
      reference_lxor = reference_lxor ? ( !h_values( i ) ) : h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    LXorFunctor f;
    f.values = values;
    Scalar init = 0;

    {
      Scalar lxor_scalar = init;
      Kokkos::Experimental::LXor< Scalar > reducer_scalar( lxor_scalar );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );

      ASSERT_EQ( lxor_scalar, reference_lxor );

      Scalar lxor_scalar_view = reducer_scalar.result_view()();
      ASSERT_EQ( lxor_scalar_view, reference_lxor );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > lxor_view( "View" );
      lxor_view() = init;
      Kokkos::Experimental::LXor< Scalar > reducer_view( lxor_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );

      Scalar lxor_view_scalar = lxor_view();
      ASSERT_EQ( lxor_view_scalar, reference_lxor );

      Scalar lxor_view_view = reducer_view.result_view()();
      ASSERT_EQ( lxor_view_view, reference_lxor );
    }
  }

  static void execute_float() {
    test_sum( 10001 );
    test_prod( 35 );
    test_min( 10003 );
    test_minloc( 10003 );
    test_max( 10007 );
    test_maxloc( 10007 );
    test_minmaxloc( 10007 );
  }

  static void execute_integer() {
    test_sum( 10001 );
    test_prod( 35 );
    test_min( 10003 );
    test_minloc( 10003 );
    test_max( 10007 );
    test_maxloc( 10007 );
    test_minmaxloc( 10007 );
    test_BAnd( 35 );
    test_BOr( 35 );
    test_BXor( 35 );
    test_LAnd( 35 );
    test_LOr( 35 );
    test_LXor( 35 );
  }

  static void execute_basic() {
    test_sum( 10001 );
    test_prod( 35 );
  }
};

} // namespace Test
