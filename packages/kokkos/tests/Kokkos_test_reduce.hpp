/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

namespace Kokkos {

//----------------------------------------------------------------------------

template< typename Scalar , unsigned Rank , class Device > struct ReduceSum ;

template< typename Scalar , class Device >
struct ReduceSum<Scalar,1,Device> {
  typedef typename Device::device_type            device_type ;
  typedef Kokkos::MDArrayView<Scalar,device_type> reduce_type ;
  typedef typename reduce_type::value_type        reduce_value_type ;
  typedef typename device_type::size_type         size_type ;

  KOKKOS_DEVICE_FUNCTION
  static void join( const reduce_type & result ,
                    const reduce_type & source )
  {
    const size_type n = result.dimension(0);
    for ( size_type i = 0 ; i < n ; ++i ) {
      result(i) += source(i);
    }
  }

  KOKKOS_DEVICE_FUNCTION
  static void init( const reduce_type & result )
  {
    const size_type n = result.dimension(0);
    for ( size_type i = 0 ; i < n ; ++i ) {
      result(i) = 0 ;
    }
  }
};

template< typename Scalar , unsigned Rank , class DeviceMap >
class MDArrayFill ;

template< typename Scalar , class DeviceMap >
class MDArrayFill<Scalar,1,DeviceMap> {
public:
  typedef typename DeviceMap::device_type device_type ;
  typedef typename DeviceMap::size_type   size_type ;
  
  MDArrayView<Scalar,DeviceMap> array ;
  Scalar                        value ;

  MDArrayFill( const MDArrayView<Scalar,DeviceMap> & a ,
               const Scalar & v )
    : array( a ), value( v ) {}

  KOKKOS_DEVICE_AND_HOST_FUNCTION
  size_type work_count() const { return array.dimension(0); }

  KOKKOS_DEVICE_FUNCTION
  void operator()( size_type iwork ) const
  {
    array(iwork) = value ;
  }
};

template< typename Scalar , class DeviceMap >
class MDArrayFill<Scalar,2,DeviceMap> {
public:
  typedef typename DeviceMap::device_type device_type ;
  typedef typename DeviceMap::size_type   size_type ;
  
  MDArrayView<Scalar,DeviceMap> array ;
  Scalar                        value ;
  size_type                     count ;

  MDArrayFill( const MDArrayView<Scalar,DeviceMap> & a ,
               const Scalar & v )
    : array( a ), value( v ), count( a.dimension(0) ) {}

  KOKKOS_DEVICE_AND_HOST_FUNCTION
  size_type work_count() const { return array.dimension(1); }

  KOKKOS_DEVICE_FUNCTION
  void operator()( size_type iwork ) const
  {
    for ( size_type i = 0 ; i < count ; ++i ) {
      array(i,iwork) = value ;
    }
  }
};

// If the array rank were part of the type information then
// utility functions such as 'fill' could have a simpler API.
// The simpler API would shift type specifics to the templated implementation.

template< typename Scalar , class DeviceMap >
void fill( const MDArrayView< Scalar , DeviceMap > & array , const Scalar & value )
{
  switch( array.rank() ) {
  case 1 :
    Kokkos::parallel_for( MDArrayFill<Scalar,1,DeviceMap>( array , value ) );
    break ;
  case 2 :
    Kokkos::parallel_for( MDArrayFill<Scalar,2,DeviceMap>( array , value ) );
    break ;
  default:
    break ;
  }
}

}

//----------------------------------------------------------------------------
/** \brief  Compute momentum and kinetic energy.
 *
 *  Inherit the commonly used summation reduction operator for rank-1 arrays.
 */
template< typename Scalar , class DeviceMap >
struct MomentumAndKineticEnergy
  : public Kokkos::ReduceSum<Scalar,1,DeviceMap>
{
  typedef Kokkos::ReduceSum<Scalar,1,DeviceMap> reduce_operator ;

  typedef          DeviceMap                    device_map_type ;
  typedef typename DeviceMap::size_type         size_type ;
  typedef typename reduce_operator::reduce_type reduce_type ;

  Kokkos::MDArrayView<Scalar,DeviceMap> velocity ; // (dimension,work)
  Kokkos::MDArrayView<Scalar,DeviceMap> mass ;     // (work)
  size_type                             dimension ;

  MomentumAndKineticEnergy(
    const Kokkos::MDArrayView<Scalar,DeviceMap> & arg_velocity ,
    const Kokkos::MDArrayView<Scalar,DeviceMap> & arg_mass )
    : velocity( arg_velocity ), mass( arg_mass ) ,
      dimension( velocity.dimension(0) ) {}

  KOKKOS_DEVICE_AND_HOST_FUNCTION
  size_type work_count() const { return mass.dimension(0); }

  KOKKOS_DEVICE_FUNCTION
  void operator()( size_type iwork , const reduce_type & result ) const
  {
    const Scalar m = mass(iwork);
    Scalar ke = 0 ;

    for ( size_type i = 0 ; i < dimension ; ++i ) {
      const Scalar v = velocity(i,iwork);
      result(i) += m * v ; // Momentum
      ke += v * v ; // Energy
    }
    result(dimension) += 0.5 * m * ke ;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename Scalar , class DeviceMap >
void test_mom_ke()
{
  typedef typename DeviceMap::device_type device_type ;

  const int parallel_work_length = 1000000 ;

  DeviceMap map( parallel_work_length );

  Kokkos::MDArrayView< Scalar , DeviceMap > mass , velocity ;
  Kokkos::MDArrayView< Scalar , device_type > result ;

  // Create arrays mapped onto the device
  mass      = map.template create_mdarray<Scalar>();
  velocity  = map.template create_mdarray<Scalar>( 3 );

  // Create unmapped result array directly on the device
  // { momentum , kinetic energy }
  result    = device_type::template create_mdarray<Scalar>( 4 );

  // Potential alternative API for creating arrays
  //
  // mass.allocate( map );
  // velocity.allocate( 3 , map );
  // result.allocate( 4 ); // Allocated on device with no map

  // Execute the parallel kernels on the arrays:

  Kokkos::fill( mass , (Scalar) 1.0 );
  Kokkos::fill( velocity , (Scalar) 1.0 );

  Kokkos::parallel_reduce( MomentumAndKineticEnergy<Scalar,DeviceMap>( velocity , mass ) , result );

/*
  std::cout << " { " << result(0)
            << " , " << result(1)
            << " , " << result(2)
            << " , " << result(3)
            << " }" << std::endl ;
*/
}


