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

#ifndef USESCASES_LINALG_BLAS_HPP
#define USESCASES_LINALG_BLAS_HPP

#include <cmath>
#include <utility>
#include <ParallelComm.hpp>
#include <KokkosArray_View.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class Scalar , class Layout , class DeviceType > struct Dot ;

template< class Scalar , class Layout , class DeviceType > struct Dot1 ;

template< typename ScalarA ,
          typename ScalarY ,
          class Layout , class Device >
struct Scale ;

template< typename ScalarA ,
          typename ScalarX ,
          typename ScalarY ,
          class Layout , class Device >
struct AXPY ;

template< typename ScalarA ,
          typename ScalarX ,
          typename ScalarB ,
          typename ScalarY ,
          typename ScalarW ,
          class Layout , class Device >
struct WAXPBY ;

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

//----------------------------------------------------------------------------

template< typename ScalarX /* Allow mix of const and non-const */ ,
          typename ScalarY /* Allow mix of const and non-const */ ,
          class Layout ,
          class Device ,
          class MX /* Allow any management type */ ,
          class MY /* Allow any management type */ >
double dot( const size_t n ,
            const View< ScalarX * , Layout , Device , MX > & x ,
            const View< ScalarY * , Layout , Device , MY > & y ,
            comm::Machine machine )
{
  double global_result = 0 ;

  Impl::Dot< ScalarX , Layout , Device >( n , x , y , global_result );

#if defined( HAVE_MPI )

  double local_result = global_result ;

  MPI_Allreduce( & local_result , & global_result , 1 ,
                 MPI_DOUBLE , MPI_SUM , machine.mpi_comm );

#endif

  return global_result ;
}

//----------------------------------------------------------------------------

template< typename ScalarX /* Allow mix of const and non-const */ ,
          class Layout ,
          class Device ,
          class MX /* Allow any management type */ >
double norm2( const size_t n ,
              const View< ScalarX * , Layout , Device , MX > & x ,
              comm::Machine machine )
{
  double global_result = 0 ;

  Impl::Dot1< ScalarX , Layout , Device >( n , x , global_result );

#if defined( HAVE_MPI )

  double local_result = global_result ;

  MPI_Allreduce( & local_result , & global_result , 1 ,
                 MPI_DOUBLE , MPI_SUM , machine.mpi_comm );

#endif

  return std::sqrt( global_result );
}

//----------------------------------------------------------------------------

template< typename ScalarA ,
          typename ScalarX ,
          class Layout ,
          class Device ,
          class MX >
void scale( const size_t n ,
            const ScalarA & alpha ,
            const View< ScalarX , Layout , Device , MX > & x )
{
  Impl::Scale< ScalarA , ScalarX , Layout , Device >( n , alpha , x );
}

//----------------------------------------------------------------------------

template< typename ScalarA ,
          typename ScalarX ,
          typename ScalarY ,
          class Layout ,
          class Device ,
          class MX ,
          class MY >
void axpy( const size_t n ,
           const ScalarA & alpha ,
           const View< ScalarX , Layout , Device , MX > & x ,
           const View< ScalarY , Layout , Device , MY > & y )
{
  Impl::AXPY< ScalarA, ScalarX, ScalarY , Layout , Device >
    ( n , alpha , x , y );
}

//----------------------------------------------------------------------------
// w = alpha * x + beta * y

template< typename ScalarA ,
          typename ScalarX ,
          typename ScalarB ,
          typename ScalarY ,
          typename ScalarW ,
          class Layout , class Device ,
          class MX , class MY , class MW >
void waxpby( const size_t n ,
             const ScalarA & alpha ,
             const View< ScalarX * , Layout , Device , MX > & x ,
             const ScalarB & beta ,
             const View< ScalarY * , Layout , Device , MY > & y ,
             const View< ScalarW * , Layout , Device , MW > & w )
{
  Impl::WAXPBY<ScalarA,ScalarX,ScalarB,ScalarY,ScalarW,Layout,Device>
    ( n , alpha , x , beta , y , w );
}

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< typename Scalar , class Layout , class DeviceType >
struct Dot
{
private:

  typedef View< const Scalar*, Layout, DeviceType , MemoryUnmanaged >  vector_const_type ;

  const vector_const_type x ;
  const vector_const_type y ;

public:

  typedef DeviceType  device_type ; // Manycore device
  typedef double      value_type ;  // Reduction value

  template< typename ScalarX , typename ScalarY , class MX , class MY >
  inline
  Dot( const size_t n ,
       const View< ScalarX * , Layout , DeviceType , MX > & arg_x ,
       const View< ScalarY * , Layout , DeviceType , MY > & arg_y ,
       double & result )
    : x( arg_x ), y( arg_y )
  {
    result = parallel_reduce( n , *this );
  }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType & i , value_type & update ) const 
  { update += x(i) * y(i); }
    
  KOKKOSARRAY_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
  { update += source;    }
    
  KOKKOSARRAY_INLINE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }
}; // Dot

//----------------------------------------------------------------------------

template< typename Scalar , class Layout , class DeviceType >
struct Dot1
{
private:

  typedef View< const Scalar*, Layout, DeviceType , MemoryUnmanaged >  vector_const_type ;

  const vector_const_type x ;

public:

  typedef DeviceType  device_type ; // Manycore device
  typedef double      value_type ;  // Reduction value

  template< typename ScalarX , class MX >
  inline
  Dot1( const size_t n ,
        const View< ScalarX * , Layout , DeviceType , MX > & arg_x ,
        double & result )
    : x( arg_x )
  {
    result = parallel_reduce( n , *this );
  }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType & i , value_type & update ) const 
  { update += x(i) * x(i) ; }
    
  KOKKOSARRAY_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
  { update += source ; }
    
  KOKKOSARRAY_INLINE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }
}; // Dot

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template < typename ScalarA ,
           typename ScalarX ,
           typename ScalarB ,
           typename ScalarY ,
           typename ScalarW ,
           class Layout , class DeviceType >
struct WAXPBY
{
private:

  const View<       ScalarW , Layout , DeviceType , MemoryUnmanaged >  w ;
  const View< const ScalarX , Layout , DeviceType , MemoryUnmanaged >  x ;
  const View< const ScalarY , Layout , DeviceType , MemoryUnmanaged >  y ;
  const ScalarA  alpha ;
  const ScalarB  beta ;

public:

  typedef DeviceType  device_type ;

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType inode ) const
  {
    w(inode) = alpha * x(inode) + beta * y(inode);
  }

  template< class MX , class MY , class MW >
  inline
  WAXPBY( const size_t  n ,
          const ScalarA & arg_alpha ,
          const View< ScalarX , Layout , DeviceType , MX > & arg_x ,
          const ScalarB & arg_beta ,
          const View< ScalarY , Layout , DeviceType , MY > & arg_y ,
          const View< ScalarW , Layout , DeviceType , MW > & arg_w )
    : w( arg_w ), x( arg_x ), y( arg_y )
    , alpha( arg_alpha ), beta( arg_beta )
  {
    parallel_for( n , *this );
  }
}; // WAXPBY

//----------------------------------------------------------------------------

template < typename ScalarB ,
           typename ScalarW ,
           class Layout , class DeviceType >
struct Scale
{
private:

  const View< ScalarW , Layout , DeviceType , MemoryUnmanaged >  w ;
  const ScalarB  beta ;

public:

  typedef DeviceType  device_type ;

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType & i ) const
  { w(i) *= beta ; }

  template< class MW >
  inline
  Scale( const size_t  n ,
         const ScalarB & arg_beta ,
         const View< ScalarW , Layout , DeviceType , MW > & arg_w )
    : w( arg_w )
    , beta( arg_beta )
  {
    parallel_for( n , *this );
  }
};

//----------------------------------------------------------------------------

template < typename ScalarA ,
           typename ScalarX ,
           typename ScalarW ,
           class Layout , class DeviceType >
struct AXPY
{
private:

  const View<       ScalarW , Layout , DeviceType , MemoryUnmanaged >  w ;
  const View< const ScalarX , Layout , DeviceType , MemoryUnmanaged >  x ;
  const ScalarA  alpha ;

public:

  typedef DeviceType  device_type ;

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType & i ) const
  { w(i) += alpha * x(i); }

  template< class MX , class MW >
  inline
  AXPY( const size_t  n ,
        const ScalarA & arg_alpha ,
        const View< ScalarX , Layout , DeviceType , MX > & arg_x ,
        const View< ScalarW , Layout , DeviceType , MW > & arg_w )
    : w( arg_w ), x( arg_x )
    , alpha( arg_alpha )
  {
    parallel_for( n , *this );
  }
}; // WAXPBY

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef USESCASES_LINALG_BLAS_HPP */


