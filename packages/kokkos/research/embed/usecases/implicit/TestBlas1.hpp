/*
//@HEADER
// ************************************************************************
// 
//    Kokkos: Manycore Performance-Portable Multidimensional Arrays
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
// Questions? Contact H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_BLAS1_HPP
#define KOKKOS_BLAS1_HPP

#include <cmath>
#include <limits>
#include <ParallelDataMap.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {

template < typename TypeScalar , typename TypeVector , class Device >
struct WAXPBY
{
private:

  const TypeScalar  alpha ;
  const TypeScalar  beta ;
  const TypeVector * const x ;
  const TypeVector * const y ;
        TypeVector * const w ;

public:

  typedef Device  execution_space ;

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType inode ) const
  { w[inode] = alpha * x[inode] + beta * y[inode]; }

  inline
  WAXPBY( const unsigned count ,
          const TypeScalar & arg_alpha ,
          const TypeVector * const arg_x ,
          const TypeScalar & arg_beta ,
          const TypeVector * const arg_y ,
                TypeVector * const arg_w )
    : alpha( arg_alpha )
    , beta(  arg_beta )
    , x( arg_x )
    , y( arg_y )
    , w( arg_w )
  {
    parallel_for( count , *this );
  }
}; // WAXPBY


template < typename TypeScalar , typename TypeVector , class Device >
struct AXPY
{
private:

  const TypeScalar  alpha ;
  const TypeVector * const x ;
        TypeVector * const y ;

public:

  typedef Device  execution_space ;

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType inode ) const
  { y[inode] = alpha * x[inode] + y[inode]; }

  inline
  AXPY( const unsigned count ,
        const TypeScalar & arg_alpha ,
        const TypeVector * const arg_x ,
              TypeVector * const arg_y )
    : alpha( arg_alpha )
    , x( arg_x )
    , y( arg_y )
  {
    parallel_for( count , *this );
  }
}; // AXPY


template < typename TypeScalar , typename TypeVector , class Device >
struct XPBY
{
private:

  const TypeScalar  beta ;
  const TypeVector * const x ;
        TypeVector * const y ;

public:

  typedef Device  execution_space ;

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType inode ) const
  { y[inode] = x[inode] + beta * y[inode]; }

  inline
  XPBY( const unsigned count ,
        const TypeVector * const arg_x ,
        const TypeScalar & arg_beta ,
              TypeVector * const arg_y )
    : beta(  arg_beta )
    , x( arg_x )
    , y( arg_y )
  {
    parallel_for( count , *this );
  }
}; // XPBY

//----------------------------------------------------------------------------

template< typename TypeX , typename TypeY , class Device >
struct Dot
{
  typedef Device execution_space ;
  typedef double value_type ;

  const TypeX * const x ;
  const TypeY * const y ;

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
  { update += source;    }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType & i , value_type & update ) const
  { update += x[i] * y[i] ; }

  Dot( const unsigned count ,
       const TypeX * const arg_x ,
       const TypeY * const arg_y ,
       double & result )
  : x( arg_x ), y( arg_y )
  {
    parallel_reduce( count , *this , result );
  }
};

template< typename TypeX , class Device >
struct Dot< TypeX , void , Device >
{
  typedef Device execution_space ;
  typedef double value_type ;

  const TypeX * const x ;

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
  { update += source;    }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType & i , value_type & update ) const
  { update += x[i] * x[i] ; }

  Dot( const unsigned count ,
       const TypeX * const arg_x ,
       double & result )
  : x( arg_x )
  {
    parallel_reduce( count , *this , result );
  }
};

}

//----------------------------------------------------------------------------

namespace Kokkos {

template< typename TypeScalar ,
          class    TypeVector ,
          class    Device ,
          class    M >
void waxpby( const ParallelDataMap & map ,
             const TypeScalar      & alpha ,
             const View<TypeVector,LayoutRight,Device,M> & x ,
             const TypeScalar      & beta ,
             const View<TypeVector,LayoutRight,Device,M> & y ,
             const View<TypeVector,LayoutRight,Device,M> & w )
{
  typedef View<TypeVector,LayoutRight,Device,M> vector_type ;
  typedef typename vector_type::scalar_type     vector_scalar_type ;

  size_t count = map.count_owned ;

  for ( unsigned r = 1 ; r < x.Rank ; ++r ) {
    count *= x.dimension(r);
  }

  WAXPBY<TypeScalar,vector_scalar_type,Device>( count , alpha , x.ptr_on_device() ,
                                                        beta ,  y.ptr_on_device() ,
                                                                w.ptr_on_device() );
}

template< typename TypeScalar ,
          class    TypeVector ,
          class    Device ,
          class    M >
void axpy( const ParallelDataMap & map ,
           const TypeScalar      & alpha ,
           const View<TypeVector,LayoutRight,Device,M> & x ,
           const View<TypeVector,LayoutRight,Device,M> & y )
{
  typedef View<TypeVector,LayoutRight,Device,M> vector_type ;
  typedef typename vector_type::scalar_type     vector_scalar_type ;

  size_t count = map.count_owned ;

  for ( unsigned r = 1 ; r < x.Rank ; ++r ) {
    count *= x.dimension(r);
  }

  AXPY<TypeScalar,vector_scalar_type,Device>( count , alpha , x.ptr_on_device() , y.ptr_on_device() );
}

template< typename TypeScalar ,
          class    TypeVector ,
          class    Device ,
          class    M >
void xpby( const ParallelDataMap & map ,
           const View<TypeVector,LayoutRight,Device,M> & x ,
           const TypeScalar      & beta ,
           const View<TypeVector,LayoutRight,Device,M> & y )
{
  typedef View<TypeVector,LayoutRight,Device,M> vector_type ;
  typedef typename vector_type::scalar_type     vector_scalar_type ;

  size_t count = map.count_owned ;

  for ( unsigned r = 1 ; r < x.Rank ; ++r ) {
    count *= x.dimension(r);
  }

  XPBY<TypeScalar,vector_scalar_type,Device>( count , x.ptr_on_device() , beta ,  y.ptr_on_device() );
}

template< class TypeVector , class Device , class M >
double dot( const ParallelDataMap & map ,
            const View<TypeVector,LayoutRight,Device,M> & x )
{
  typedef View<TypeVector,LayoutRight,Device,M> vector_type ;
  typedef typename vector_type::scalar_type     vector_scalar_type ;

  size_t count = map.count_owned ;

  for ( unsigned r = 1 ; r < x.Rank ; ++r ) {
    count *= x.dimension(r);
  }

  double local_result , result ;

  Dot<vector_scalar_type,void,Device>( count , x.ptr_on_device() , local_result );

#if defined( KOKKOS_HAVE_MPI )
  MPI_Allreduce( & local_result , & result , 1 , MPI_DOUBLE , MPI_SUM , map.machine.mpi_comm );
#else
  result = local_result ;
#endif

  return result ;
}

template< class TypeVector , class Device , class M >
double dot( const ParallelDataMap & map ,
            const View<TypeVector,LayoutRight,Device,M> & x ,
            const View<TypeVector,LayoutRight,Device,M> & y )
{
  typedef View<TypeVector,LayoutRight,Device,M> vector_type ;
  typedef typename vector_type::scalar_type     vector_scalar_type ;

  size_t count = map.count_owned ;

  for ( unsigned r = 1 ; r < x.Rank ; ++r ) {
    count *= x.dimension(r);
  }

  double local_result , result ;

  Dot<vector_scalar_type,vector_scalar_type,Device>( count , x.ptr_on_device() , y.ptr_on_device() , local_result );

#if defined( KOKKOS_HAVE_MPI )
  MPI_Allreduce( & local_result , & result , 1 , MPI_DOUBLE , MPI_SUM , map.machine.mpi_comm );
#else
  result = local_result ;
#endif

  return result ;
}

}

//----------------------------------------------------------------------------

#endif

