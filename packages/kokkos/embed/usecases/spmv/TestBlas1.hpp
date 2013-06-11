/*
//@HEADER
// ************************************************************************
// 
//    KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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

#ifndef KOKKOSARRAY_BLAS1_HPP
#define KOKKOSARRAY_BLAS1_HPP

#include <cmath>
#include <limits>
#include <TestCrsMatrix.hpp>

//----------------------------------------------------------------------------

namespace KokkosArray {

template < typename TypeAlpha ,
           typename TypeX ,
           typename TypeBeta ,
           typename TypeY ,
           typename TypeW >
struct WAXPBY ;

template< class TypeX , class TypeY >
struct Dot ;

template< typename Scalar , unsigned N , class PX , class PY >
KOKKOSARRAY_INLINE_FUNCTION
double dot( const Array< Scalar , N , PX > & x ,
            const Array< Scalar , N , PY > & y )
{
  double result = 0 ;
  for ( unsigned i = 0 ; i < N ; ++i ) {
    result += x[i] * y[i] ;
  }
  return result ;
}

template< typename Scalar , unsigned N , class PX >
KOKKOSARRAY_INLINE_FUNCTION
double dot( const Array< Scalar , N , PX > & x )
{
  double result = 0 ;
  for ( int i = 0 ; i < N ; ++i ) {
    result += x[i] * x[i] ;
  }
  return result ;
}

KOKKOSARRAY_INLINE_FUNCTION
double dot( const double & x , const double & y )
{ return x * y ; }

KOKKOSARRAY_INLINE_FUNCTION
double dot( const double & x )
{ return x * x ; }


}

//----------------------------------------------------------------------------

namespace KokkosArray {

template < typename Scalar , unsigned N , class Device >
struct WAXPBY< Scalar ,
               View< Array<Scalar,N> * , LayoutRight , Device > ,
               Scalar ,
               View< Array<Scalar,N> * , LayoutRight , Device > ,
               View< Array<Scalar,N> * , LayoutRight , Device > >
{
private:

  const View< Array<Scalar,N> *, LayoutRight , Device , MemoryUnmanaged >  w ;
  const View< Array<Scalar,N> *, LayoutRight , Device , MemoryUnmanaged >  x ;
  const View< Array<Scalar,N> *, LayoutRight , Device , MemoryUnmanaged >  y ;
  const Scalar  alpha ;
  const Scalar  beta ;

public:

  typedef Device  device_type ;

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType inode ) const
  { w(inode) = alpha * x(inode) + beta * y(inode); }

  inline
  WAXPBY( const size_t  n ,
          const Scalar & arg_alpha ,
          const View< Array<Scalar,N> *, LayoutRight , Device > & arg_x ,
          const Scalar & arg_beta ,
          const View< Array<Scalar,N> *, LayoutRight , Device > & arg_y ,
          const View< Array<Scalar,N> *, LayoutRight , Device > & arg_w )
    : w( arg_w ), x( arg_x ), y( arg_y )
    , alpha( arg_alpha ), beta( arg_beta )
  {
    parallel_for( n , *this );
  }
}; // WAXPBY

template < typename Scalar , unsigned N , class Device >
struct WAXPBY< Scalar ,
               View< Array<Scalar,N> * , LayoutRight , Device > ,
               void , void ,
               View< Array<Scalar,N> * , LayoutRight , Device > >
{
private:

  const View< Array<Scalar,N> *, LayoutRight , Device , MemoryUnmanaged >  w ;
  const View< Array<Scalar,N> *, LayoutRight , Device , MemoryUnmanaged >  x ;
  const Scalar  alpha ;

public:

  typedef Device  device_type ;

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType inode ) const
    { w(inode) = w(inode) + alpha * x(inode); }

  inline
  WAXPBY( const size_t  n ,
          const Scalar & arg_alpha ,
          const View< Array<Scalar,N> *, LayoutRight , Device > & arg_x ,
          const View< Array<Scalar,N> *, LayoutRight , Device > & arg_w )
    : w( arg_w ), x( arg_x )
    , alpha( arg_alpha )
  {
    parallel_for( n , *this );
  }
}; // WAXPBY

template < typename Scalar , unsigned N , class Device >
struct WAXPBY< void ,
               View< Array<Scalar,N> * , LayoutRight , Device > ,
               Scalar ,
               void ,
               View< Array<Scalar,N> * , LayoutRight , Device > >
{
private:

  const View< Array<Scalar,N> *, LayoutRight , Device , MemoryUnmanaged >  w ;
  const View< Array<Scalar,N> *, LayoutRight , Device , MemoryUnmanaged >  x ;
  const Scalar  beta ;

public:

  typedef Device  device_type ;

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType inode ) const
  { w(inode) = x(inode) + beta * w(inode); }

  inline
  WAXPBY( const size_t  n ,
          const View< Array<Scalar,N> *, LayoutRight , Device > & arg_x ,
          const Scalar & arg_beta ,
          const View< Array<Scalar,N> *, LayoutRight , Device > & arg_w )
    : w( arg_w ), x( arg_x )
    , beta( arg_beta )
  {
    parallel_for( n , *this );
  }
}; // WAXPBY

template < class Device >
struct WAXPBY< double ,
               View< double * , LayoutRight , Device > ,
               double ,
               View< double * , LayoutRight , Device > ,
               View< double * , LayoutRight , Device > >
{
private:

  const View< double *, LayoutRight , Device , MemoryUnmanaged >  w ;
  const View< double *, LayoutRight , Device , MemoryUnmanaged >  x ;
  const View< double *, LayoutRight , Device , MemoryUnmanaged >  y ;
  const double  alpha ;
  const double  beta ;

public:

  typedef Device  device_type ;

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType inode ) const
  { w(inode) = alpha * x(inode) + beta * y(inode); }

  inline
  WAXPBY( const size_t  n ,
          const double & arg_alpha ,
          const View< double *, LayoutRight , Device > & arg_x ,
          const double & arg_beta ,
          const View< double *, LayoutRight , Device > & arg_y ,
          const View< double *, LayoutRight , Device > & arg_w )
    : w( arg_w ), x( arg_x ), y( arg_y )
    , alpha( arg_alpha ), beta( arg_beta )
  {
    parallel_for( n , *this );
  }
}; // WAXPBY

template < class Device >
struct WAXPBY< double ,
               View< double * , LayoutRight , Device > ,
               void , void ,
               View< double * , LayoutRight , Device > >
{
private:

  const View< double *, LayoutRight , Device , MemoryUnmanaged >  w ;
  const View< double *, LayoutRight , Device , MemoryUnmanaged >  x ;
  const double  alpha ;

public:

  typedef Device  device_type ;

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType inode ) const
  { w(inode) += alpha * x(inode); }

  inline
  WAXPBY( const size_t  n ,
          const double & arg_alpha ,
          const View< double *, LayoutRight , Device > & arg_x ,
          const View< double *, LayoutRight , Device > & arg_w )
    : w( arg_w ), x( arg_x )
    , alpha( arg_alpha )
  {
    parallel_for( n , *this );
  }
}; // WAXPBY

template < class Device >
struct WAXPBY< void ,
               View< double * , LayoutRight , Device > ,
               double ,
               void ,
               View< double * , LayoutRight , Device > >
{
private:

  const View< double *, LayoutRight , Device , MemoryUnmanaged >  w ;
  const View< double *, LayoutRight , Device , MemoryUnmanaged >  x ;
  const double  beta ;

public:

  typedef Device  device_type ;

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType inode ) const
  { w(inode) = x(inode) + beta * w(inode); }

  inline
  WAXPBY( const size_t  n ,
          const View< double *, LayoutRight , Device > & arg_x ,
          const double & arg_beta ,
          const View< double *, LayoutRight , Device > & arg_w )
    : w( arg_w ), x( arg_x )
    , beta( arg_beta )
  {
    parallel_for( n , *this );
  }
}; // WAXPBY

//----------------------------------------------------------------------------

template< typename Scalar , unsigned N , class Device >
struct Dot< View< Array< Scalar , N > * , LayoutRight , Device > ,
            View< Array< Scalar , N > * , LayoutRight , Device > >
{
  typedef Device device_type ;
  typedef double value_type ;

  const View< Array< Scalar , N > *, LayoutRight , Device , MemoryUnmanaged >  x ;
  const View< Array< Scalar , N > *, LayoutRight , Device , MemoryUnmanaged >  y ;

  KOKKOSARRAY_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
  { update += source;    }

  KOKKOSARRAY_INLINE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType & i , value_type & update ) const
  { update += dot( x(i) , y(i) ); }

  Dot( const size_t n ,
       const View< Array< Scalar , N > * , LayoutRight , Device > & arg_x ,
       const View< Array< Scalar , N > * , LayoutRight , Device > & arg_y ,
       double & result )
  : x( arg_x ), y( arg_y )
  {
    result = parallel_reduce( n , *this );
  }
};

template< typename Scalar , unsigned N , class Device >
struct Dot< View< Array< Scalar , N > * , LayoutRight , Device > , void >
{
  typedef Device device_type ;
  typedef double value_type ;

  const View< Array< Scalar , N > *, LayoutRight , Device , MemoryUnmanaged >  x ;

  KOKKOSARRAY_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
  { update += source;    }

  KOKKOSARRAY_INLINE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType & i , value_type & update ) const
  { update += dot( x(i) , x(i) ); }

  Dot( const size_t n ,
       const View< Array< Scalar , N > * , LayoutRight , Device > & arg_x ,
       double & result )
  : x( arg_x )
  {
    result = parallel_reduce( n , *this );
  }
};

template< class Device >
struct Dot< View< double * , LayoutRight , Device > ,
            View< double * , LayoutRight , Device > >
{
  typedef Device device_type ;
  typedef double value_type ;

  const View< double *, LayoutRight , Device , MemoryUnmanaged >  x ;
  const View< double *, LayoutRight , Device , MemoryUnmanaged >  y ;

  KOKKOSARRAY_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
  { update += source;    }

  KOKKOSARRAY_INLINE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType & i , value_type & update ) const
  { update += dot( x(i) , y(i) ); }

  Dot( const size_t n ,
       const View< double * , LayoutRight , Device > & arg_x ,
       const View< double * , LayoutRight , Device > & arg_y ,
       double & result )
  : x( arg_x ), y( arg_y )
  {
    result = parallel_reduce( n , *this );
  }
};

template< class Device >
struct Dot< View< double * , LayoutRight , Device > , void >
{
  typedef Device device_type ;
  typedef double value_type ;

  const View< double *, LayoutRight , Device , MemoryUnmanaged >  x ;

  KOKKOSARRAY_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
  { update += source;    }

  KOKKOSARRAY_INLINE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const iType & i , value_type & update ) const
  { update += dot( x(i) , x(i) ); }

  Dot( const size_t n ,
       const View< double * , LayoutRight , Device > & arg_x ,
       double & result )
  : x( arg_x )
  {
    result = parallel_reduce( n , *this );
  }
};

}

//----------------------------------------------------------------------------

#if defined( __CUDACC__ )

namespace KokkosArray {

template < typename Scalar , unsigned N >
struct WAXPBY< Scalar ,
               View< Array<Scalar,N> * , LayoutRight , Cuda > ,
               Scalar ,
               View< Array<Scalar,N> * , LayoutRight , Cuda > ,
               View< Array<Scalar,N> * , LayoutRight , Cuda > >
{
private:

  const View< Array<Scalar,N> *, LayoutRight , Cuda , MemoryUnmanaged >  w ;
  const View< Array<Scalar,N> *, LayoutRight , Cuda , MemoryUnmanaged >  x ;
  const View< Array<Scalar,N> *, LayoutRight , Cuda , MemoryUnmanaged >  y ;
  const Scalar  alpha ;
  const Scalar  beta ;

public:

  __device__
  void operator()() const
  {
    if ( threadIdx.x < N ) {
      for ( int iRow = threadIdx.y + blockDim.y * blockIdx.x ;
                iRow < x.dimension_0() ; iRow += blockDim.y * gridDim.x ) {
        w(iRow,threadIdx.x) = alpha * x(iRow,threadIdx.x) + beta * y(iRow,threadIdx.x);
      }
    }
  }

  inline
  WAXPBY( const size_t  n ,
          const Scalar & arg_alpha ,
          const View< Array<Scalar,N> *, LayoutRight , Cuda > & arg_x ,
          const Scalar & arg_beta ,
          const View< Array<Scalar,N> *, LayoutRight , Cuda > & arg_y ,
          const View< Array<Scalar,N> *, LayoutRight , Cuda > & arg_w )
    : w( arg_w ), x( arg_x ), y( arg_y )
    , alpha( arg_alpha ), beta( arg_beta )
  {
    enum { W = Impl::CudaTraits::WarpSize };
    enum { NX = ( N + W - 1 ) / W };
    enum { NY = NX < 4 ? 16 : 1 };
    const size_t gMax = Impl::CudaTraits::UpperBoundGridCount ;

    const dim3 dBlock( W * NX , NY , 1 );
    const dim3 dGrid( std::min( ( n + NY - 1 ) / NY , gMax ) , 1 , 1 );

    Impl::cuda_parallel_launch_local_memory< WAXPBY ><<< dGrid , dBlock >>>( *this );
  }
}; // WAXPBY

template < typename Scalar , unsigned N >
struct WAXPBY< Scalar ,
               View< Array<Scalar,N> * , LayoutRight , Cuda > ,
               void , void ,
               View< Array<Scalar,N> * , LayoutRight , Cuda > >
{
private:

  const View< Array<Scalar,N> *, LayoutRight , Cuda , MemoryUnmanaged >  w ;
  const View< Array<Scalar,N> *, LayoutRight , Cuda , MemoryUnmanaged >  x ;
  const Scalar  alpha ;

public:

  __device__
  void operator()() const
  {
    if ( threadIdx.x < N ) {
      for ( int iRow = threadIdx.y + blockDim.y * blockIdx.x ;
                iRow < x.dimension_0() ; iRow += blockDim.y * gridDim.x ) {
        w(iRow,threadIdx.x) = alpha * x(iRow,threadIdx.x) + w(iRow,threadIdx.x);
      }
    }
  }

  inline
  WAXPBY( const size_t  n ,
          const Scalar & arg_alpha ,
          const View< Array<Scalar,N> *, LayoutRight , Cuda > & arg_x ,
          const View< Array<Scalar,N> *, LayoutRight , Cuda > & arg_w )
    : w( arg_w ), x( arg_x )
    , alpha( arg_alpha )
  {
    enum { W = Impl::CudaTraits::WarpSize };
    enum { NX = ( N + W - 1 ) / W };
    enum { NY = NX < 4 ? 16 : 1 };
    const size_t gMax = Impl::CudaTraits::UpperBoundGridCount ;

    const dim3 dBlock( W * NX , NY , 1 );
    const dim3 dGrid( std::min( ( n + NY - 1 ) / NY , gMax ) , 1 , 1 );

    Impl::cuda_parallel_launch_local_memory< WAXPBY ><<< dGrid , dBlock >>>( *this );
  }
}; // WAXPBY

template < typename Scalar , unsigned N >
struct WAXPBY< void ,
               View< Array<Scalar,N> * , LayoutRight , Cuda > ,
               Scalar ,
               void ,
               View< Array<Scalar,N> * , LayoutRight , Cuda > >
{
private:

  const View< Array<Scalar,N> *, LayoutRight , Cuda , MemoryUnmanaged >  w ;
  const View< Array<Scalar,N> *, LayoutRight , Cuda , MemoryUnmanaged >  x ;
  const Scalar  beta ;

public:

  __device__
  void operator()() const
  {
    if ( threadIdx.x < N ) {
      for ( int iRow = threadIdx.y + blockDim.y * blockIdx.x ;
                iRow < x.dimension_0() ; iRow += blockDim.y * gridDim.x ) {
        w(iRow,threadIdx.x) = x(iRow,threadIdx.x) + beta * w(iRow,threadIdx.x);
      }
    }
  }

  inline
  WAXPBY( const size_t  n ,
          const View< Array<Scalar,N> *, LayoutRight , Cuda > & arg_x ,
          const Scalar & arg_beta ,
          const View< Array<Scalar,N> *, LayoutRight , Cuda > & arg_w )
    : w( arg_w ), x( arg_x )
    , beta( arg_beta )
  {
    enum { W = Impl::CudaTraits::WarpSize };
    enum { NX = ( N + W - 1 ) / W };
    enum { NY = NX < 4 ? 16 : 1 };
    const size_t gMax = Impl::CudaTraits::UpperBoundGridCount ;

    const dim3 dBlock( W * NX , NY , 1 );
    const dim3 dGrid( std::min( ( n + NY - 1 ) / NY , gMax ) , 1 , 1 );

    Impl::cuda_parallel_launch_local_memory< WAXPBY ><<< dGrid , dBlock >>>( *this );
  }
}; // WAXPBY

//----------------------------------------------------------------------------

template< typename Scalar , unsigned N >
struct Dot< View< Array< Scalar , N > * , LayoutRight , Cuda > ,
            View< Array< Scalar , N > * , LayoutRight , Cuda > >
{
  typedef double value_type ;

  struct ValueOper {
    typedef double value_type ;

    KOKKOSARRAY_INLINE_FUNCTION
    static void join( volatile value_type & update ,
                      const volatile value_type & source )
      { update += source;    }

    KOKKOSARRAY_INLINE_FUNCTION
    static void init( value_type & update )
      { update = 0 ; }

    KOKKOSARRAY_INLINE_FUNCTION
    void operator()( const value_type & value ) const
      { *m_dev = value ; }

    ValueOper()
      : m_dev( Impl::CudaReduceResult<value_type>::device_pointer(1) ) {}

  private:
    value_type * const m_dev ;
  };

  const Impl::CudaReduceShared< ValueOper , ValueOper > m_reduce ;
  const View< Array< Scalar , N > *, LayoutRight , Cuda , MemoryUnmanaged >  x ;
  const View< Array< Scalar , N > *, LayoutRight , Cuda , MemoryUnmanaged >  y ;

  __device__
  void operator()(void) const
  {
    const int tid = threadIdx.x + blockDim.x * threadIdx.y ;

    value_type & update = m_reduce.value( tid );

    const int work_stride = blockDim.y * gridDim.x ;

    for ( int i = threadIdx.y + blockDim.y * blockIdx.x ; 
              i < x.dimension_0() ; i += work_stride ) {

      update += x(i,threadIdx.x) * y(i,threadIdx.x);
    }

    m_reduce( tid , blockIdx.x );
  }

  static inline
  int warp_count()
  {
    const Impl::CudaReduceSharedSizes<value_type> tmp( sizeof(value_type) );
    return tmp.warp_count ;
  }

  inline
  int block_count( const size_t n )
  {
    enum { W = Impl::CudaTraits::WarpSize };
    enum { NX = ( N + W - 1 ) / W };
    const size_t NY = warp_count() / NX ;
    const size_t gMax = Impl::CudaTraits::UpperBoundGridCount ;

    return std::min( ( n + NY - 1 ) / NY , gMax );
  }

  inline
  Dot( const size_t n ,
       const View< Array< Scalar , N > * , LayoutRight , Cuda > & arg_x ,
       const View< Array< Scalar , N > * , LayoutRight , Cuda > & arg_y ,
       double & result )
  : m_reduce( ValueOper() , 0 , block_count(n)), x( arg_x ), y( arg_y )
  {
    enum { W = Impl::CudaTraits::WarpSize };
    enum { NX = ( N + W - 1 ) / W };

    const size_t NY = warp_count() / NX ;
    const dim3 dBlock( W * NX , NY , 1 );
    const dim3 dGrid( block_count(n) , 1 , 1 );
    const int shmem = m_reduce.shmem_size();

    Impl::cuda_parallel_launch_local_memory< Dot ><<< dGrid , dBlock , shmem >>>( *this );

    result = Impl::CudaReduceResult<value_type>::return_to_host();
  }
};

template< typename Scalar , unsigned N >
struct Dot< View< Array< Scalar , N > * , LayoutRight , Cuda > , void >
{
  typedef View< Array< Scalar , N > * , LayoutRight , Cuda > vector_type ;
  Dot( const size_t n , const vector_type & arg_x , double & result )
  {
    Dot<vector_type,vector_type>( n , arg_x , arg_x , result );
  }
};

}

#endif

//----------------------------------------------------------------------------

namespace KokkosArray {

template< typename TypeAlpha ,
          class    TypeX ,
          typename TypeBeta ,
          class    TypeY ,
          class    TypeW >
void waxpby( const size_t count ,
             const TypeAlpha & alpha ,
             const TypeX     & x ,
             const TypeBeta  & beta ,
             const TypeY     & y ,
             const TypeW     & w )
{
  WAXPBY<TypeAlpha,TypeX,TypeBeta,TypeY,TypeW>(count,alpha,x,beta,y,w);
}

template< typename TypeAlpha ,
          class    TypeX ,
          class    TypeY >
void axpy( const size_t count ,
             const TypeAlpha & alpha ,
             const TypeX     & x ,
             const TypeY     & y )
{
  WAXPBY<TypeAlpha,TypeX,void,void,TypeY>(count,alpha,x,y);
}

template< class    TypeX ,
          typename TypeBeta ,
          class    TypeY >
void xpby( const size_t count ,
             const TypeX     & x ,
             const TypeBeta  & beta ,
             const TypeY     & y )
{
  WAXPBY<void,TypeX,TypeBeta,void,TypeY>(count,x,beta,y);
}

template< class TypeX >
double dot( const size_t count , const TypeX & x )
{
  double result ;
  Dot<TypeX,void>( count , x , result );
  return result ;
}

template< class TypeX , class TypeY >
double dot( const size_t count , const TypeX & x , const TypeY & y )
{
  double result ;
  Dot<TypeX,TypeY>( count , x , y , result );
  return result ;
}

}

//----------------------------------------------------------------------------

#endif

