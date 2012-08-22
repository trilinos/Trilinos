/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_BLAS_KERNELS_HPP
#define KOKKOSARRAY_BLAS_KERNELS_HPP

namespace KokkosArray {

template< class ConstVectorType ,
          class Device = typename ConstVectorType::device_type >
struct Dot ;

template< class ConstVectorType ,
          class Device = typename ConstVectorType::device_type >
struct DotSingle ;

template< class ConstScalarType ,
          class VectorType ,
          class Device = typename VectorType::device_type >
struct Scale ;

template< class ConstScalarType ,
          class ConstVectorType ,
          class VectorType ,
          class Device = typename VectorType::device_type >
struct AXPBY ;

/** \brief  Y = alpha * X + beta * Y */
template< class ConstScalarType ,
          class ConstVectorType ,
          class      VectorType >
void axpby( const ConstScalarType & alpha ,
            const ConstVectorType & X ,
            const ConstScalarType & beta ,
            const      VectorType & Y )
{
  typedef AXPBY< ConstScalarType , ConstVectorType , VectorType > functor ;

  parallel_for( Y.dimension_0() , functor( alpha , X , beta , Y ) );
}

/** \brief  Y *= alpha */
template< class ConstScalarType ,
          class      VectorType >
void scale( const ConstScalarType & alpha , const VectorType & Y )
{
  typedef Scale< ConstScalarType , VectorType > functor ;

  parallel_for( Y.dimension_0() , functor( alpha , Y ) );
}

template< class ConstVectorType ,
          class Finalize >
void dot( const ConstVectorType & X ,
          const ConstVectorType & Y ,
          const Finalize & finalize )
{
  typedef Dot< ConstVectorType >  functor ;

  parallel_reduce( X.dimension_0() , functor( X , Y ) , finalize );
}

template< class ConstVectorType ,
          class Finalize >
void dot( const ConstVectorType & X ,
          const Finalize & finalize )
{
  typedef DotSingle< ConstVectorType >  functor ;

  parallel_reduce( X.dimension_0() , functor( X ) , finalize );
}

} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_BLAS_KERNELS_HPP */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOSARRAY_MACRO_DEVICE ) && defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION )

namespace KokkosArray {

template< class Type >
struct Dot< Type , KOKKOSARRAY_MACRO_DEVICE >
{
  typedef KOKKOSARRAY_MACRO_DEVICE device_type ;

  typedef typename 
    Impl::StaticAssertSame< Impl::unsigned_< 1 > ,
                            Impl::unsigned_< Type::Rank > >::type ok_rank ;


  typedef typename 
    Impl::StaticAssertSame< device_type ,
                            typename Type::device_type >::type ok_device ;

  typedef double value_type ;

#if 1
  typename Type::const_type X ;
  typename Type::const_type Y ;
#else
  Type X ;
  Type Y ;
#endif
  
  Dot( const Type & arg_x , const Type & arg_y )
    : X(arg_x) , Y(arg_y) { }

  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  void operator()( int i , value_type & update ) const 
    { update += X[i] * Y[i]; }
    
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
    { update += source; }
    
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  static void init( value_type & update )
    { update = 0 ; }
};

template< class Type >
struct DotSingle< Type , KOKKOSARRAY_MACRO_DEVICE >
{
  typedef KOKKOSARRAY_MACRO_DEVICE device_type ;

  typedef typename 
    Impl::StaticAssertSame< Impl::unsigned_< 1 > ,
                            Impl::unsigned_< Type::Rank > >::type ok_rank ;

  typedef typename
    Impl::StaticAssertSame< device_type ,
                            typename Type::device_type >::type ok_device ;

  typedef double value_type ;

#if 1
  typename Type::const_type X ;
#else
  Type X ;
#endif

  DotSingle( const Type & arg_x ) : X(arg_x) {}

  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  void operator()( int i , value_type & update ) const 
    {
      const typename Type::value_type x = X[i]; update += x * x ;
    }
    
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
    { update += source; }
    
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  static void init( value_type & update )
    { update = 0 ; }
};


template< class ScalarType , class VectorType >
struct Scale< ScalarType , VectorType , KOKKOSARRAY_MACRO_DEVICE >
{
  typedef KOKKOSARRAY_MACRO_DEVICE device_type ;

  typedef typename
    Impl::StaticAssertSame< device_type ,
                            typename ScalarType::device_type >::type
      ok_scalar_device ;

  typedef typename
    Impl::StaticAssertSame< device_type ,
                            typename VectorType::device_type >::type
      ok_vector_device ;

  typedef typename
    Impl::StaticAssertSame< Impl::unsigned_< 0 > ,
                            Impl::unsigned_< ScalarType::Rank > >::type
      ok_scalar_rank ;

  typedef typename 
    Impl::StaticAssertSame< Impl::unsigned_< 1 > ,
                            Impl::unsigned_< VectorType::Rank > >::type
      ok_vector_rank ;

#if 1
  typename ScalarType::const_type alpha ;
#else
  ScalarType alpha ;
#endif

  VectorType Y ;

  Scale( const ScalarType & arg_alpha , const VectorType & arg_Y )
    : alpha( arg_alpha ), Y( arg_Y ) {}

  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  void operator()( int i ) const 
    {
      Y[i] *= *alpha ;
    }
};


template< class ScalarType ,
          class ConstVectorType ,
          class VectorType >
struct AXPBY< ScalarType , ConstVectorType , VectorType , KOKKOSARRAY_MACRO_DEVICE >
{
  typedef KOKKOSARRAY_MACRO_DEVICE device_type ;

  typedef typename
    Impl::StaticAssertSame< device_type ,
                            typename ScalarType::device_type >::type
      ok_scalar_device ;

  typedef typename
    Impl::StaticAssertSame< device_type ,
                            typename ConstVectorType::device_type >::type
      ok_const_vector_device ;

  typedef typename
    Impl::StaticAssertSame< device_type ,
                            typename VectorType::device_type >::type
      ok_vector_device ;

  typedef typename
    Impl::StaticAssertSame< Impl::unsigned_< 0 > ,
                            Impl::unsigned_< ScalarType::Rank > >::type
      ok_scalar_rank ;

  typedef typename
    Impl::StaticAssertSame< Impl::unsigned_< 1 > ,
                            Impl::unsigned_< ConstVectorType::Rank > >::type
      ok_const_vector_rank ;

  typedef typename
    Impl::StaticAssertSame< Impl::unsigned_< 1 > ,
                            Impl::unsigned_< VectorType::Rank > >::type
      ok_vector_rank ;

#if 1
  typename ScalarType::const_type alpha , beta ;
  typename ConstVectorType::const_type X ;
#else
  ScalarType alpha , beta ;
  ConstVectorType X ;
#endif

  VectorType Y ;

  AXPBY( const ScalarType      & arg_alpha ,
         const ConstVectorType & arg_X ,
         const ScalarType      & arg_beta ,
         const VectorType      & arg_Y )
    : alpha( arg_alpha ), beta( arg_beta ), X( arg_X ), Y( arg_Y ) {}

  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  void operator()( int i ) const 
    {
      Y[i] = *alpha * X[i] + *beta * Y[i] ;
    }
};

} /* namespace KokkosArray */

#endif /* #if defined( KOKKOSARRAY_MACRO_DEVICE ) &&
              defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION ) */

