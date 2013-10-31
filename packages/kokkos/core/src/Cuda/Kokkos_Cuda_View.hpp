/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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

#ifndef KOKKOS_CUDA_VIEW_HPP
#define KOKKOS_CUDA_VIEW_HPP

#include <cstring>

#if defined( __CUDACC__ )
#include <cuda_runtime.h>
#endif

#include <Kokkos_View.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Kokkos_CudaSpace.hpp>
#include <Kokkos_CudaTypes.hpp>
#include <Cuda/Kokkos_Cuda_abort.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<>
struct AssertShapeBoundsAbort< CudaSpace >
{
  KOKKOS_INLINE_FUNCTION
  static void apply( const size_t /* rank */ ,
                     const size_t /* n0 */ , const size_t /* n1 */ ,
                     const size_t /* n2 */ , const size_t /* n3 */ ,
                     const size_t /* n4 */ , const size_t /* n5 */ ,
                     const size_t /* n6 */ , const size_t /* n7 */ ,

                     const size_t /* arg_rank */ ,
                     const size_t /* i0 */ , const size_t /* i1 */ ,
                     const size_t /* i2 */ , const size_t /* i3 */ ,
                     const size_t /* i4 */ , const size_t /* i5 */ ,
                     const size_t /* i6 */ , const size_t /* i7 */ )
    {
      Kokkos::cuda_abort("Kokkos::View array bounds violation");
    }
};

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Cuda 5.0 <texture_types.h> defines 'cudaTextureObject_t'
// to be an 'unsigned long long'.  This chould change with
// future version of Cuda and this typedef would have to
// change accordingly.

#if defined( CUDA_VERSION ) && ( 500 <= CUDA_VERSION )

typedef enable_if<
  sizeof(::cudaTextureObject_t) == sizeof(const void *) ,
  ::cudaTextureObject_t >::type cuda_texture_object_type ;

cuda_texture_object_type
cuda_texture_object_attach(
  const cudaChannelFormatDesc & ,
  const void * const );

template< typename TextureType >
inline
cuda_texture_object_type
cuda_texture_object_attach( const void * const base_view_ptr )
{
  return cuda_texture_object_attach( cudaCreateChannelDesc<TextureType>() , base_view_ptr );
}

#else

typedef const void * cuda_texture_object_type ;

template< typename TextureType >
inline
cuda_texture_object_type
cuda_texture_object_attach( const void * const )
{ return 0 ; }

#endif

//----------------------------------------------------------------------------

template< typename ValueType >
struct CudaTextureFetch ;

/** \brief  Cuda texture fetch is limited to a subset of Cuda types.
 *          Map commonly used types to the required subset of Cuda types.
 */

template< typename ValueType >
struct CudaTextureFetch< const ValueType > {
private:

  cuda_texture_object_type  obj ;

public:

  const ValueType * ptr ;

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch() : obj( 0 ) , ptr( 0 ) {}

  KOKKOS_INLINE_FUNCTION
  ~CudaTextureFetch() {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( const CudaTextureFetch & rhs )
    : obj( rhs.obj ) , ptr( rhs.ptr ) {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = ( const CudaTextureFetch & rhs )
    { obj = rhs.obj ; ptr = rhs.ptr ; return *this ; }

  explicit
  CudaTextureFetch( ValueType * const base_view_ptr )
    : obj( cuda_texture_object_attach<ValueType>( base_view_ptr ) )
    , ptr( base_view_ptr ) {}

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  ValueType operator[]( const iType & i ) const
  {
    return ptr[ i ];
  }
};

template<>
struct CudaTextureFetch< const int > {
private:

  cuda_texture_object_type  obj ;

public:

  const int * ptr ;

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch() : obj( 0 ) , ptr( 0 ) {}

  KOKKOS_INLINE_FUNCTION
  ~CudaTextureFetch() {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( const CudaTextureFetch & rhs )
    : obj( rhs.obj ) , ptr( rhs.ptr ) {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = ( const CudaTextureFetch & rhs )
    { obj = rhs.obj ; ptr = rhs.ptr ; return *this ; }

  explicit
  CudaTextureFetch( const int * const base_view_ptr )
    : obj( cuda_texture_object_attach<int>( base_view_ptr ) )
    , ptr( base_view_ptr ) {}

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  int operator[]( const iType & i ) const
  {
#if defined( __CUDA_ARCH__ ) && ( 300 <= __CUDA_ARCH__ )
#ifdef KOKKOS_USE_LDG_INTRINSIC
    return _ldg(&ptr[i]);
#else
    return tex1Dfetch<int>( obj , i );
#endif
#else
    return ptr[ i ];
#endif
  }
};

template<>
struct CudaTextureFetch< const unsigned int > {
private:

  cuda_texture_object_type  obj ;

public:

  const unsigned int * ptr ;

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch() : obj( 0 ) , ptr( 0 ) {}

  KOKKOS_INLINE_FUNCTION
  ~CudaTextureFetch() {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( const CudaTextureFetch & rhs )
    : obj( rhs.obj ) , ptr( rhs.ptr ) {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = ( const CudaTextureFetch & rhs )
    { obj = rhs.obj ; ptr = rhs.ptr ; return *this ; }

  explicit
  CudaTextureFetch( const unsigned int * const base_view_ptr )
    : obj( cuda_texture_object_attach<unsigned int>( base_view_ptr ) )
    , ptr( base_view_ptr ) {}

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  unsigned int operator[]( const iType & i ) const
  {
#if defined( __CUDA_ARCH__ ) && ( 300 <= __CUDA_ARCH__ )
#ifdef KOKKOS_USE_LDG_INTRINSIC
    return _ldg(&ptr[i]);
#else
    return tex1Dfetch<unsigned int>( obj , i );
#endif
#else
    return ptr[ i ];
#endif
  }
};

template<>
struct CudaTextureFetch< const float > {
private:

  cuda_texture_object_type  obj ;

public:

  const float * ptr ;

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch() : obj( 0 ) , ptr( 0 ) {}

  KOKKOS_INLINE_FUNCTION
  ~CudaTextureFetch() {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( const CudaTextureFetch & rhs )
    : obj( rhs.obj ) , ptr( rhs.ptr ) {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = ( const CudaTextureFetch & rhs )
    { obj = rhs.obj ; ptr = rhs.ptr ; return *this ; }

  explicit
  CudaTextureFetch( const float * const base_view_ptr )
    : obj( cuda_texture_object_attach<float>( base_view_ptr ) )
    , ptr( base_view_ptr ) {}

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  float operator[]( const iType & i ) const
  {
#if defined( __CUDA_ARCH__ ) && ( 300 <= __CUDA_ARCH__ )
#ifdef KOKKOS_USE_LDG_INTRINSIC
    return _ldg(&ptr[i]);
#else
    return tex1Dfetch<float>( obj , i );
#endif
#else
    return ptr[ i ];
#endif
  }
};

template<>
struct CudaTextureFetch< const double > {
private:

  cuda_texture_object_type  obj ;

public:

  const double * ptr ;

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch() : obj( 0 ) , ptr( 0 ) {}

  KOKKOS_INLINE_FUNCTION
  ~CudaTextureFetch() {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( const CudaTextureFetch & rhs )
    : obj( rhs.obj ) , ptr( rhs.ptr ) {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = ( const CudaTextureFetch & rhs )
    { obj = rhs.obj ; ptr = rhs.ptr ; return *this ; }

  explicit
  CudaTextureFetch( const double * const base_view_ptr )
    : obj( cuda_texture_object_attach<int2>( base_view_ptr ) )
    , ptr( base_view_ptr ) {}

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  double operator[]( const iType & i ) const
  {
#if defined( __CUDA_ARCH__ ) && ( 300 <= __CUDA_ARCH__ )
#ifdef KOKKOS_USE_LDG_INTRINSIC
    return _ldg(&ptr[i]);
#else
    int2 v = tex1Dfetch<int2>( obj , i );
    return __hiloint2double(v.y, v.x);
#endif
#else
    return ptr[ i ];
#endif
  }
};

template<>
struct CudaTextureFetch< const double2 > {
private:

  cuda_texture_object_type  obj ;

public:

  const double2 * ptr ;

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch() : obj( 0 ) , ptr( 0 ) {}

  KOKKOS_INLINE_FUNCTION
  ~CudaTextureFetch() {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( const CudaTextureFetch & rhs )
    : obj( rhs.obj ) , ptr( rhs.ptr ) {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = ( const CudaTextureFetch & rhs )
    { obj = rhs.obj ; ptr = rhs.ptr ; return *this ; }

  explicit
  CudaTextureFetch( const double2 * const base_view_ptr )
    : obj( cuda_texture_object_attach<int4>( base_view_ptr ) )
    , ptr( base_view_ptr ) {}

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  double2 operator[]( const iType & i ) const
  {
#if defined( __CUDA_ARCH__ ) && ( 300 <= __CUDA_ARCH__ )
#ifdef KOKKOS_USE_LDG_INTRINSIC
    return _ldg(&ptr[i]);
#else
    int4 v = tex1Dfetch<int4>(tex_obj , idx);
    double2 retval = { __hiloint2double(v.y, v.x) , __hiloint2double(v.w, v.z) };
    return retval ;
#endif
#else
    return ptr[ i ];
#endif
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

struct CudaTexture {};

#if defined( CUDA_VERSION ) && ( 500 <= CUDA_VERSION )

/** \brief  Replace LayoutDefault specialization */
template< typename ScalarType , class Rank , class RankDynamic >
struct ViewSpecialize< const ScalarType , const ScalarType ,
                       LayoutLeft , Rank , RankDynamic ,
                       CudaSpace , MemoryRandomRead >
{ typedef CudaTexture type ; };

template< typename ScalarType , class Rank , class RankDynamic >
struct ViewSpecialize< const ScalarType , const ScalarType ,
                       LayoutRight , Rank , RankDynamic ,
                       CudaSpace , MemoryRandomRead >
{ typedef CudaTexture type ; };

/** \brief Scalar View matching **/
template< typename ScalarType >
struct ViewSpecialize< const ScalarType , const ScalarType ,
                       LayoutLeft , unsigned_<0> , unsigned_<0> ,
                       CudaSpace , MemoryRandomRead >
{ typedef CudaTexture type ; };

template< typename ScalarType >
struct ViewSpecialize< const ScalarType , const ScalarType ,
                       LayoutRight , unsigned_<0> , unsigned_<0> ,
                       CudaSpace , MemoryRandomRead >
{ typedef CudaTexture type ; };

#endif

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< CudaTexture , CudaTexture , void >
{
  /** \brief Assign compatible views */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,CudaTexture> & dst ,
                  const View<ST,SL,SD,SM,CudaTexture> & src ,
                  const typename enable_if<(
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::value
                  ) >::type * = 0 )
  {
    typedef View<DT,DL,DD,DM,CudaTexture> DstViewType ;

    typedef typename DstViewType::shape_type    shape_type ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    dst.m_texture  = src.m_texture ;
    dst.m_stride   = src.m_stride ;

    shape_type::assign( dst.m_shape,
                        src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                        src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );
  }
};


template<>
struct ViewAssignment< CudaTexture , LayoutDefault , void >
{
  /** \brief Assign compatible views */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  inline
  ViewAssignment(       View<DT,DL,DD,DM,CudaTexture> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
                  const typename enable_if<(
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::value
                  )>::type * = 0 )
  {
    typedef View<DT,DL,DD,DM,CudaTexture> DstViewType ;

    typedef typename DstViewType::shape_type  shape_type ;
    typedef typename DstViewType::scalar_type scalar_type ;
    typedef typename DstViewType::stride_type stride_type ;

    dst.m_texture = CudaTextureFetch< scalar_type >( src.m_ptr_on_device );

    shape_type::assign( dst.m_shape,
                        src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                        src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

    stride_type::assign( dst.m_stride , src.m_stride.value );
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
template< class T , class L, class D , class M >
class View< T , L , D , M , Impl::CudaTexture >
  : public ViewTraits< T , L , D , M >
{
public:

  typedef ViewTraits< T , L , D , M > traits ;

private:

  template< class , class , class > friend struct Impl::ViewAssignment ;

  typedef Impl::LayoutStride< typename traits::shape_type ,
                              typename traits::array_layout > stride_type ;

  Impl::CudaTextureFetch<typename traits::scalar_type > m_texture ;
  typename traits::shape_type           m_shape ;
  stride_type                           m_stride ;

public:

  typedef Impl::CudaTexture specialize ;

  typedef View< typename traits::const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type::host_mirror_device_type ,
                void > HostMirror ;

  enum { Rank = traits::rank };

  KOKKOS_INLINE_FUNCTION typename traits::shape_type shape() const { return m_shape ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_0() const { return m_shape.N0 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_1() const { return m_shape.N1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_2() const { return m_shape.N2 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_3() const { return m_shape.N3 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_4() const { return m_shape.N4 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_5() const { return m_shape.N5 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_6() const { return m_shape.N6 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_7() const { return m_shape.N7 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type size() const
  {
    return   m_shape.N0
           * m_shape.N1
           * m_shape.N2
           * m_shape.N3
           * m_shape.N4
           * m_shape.N5
           * m_shape.N6
           * m_shape.N7
           ;
  }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  typename traits::size_type dimension( const iType & i ) const
    { return Impl::dimension( m_shape , i ); }

  //------------------------------------

  View() : m_texture()
   {
     traits::shape_type::assign(m_shape,0,0,0,0,0,0,0,0);
     stride_type::assign( m_stride , 0 );
   }

  ~View() {}

  View( const View & rhs )
    : m_texture( rhs.m_texture )
    , m_stride(  rhs.m_stride )
    { m_shape = rhs.m_shape ; }

  View & operator = ( const View & rhs )
    {
      (void)Impl::ViewAssignment< Impl::CudaTexture , Impl::CudaTexture >( *this , rhs );
      return *this ;
    }

  template< class RT , class RL, class RD , class RM , class RS >
  View( const View<RT,RL,RD,RM,RS> & rhs )
    : m_texture(0)
    {
      Impl::ViewAssignment< Impl::CudaTexture , RS >( *this , rhs );
    }

  template< class RT , class RL, class RD, class RM , class RS >
  View & operator = ( const View<RT,RL,RD,RM,RS> & rhs )
    {
      Impl::ViewAssignment< Impl::CudaTexture , RS >( *this , rhs );
      return *this ;
    }

  //------------------------------------

  KOKKOS_INLINE_FUNCTION
  bool is_null() const { return 0 == m_texture.ptr ; }

  //------------------------------------
  // Rank = 1 access operators:

  template < typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type , traits , LayoutLeft , 1 , iType0 >::type
    operator[] ( const iType0 & i0 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      return m_texture[ i0 ];
    }

  template < typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type , traits , LayoutRight , 1 , iType0 >::type
    operator[] ( const iType0 & i0 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      return m_texture[ i0 ];
    }

  template < typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type , traits , LayoutLeft , 1 , iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      return m_texture[ i0 ];
    }

  template < typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type , traits , LayoutRight , 1 , iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      return m_texture[ i0 ];
    }

  //------------------------------------
  // Layout left:


  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type , traits, LayoutLeft, 2, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ i0 + m_stride.value * i1 ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type ,
                                      traits, LayoutLeft, 3, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0,i1,i2 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ i0 + m_stride.value * (
                        i1 + m_shape.N1 * i2 ) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type ,
                                      traits, LayoutLeft, 4, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0,i1,i2,i3 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ i0 + m_stride.value * (
                        i1 + m_shape.N1 * (
                        i2 + m_shape.N2 * i3 )) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type ,
                                      traits, LayoutLeft, 5, iType0, iType1, iType2, iType3, iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0,i1,i2,i3,i4 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ i0 + m_stride.value * (
                        i1 + m_shape.N1 * (
                        i2 + m_shape.N2 * (
                        i3 + m_shape.N3 * i4 ))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type ,
                                      traits, LayoutLeft, 6, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0,i1,i2,i3,i4,i5 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ i0 + m_stride.value * (
                        i1 + m_shape.N1 * (
                        i2 + m_shape.N2 * (
                        i3 + m_shape.N3 * (
                        i4 + m_shape.N4 * i5 )))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type ,
                                      traits, LayoutLeft, 7, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_shape, i0,i1,i2,i3,i4,i5,i6 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ i0 + m_stride.value * (
                        i1 + m_shape.N1 * (
                        i2 + m_shape.N2 * (
                        i3 + m_shape.N3 * (
                        i4 + m_shape.N4 * (
                        i5 + m_shape.N5 * i6 ))))) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 , typename iType7 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type ,
                                      traits, LayoutLeft, 8, iType0, iType1, iType2, iType3, iType4, iType5, iType6, iType7 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , const iType7 & i7 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_8( m_shape, i0,i1,i2,i3,i4,i5,i6,i7 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ i0 + m_stride.value * (
                        i1 + m_shape.N1 * (
                        i2 + m_shape.N2 * (
                        i3 + m_shape.N3 * (
                        i4 + m_shape.N4 * (
                        i5 + m_shape.N5 * (
                        i6 + m_shape.N6 * i7 )))))) ];
    }


  //------------------------------------
  // Layout right:


  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type ,
                                      traits, LayoutRight, 2, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ i1 + i0 * m_stride.value ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type ,
                                      traits, LayoutRight, 3, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_shape, i0,i1,i2 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ i2 + m_shape.N2 * i1 + i0 * m_stride.value ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type ,
                                      traits, LayoutRight, 4, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_shape, i0,i1,i2,i3 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ i3 + m_shape.N3 * (
                        i2 + m_shape.N2 * (
                        i1 )) + i0 * m_stride.value ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type ,
                                      traits, LayoutRight, 5, iType0, iType1, iType2, iType3, iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_shape, i0,i1,i2,i3,i4 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ i4 + m_shape.N4 * (
                        i3 + m_shape.N3 * (
                        i2 + m_shape.N2 * (
                        i1 ))) + i0 * m_stride.value ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type ,
                                      traits, LayoutRight, 6, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_shape, i0,i1,i2,i3,i4,i5 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ i5 + m_shape.N5 * (
                        i4 + m_shape.N4 * (
                        i3 + m_shape.N3 * (
                        i2 + m_shape.N2 * (
                        i1 )))) + i0 * m_stride.value ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type ,
                                      traits, LayoutRight, 7, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_shape, i0,i1,i2,i3,i4,i5,i6 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ i6 + m_shape.N6 * (
                        i5 + m_shape.N5 * (
                        i4 + m_shape.N4 * (
                        i3 + m_shape.N3 * (
                        i2 + m_shape.N2 * (
                        i1 ))))) + i0 * m_stride.value ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 , typename iType7 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::scalar_type ,
                                      traits, LayoutRight, 8, iType0, iType1, iType2, iType3, iType4, iType5, iType6, iType7 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , const iType7 & i7 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_8( m_shape, i0,i1,i2,i3,i4,i5,i6,i7 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ i7 + m_shape.N7 * (
                        i6 + m_shape.N6 * (
                        i5 + m_shape.N5 * (
                        i4 + m_shape.N4 * (
                        i3 + m_shape.N3 * (
                        i2 + m_shape.N2 * (
                        i1 )))))) + i0 * m_stride.value ];
    }

  //------------------------------------

  KOKKOS_INLINE_FUNCTION
  typename traits::scalar_type * ptr_on_device() const { return m_texture.ptr ; }

  // Stride of physical storage, dimensioned to at least Rank
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
  {
    enum { is_left = Impl::is_same< typename traits::array_layout , LayoutLeft >::value };

    if ( 1 == Rank ) {
      s[0] = 1 ;
    }
    else if ( is_left ) {
      s[0] = 1 ;
      s[1] = m_stride.value ;
      for ( int i = 2 ; i < Rank ; ++i ) { s[i] = s[i-1] * dimension(i-1); }
    }
    else {
      s[0] = m_stride.value ;
      s[Rank-1] = 1 ;
      for ( int i = Rank - 2 ; 0 < i ; --i ) { s[i] = s[i+1] * dimension(i+1); }
    }
  }
};

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_CUDA_VIEW_HPP */

