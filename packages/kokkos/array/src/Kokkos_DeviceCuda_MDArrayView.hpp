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

#ifndef KOKKOS_DEVICECUDA_MDARRAYVIEW_HPP
#define KOKKOS_DEVICECUDA_MDARRAYVIEW_HPP

#include <Kokkos_MDArrayView.hpp>
#include <Kokkos_DeviceCuda.hpp>

#include <impl/Kokkos_ArrayBounds.hpp>
#include <DeviceCuda/Kokkos_DeviceCuda_ParallelDriver.hpp>
#include <DeviceCuda/Kokkos_DeviceCuda_DeepCopy.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {

template< typename ValueType >
class MDArrayView< ValueType , DeviceCuda , DeviceCuda::mdarray_map >
{
public:
  typedef ValueType                value_type ;
  typedef DeviceCuda               device_type ;
  typedef DeviceCuda::mdarray_map  mdarray_map ;
  typedef DeviceCuda::size_type    size_type ;

  typedef MDArrayView< value_type , Serial< HostMemory , mdarray_map > > HostView ;

  /*------------------------------------------------------------------*/
  /** \brief  Not contiguous due to the need to pad for memory alignment */
  enum { Contiguous = false };

  /** \brief  Query rank of the array */
  inline
  __device__ __host__
  size_type rank() const { return m_rank ; }

  /** \brief  Query dimension of the given ordinate of the array */
  template < typename iType >
  inline
  __device__ __host__
  size_type dimension( const iType & rank_ordinate ) const
  { return m_dims[ rank_ordinate ]; }

  template < typename iType >
  inline
  __device__ __host__
  void dimensions( iType * const dims ) const
  { for ( size_type i = 0 ; i < m_rank ; ++i ) { dims[i] = m_dims[i] ; } }

  /** \brief  Query cardinality of the array */
  inline
  __device__ __host__
  size_type size() const
  {
    size_type n = m_dims[0] ;
    for ( size_type i = 1 ; i < m_rank ; ++i ) { n *= m_dims[i] ; }
    return n ;
  }

  /*------------------------------------------------------------------*/

#if 0
  /** \brief  Because memory is contiguous this is exposed */
  inline
  __device__
  value_type * ptr_on_device() const { return m_memory.ptr_on_device(); }
#endif

  /*------------------------------------------------------------------*/
  /** \brief  Query value of a rank 8 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 , typename iType7 >
  inline
  __device__
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 , const iType7 & i7 ) const
  { return m_memory.ptr_on_device()[
           ( iP + m_stride  * ( i1 + m_dims[1] *
           ( i2 + m_dims[2] * ( i3 + m_dims[3] *
           ( i4 + m_dims[4] * ( i5 + m_dims[5] *
           ( i6 + m_dims[6] * ( i7 )))))))) ]; }

  /** \brief  Query value of a rank 7 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 >
  inline
  __device__
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 ) const
  { return m_memory.ptr_on_device()[
           ( iP + m_stride  * ( i1 + m_dims[1] *
           ( i2 + m_dims[2] * ( i3 + m_dims[3] *
           ( i4 + m_dims[4] * ( i5 + m_dims[5] *
           ( i6 ))))))) ]; }

  /** \brief  Query value of a rank 6 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  inline
  __device__
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ) const
  { return m_memory.ptr_on_device()[
           ( iP + m_stride  * ( i1 + m_dims[1] *
           ( i2 + m_dims[2] * ( i3 + m_dims[3] *
           ( i4 + m_dims[4] * ( i5 )))))) ]; }

  /** \brief  Query value of a rank 5 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 >
  inline
  __device__
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 ) const
  { return m_memory.ptr_on_device()[
           ( iP + m_stride  * ( i1 + m_dims[1] *
           ( i2 + m_dims[2] * ( i3 + m_dims[3] *
           ( i4 ))))) ]; }

  /** \brief  Query value of a rank 4 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 >
  inline
  __device__
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ) const
  { return m_memory.ptr_on_device()[
           ( iP + m_stride  * ( i1 + m_dims[1] *
           ( i2 + m_dims[2] * ( i3 )))) ]; }

  /** \brief  Query value of a rank 3 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 >
  inline
  __device__
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 ) const
  { return m_memory.ptr_on_device()[
           ( iP + m_stride * ( i1 + m_dims[1] * ( i2 ))) ]; }

  /** \brief  Query value of a rank 2 array */
  template< typename iTypeP , typename iType1 >
  inline
  __device__
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ) const
  { return m_memory.ptr_on_device()[ ( iP + m_stride * i1 ) ]; }

  /** \brief  Query value of a rank 1 array */
  template< typename iTypeP >
  inline
  __device__
  value_type & operator()( const iTypeP & iP ) const
  { return m_memory.ptr_on_device()[ iP ]; }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  inline
  __device__ __host__
  MDArrayView() : m_memory(), m_rank(0), m_stride(0)
  {
    m_dims[0] = 0 ; m_dims[1] = 0 ; m_dims[2] = 0 ; m_dims[3] = 0 ;
    m_dims[4] = 0 ; m_dims[5] = 0 ; m_dims[6] = 0 ; m_dims[7] = 0 ;
  }

  /** \brief  Construct a view of the array */
  inline
  __device__ __host__
  MDArrayView( const MDArrayView & rhs )
    : m_memory(), m_rank( rhs.m_rank ), m_stride( rhs.m_stride )
    {
      m_dims[0] = rhs.m_dims[0] ; m_dims[1] = rhs.m_dims[1] ;
      m_dims[2] = rhs.m_dims[2] ; m_dims[3] = rhs.m_dims[3] ;
      m_dims[4] = rhs.m_dims[4] ; m_dims[5] = rhs.m_dims[5] ;
      m_dims[6] = rhs.m_dims[6] ; m_dims[7] = rhs.m_dims[7] ;
      device_type::assign_memory_view( m_memory , rhs.m_memory);
    }

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  inline
  __device__ __host__
  MDArrayView & operator = ( const MDArrayView & rhs )
    {
      m_rank   = rhs.m_rank ;
      m_stride = rhs.m_stride ;
      m_dims[0] = rhs.m_dims[0] ; m_dims[1] = rhs.m_dims[1] ;
      m_dims[2] = rhs.m_dims[2] ; m_dims[3] = rhs.m_dims[3] ;
      m_dims[4] = rhs.m_dims[4] ; m_dims[5] = rhs.m_dims[5] ;
      device_type::assign_memory_view( m_memory , rhs.m_memory );
      return *this ;
    }
  
  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  inline
  __device__ __host__
  ~MDArrayView() {}

  /*------------------------------------------------------------------*/
  inline
  __device__ __host__
  operator bool () const
  { return m_memory.operator bool(); }

  inline
  __device__ __host__
  bool operator == ( const MDArrayView & rhs ) const
  { return m_memory.operator == ( rhs.m_memory ); }

  inline
  __device__ __host__
  bool operator != ( const MDArrayView & rhs ) const
  { return m_memory.operator != ( rhs.m_memory ); }

  /*------------------------------------------------------------------*/

private:

  MemoryView< value_type , device_type >  m_memory ;
  size_type                               m_rank ;
  size_type                               m_stride ;
  size_type                               m_dims[ MDArrayMaxRank ];

  inline
  MDArrayView( const std::string & label ,
               size_t nP , size_t n1 , size_t n2 , size_t n3 ,
               size_t n4 , size_t n5 , size_t n6 , size_t n7 )
    : m_memory()
    , m_rank( Impl::mdarray_deduce_rank( nP, n1,n2, n3, n4, n5, n6, n7 ) )
    , m_stride( 0 )
    {
      m_dims[0] = nP ; m_dims[1] = n1 ; m_dims[2] = n2 ; m_dims[3] = n3 ;
      m_dims[4] = n4 ; m_dims[5] = n5 ; m_dims[6] = n6 ; m_dims[7] = n7 ;

      // For optimal global memory access performance align on the word * warp boundary.
      // Could use half warp size for slightly degraded performance.
      // Other memory alignments can result in severe performance degradation.
      // See section 5.3.2.1.2 of the CUDA C Programming Guide, Version 3.2.
      enum { NP_SIZE_ALIGN = sizeof(size_type) * Impl::DeviceCudaTraits::WarpSize };

      // Round up nP for optimal global memory access alignment 
      //  0 == ( sizeof(value_type) * nP ) % NP_SIZE_ALIGN

      size_type nP_size = sizeof(value_type) * nP ;

      if ( nP_size % NP_SIZE_ALIGN ) {
        nP_size += NP_SIZE_ALIGN - nP_size % NP_SIZE_ALIGN ;
      }
      m_stride = nP_size / sizeof(value_type);

      size_type size = m_stride ;
      for ( size_type r = 1 ; r < rank() ; ++r ) {
        size *= dimension(r);
      }

      device_type::allocate_memory_view( m_memory , size , label );

      parallel_for( size , Impl::AssignContiguous<value_type,device_type>( m_memory.ptr_on_device() , 0 ) );
    }

  template< typename V , class D , class M >
  friend
  MDArrayView< V , D , M >
  create_labeled_mdarray( const std::string & label ,
                          size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                          size_t n4 , size_t n5 , size_t n6 , size_t n7 );

  template< typename V , class DeviceDst , class MapDst , bool ,
                         class DeviceSrc , class MapSrc , bool >
  friend
  class Impl::MDArrayDeepCopy ;

  template< typename V , class D , class MapDst , class MapSrc , unsigned R >
  friend
  class Impl::MDArrayDeepCopyFunctor ;
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  Deep copy device to device */
template< typename ValueType >
class MDArrayDeepCopy< ValueType , DeviceCuda , DeviceCuda::mdarray_map , false ,
                                   DeviceCuda , DeviceCuda::mdarray_map , false >
{
public:
  typedef DeviceCuda            device_type ;
  typedef DeviceCuda::size_type size_type ;

  typedef MDArrayView< ValueType , DeviceCuda , DeviceCuda::mdarray_map > array_type ;

        ValueType * const dst ;
  const ValueType * const src ;

  MDArrayDeepCopy( ValueType * arg_dst , const ValueType * arg_src )
    : dst( arg_dst ), src( arg_src ) {}

  inline
  __device__
  void operator()( size_type iwork ) const
  { dst[iwork] = src[iwork] ; }

  static void run( const array_type & dst , const array_type & src )
  {
    Impl::mdarray_require_equal_dimension( dst , src );

    size_type n = dst.m_stride ;
    for ( size_type r = 1 ; r < dst.m_rank ; ++r ) { n *= dst.m_dims[r] ; }

    parallel_for( n , MDArrayDeepCopy( dst.m_memory.ptr_on_device() ,
                                       src.m_memory.ptr_on_device() ) );

  }
};

/** \brief  Copy Host to Cuda specialization */
template< typename ValueType >
class MDArrayDeepCopy< ValueType ,
                       DeviceCuda , DeviceCuda::mdarray_map , false ,
                       DeviceHost , DeviceCuda::mdarray_map , true >
{
public:
  typedef DeviceCuda::size_type                          size_type ;
  typedef MDArrayView< ValueType , DeviceCuda , DeviceCuda::mdarray_map > dst_type ;
  typedef MDArrayView< ValueType , DeviceHost , DeviceCuda::mdarray_map > src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    Impl::mdarray_require_equal_dimension( dst , src );

    const size_type length = dst.m_dims[0] ;
    size_type count = 1 ;
    for ( size_type r = 1 ; r < dst.m_rank ; ++r ) {
      count *= dst.m_dims[r] ;
    }

    for ( size_type i = 0 ; i < count ; ++i ) {
      Impl::copy_to_cuda_from_host(
        dst.m_memory.ptr_on_device() + dst.m_stride * i ,
        src.m_memory.ptr_on_device() + length * i ,
        sizeof(ValueType) , length );
    }
  }
};


/** \brief  Copy Cuda to Host specialization */
template< typename ValueType >
class MDArrayDeepCopy< ValueType ,
                       DeviceHost , DeviceCuda::mdarray_map , true ,
                       DeviceCuda , DeviceCuda::mdarray_map , false >
{
public:
  typedef DeviceCuda::size_type                          size_type ;
  typedef MDArrayView< ValueType , DeviceHost , DeviceCuda::mdarray_map > dst_type ;
  typedef MDArrayView< ValueType , DeviceCuda , DeviceCuda::mdarray_map > src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    Impl::mdarray_require_equal_dimension( dst , src );

    const size_type length = src.m_dims[0] ;
    size_type count = 1 ;
    for ( size_type r = 1 ; r < src.m_rank ; ++r ) {
      count *= src.m_dims[r] ;
    }

    for ( size_type i = 0 ; i < count ; ++i ) {
      Impl::copy_to_host_from_cuda(
        dst.m_memory.ptr_on_device() + length * i ,
        src.m_memory.ptr_on_device() + src.m_stride * i ,
        sizeof(ValueType) , length );
    }
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_DEVICECUDA_MDARRAYVIEW_HPP */


