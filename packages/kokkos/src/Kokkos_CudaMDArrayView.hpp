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

#ifndef KOKKOS_CUDAMDARRAYVIEW_HPP
#define KOKKOS_CUDAMDARRAYVIEW_HPP

#include <Kokkos_MDArrayViewHelper.hpp>
#include <Kokkos_CudaDevice.hpp>
#include <Kokkos_CudaMap.hpp>

namespace Kokkos {

template< typename ValueType , class DeviceMapType > class MDArrayView ;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** \brief  Multidimensional array allocated and mapped
 *          onto a compute device.
 *
 *  The array is a simple rank-N container of simple scalar values
 *  where 1 <= N <= 8.
 *
 *  The first N-1 ranks corresponds to a block of scalars.
 *  The last rank corresponds to the parallel work index.
 *
 *  No assumptions should be made as to the mapping, contiguity, or strides
 *  of the storage of these arrays.  The mapping will vary according to the
 *  underlying device.  The usage model is for algorithms to be parameterized
 *  with respect to the type of the mapped array and thus achieve portability
 *  across compute devices.
 */

template< typename ValueType >
class MDArrayView< ValueType , CudaMap > {
public:
  typedef CudaMap                        device_map_type ;
  typedef typename CudaMap::device_type  device_type ;
  typedef typename CudaMap::size_type    size_type ;
  typedef ValueType                      value_type ;

  /*------------------------------------------------------------------*/
  /** \brief  Query rank of the array */
  KOKKOS_DEVICE_AND_HOST_FUNCTION
  size_type rank() const { return m_data.m_rank ; }

  /** \brief  Query dimension of the given ordinate of the array */
  template < typename iType >
  KOKKOS_DEVICE_AND_HOST_FUNCTION
  size_type dimension( const iType & rank_ordinate ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_rank(rank_ordinate) );
      return m_data.m_dimension[rank_ordinate] ;
    }

  /*------------------------------------------------------------------*/
  /** \brief  Query value of a rank 8 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 , typename iTypeP >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 , const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_bounds(i0,i1,i2,i3,i4,i5,i6,iP) );

      return ( (value_type*) m_data.m_ptr_on_device )
             [ iP + m_data.m_dimension[7] *
               MULTI_INDEX_LEFT_7(i0,i1,i2,i3,i4,i5,i6,m_data.m_dimension) ];
    }

  /** \brief  Query value of a rank 7 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iTypeP >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_bounds(i0,i1,i2,i3,i4,i5,iP) );

      return ( (value_type*) m_data.m_ptr_on_device )
             [ iP + m_data.m_dimension[6] *
               MULTI_INDEX_LEFT_6(i0,i1,i2,i3,i4,i5,m_data.m_dimension) ];
    }

  /** \brief  Query value of a rank 6 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iTypeP >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_bounds(i0,i1,i2,i3,i4,iP) );

      return ( (value_type*) m_data.m_ptr_on_device )
             [ iP + m_data.m_dimension[5] *
               MULTI_INDEX_LEFT_5(i0,i1,i2,i3,i4,m_data.m_dimension) ];
    }

  /** \brief  Query value of a rank 5 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iTypeP >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_bounds(i0,i1,i2,i3,iP) );

      return ( (value_type*) m_data.m_ptr_on_device )
             [ iP + m_data.m_dimension[4] *
               MULTI_INDEX_LEFT_4(i0,i1,i2,i3,m_data.m_dimension) ];
    }

  /** \brief  Query value of a rank 4 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iTypeP >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_bounds(i0,i1,i2,iP) );

      return ( (value_type*) m_data.m_ptr_on_device )
             [ iP + m_data.m_dimension[3] *
               MULTI_INDEX_LEFT_3(i0,i1,i2,m_data.m_dimension) ];
    }

  /** \brief  Query value of a rank 3 array */
  template< typename iType0 , typename iType1 ,
            typename iTypeP >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_bounds(i0,i1,iP) );

      return ( (value_type*) m_data.m_ptr_on_device )
             [ iP + m_data.m_dimension[2] *
               MULTI_INDEX_LEFT_2(i0,i1,m_data.m_dimension) ];
    }

  /** \brief  Query value of a rank 2 array */
  template< typename iType0 , typename iTypeP >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_bounds(i0,iP) );

      return ( (value_type*) m_data.m_ptr_on_device )
             [ iP + m_data.m_dimension[1] *
               MULTI_INDEX_LEFT_1(i0,m_data.m_dimension) ];
    }

  /** \brief  Query value of a rank 1 array */
  template< typename iTypeP >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_bounds(iP) );

      return ( (value_type*) m_data.m_ptr_on_device )
             [ ( iP ) ];
    }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  inline
  MDArrayView() : m_data() {}

  /** \brief  Construct a view of the array */
  inline
  MDArrayView( const MDArrayView & rhs ) : m_data( rhs.m_data ) {}

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  inline
  MDArrayView & operator = ( const MDArrayView & rhs )
    { m_data.operator=( rhs.m_data ); return *this ; }
  
  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  inline
  ~MDArrayView() {}

  void clear_view()
    { m_data.clear_view(); }

  /*------------------------------------------------------------------*/

private:

  friend class CudaMap ;

  /** \brief  Separate view and data management details from
   *          multi-index mapping details.
   */
  MDArrayViewRawData<value_type,CudaDevice> m_data ;

  // Constructors that allocate the view on the device,
  // only callable by the CudaMap class.

  MDArrayView( size_type n0 , size_type n1 ,
               size_type n2 , size_type n3 ,
               size_type n4 , size_type n5 ,
               size_type n6 , size_type n7 ,
               const std::string & label )
    : m_data( n0, n1, n2, n3, n4, n5, n6, n7, label ) {}

  MDArrayView( size_type n0 , size_type n1 ,
               size_type n2 , size_type n3 ,
               size_type n4 , size_type n5 ,
               size_type n6 ,
               const std::string & label )
    : m_data( n0, n1, n2, n3, n4, n5, n6, label ) {}

  MDArrayView( size_type n0 , size_type n1 ,
               size_type n2 , size_type n3 ,
               size_type n4 , size_type n5 ,
               const std::string & label )
    : m_data( n0, n1, n2, n3, n4, n5, label ) {}

  MDArrayView( size_type n0 , size_type n1 ,
               size_type n2 , size_type n3 ,
               size_type n4 ,
               const std::string & label )
    : m_data( n0, n1, n2, n3, n4, label ) {}

  MDArrayView( size_type n0 , size_type n1 ,
               size_type n2 , size_type n3 ,
               const std::string & label )
    : m_data( n0, n1, n2, n3, label ) {}

  MDArrayView( size_type n0 , size_type n1 ,
               size_type n2 ,
               const std::string & label )
    : m_data( n0, n1, n2, label ) {}

  MDArrayView( size_type n0 , size_type n1 ,
               const std::string & label )
    : m_data( n0, n1, label ) {}

  MDArrayView( size_type n0 ,
               const std::string & label )
    : m_data( n0, label ) {}

  MDArrayView( const std::string & label )
    : m_data( label ) {}
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  Multidimensional array allocated on a compute device.
 *
 *  The array is a simple rank-N container of simple scalar values
 *  where 1 <= N <= 8.
 *
 *  No assumptions should be made as to the mapping, contiguity, or strides
 *  of the storage of these arrays.
 *
 *  Unmapped arrays are typically used for shared values such as
 *  the results of a reduction operation.
 */

template< typename ValueType >
class MDArrayView< ValueType , CudaDevice > {
public:
  typedef CudaDevice  device_map_type ;
  typedef CudaDevice  device_type ;
  typedef CudaDevice  size_type ;
  typedef ValueType   value_type ;

  /*------------------------------------------------------------------*/
  /** \brief  Query rank of the array */
  KOKKOS_DEVICE_AND_HOST_FUNCTION
  size_type rank() const { return m_data.m_rank ; }

  /** \brief  Query dimension of the given ordinate of the array */
  template < typename iType >
  KOKKOS_DEVICE_AND_HOST_FUNCTION
  size_type dimension( const iType & rank_ordinate ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_rank(rank_ordinate) );
      return m_data.m_dimension[rank_ordinate] ;
    }

  /*------------------------------------------------------------------*/
  /** \brief  Query value of a rank 7 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_bounds(i0,i1,i2,i3,i4,i5,i6) );

      return ( (value_type*) m_data.m_ptr_on_device )
             [ m_data.m_dimesion[7] *
               MULTI_INDEX_LEFT_7(i0,i1,i2,i3,i4,i5,i6,m_data.m_dimension) ];
    }

  /** \brief  Query value of a rank 6 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_bounds(i0,i1,i2,i3,i4,i5) );

      return ( (value_type*) m_data.m_ptr_on_device )
             [ m_data.m_dimension[6] *
               MULTI_INDEX_LEFT_6(i0,i1,i2,i3,i4,i5,m_data.m_dimension) ];
    }

  /** \brief  Query value of a rank 5 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_bounds(i0,i1,i2,i3,i4) );

      return ( (value_type*) m_data.m_ptr_on_device )
             [ m_data.m_dimension[5] *
               MULTI_INDEX_LEFT_5(i0,i1,i2,i3,i4,m_data.m_dimension) ];
    }

  /** \brief  Query value of a rank 4 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_bounds(i0,i1,i2,i3) );

      return ( (value_type*) m_data.m_ptr_on_device )
             [ m_data.m_dimension[4] *
               MULTI_INDEX_LEFT_4(i0,i1,i2,i3,m_data.m_dimension) ];
    }

  /** \brief  Query value of a rank 3 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_bounds(i0,i1,i2) );

      return ( (value_type*) m_data.m_ptr_on_device )
             [ m_data.m_dimension[3] *
               MULTI_INDEX_LEFT_3(i0,i1,i2,m_data.m_dimension) ];
    }

  /** \brief  Query value of a rank 2 array */
  template< typename iType0 , typename iType1 >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_bounds(i0,i1) );

      return ( (value_type*) m_data.m_ptr_on_device )
             [ m_data.m_dimension[2] *
               MULTI_INDEX_LEFT_2(i0,i1,m_data.m_dimension) ];
    }

  /** \brief  Query value of a rank 1 array */
  template< typename iType0 >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 ) const
    {
      KOKKOS_MDARRAY_CHECK( m_data.require_in_bounds(i0) );

      return ( (value_type*) m_data.m_ptr_on_device )
             [ m_data.m_dimension[1] * ( i0 ) ];
    }

  /** \brief  Query value of a rank 0 array */
  template< typename iType0 >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()() const
    {
      return *( (value_type*) m_data.m_ptr_on_device );
    }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  inline
  MDArrayView() : m_data() {}

  /** \brief  Construct a view of the array */
  inline
  MDArrayView( const MDArrayView & rhs ) : m_data( rhs.m_data ) {}

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  inline
  MDArrayView & operator = ( const MDArrayView & rhs )
    { m_data.operator=( rhs.m_data ); return *this ; }
  
  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  inline
  ~MDArrayView() {}

  void clear_view()
    { m_data.clear_view(); }

  /*------------------------------------------------------------------*/
  /** \brief  Wrap MDArrayView around CUDA memory on device. */

  __device__
  void assign_on_device( const MDArrayView & rhs )
    {
      m_data.m_ptr_on_device = rhs.m_data.m_ptr_on_device ;
      m_data.m_rank          = rhs.m_data.m_rank ;
      m_data.m_dimension[0]  = rhs.m_data.m_dimension[0] ;
      m_data.m_dimension[1]  = rhs.m_data.m_dimension[1] ;
      m_data.m_dimension[2]  = rhs.m_data.m_dimension[2] ;
      m_data.m_dimension[3]  = rhs.m_data.m_dimension[3] ;
      m_data.m_dimension[4]  = rhs.m_data.m_dimension[4] ;
      m_data.m_dimension[5]  = rhs.m_data.m_dimension[5] ;
      m_data.m_dimension[6]  = rhs.m_data.m_dimension[6] ;
      m_data.m_dimension[7]  = rhs.m_data.m_dimension[7] ;
    }

  __device__
  void assign_on_device( value_type * ptr_on_device , size_type thread_stride )
    {
      m_data.m_ptr_on_device              = ptr_on_device ;
      m_data.m_dimension[ m_data.m_rank ] = thread_stride ;
    }

  __device__
  void assign_on_device( value_type * ptr_on_device )
    { m_data.m_ptr_on_device = ptr_on_device ; }

  __device__
  value_type * address_on_device() const { return m_data.m_ptr_on_device ; }

  /*------------------------------------------------------------------*/

private:

  friend class CudaDevice ;

  /** \brief  Separate view and data management details from
   *          multi-index mapping details.
   */
  MDArrayViewRawData<value_type,CudaDevice> m_data ;

  // Constructors that allocate the view on the device

  MDArrayView( size_type n0 , size_type n1 ,
               size_type n2 , size_type n3 ,
               size_type n4 , size_type n5 ,
               size_type n6 ,
               const std::string & label )
    : m_data( n0, n1, n2, n3, n4, n5, n6, label )
    { m_data.m_dimension[7] = 1 ; }

  MDArrayView( size_type n0 , size_type n1 ,
               size_type n2 , size_type n3 ,
               size_type n4 , size_type n5 ,
               const std::string & label )
    : m_data( n0, n1, n2, n3, n4, n5, label )
    { m_data.m_dimension[6] = 1 ; }

  MDArrayView( size_type n0 , size_type n1 ,
               size_type n2 , size_type n3 ,
               size_type n4 ,
               const std::string & label )
    : m_data( n0, n1, n2, n3, n4, label )
    { m_data.m_dimension[5] = 1 ; }

  MDArrayView( size_type n0 , size_type n1 ,
               size_type n2 , size_type n3 ,
               const std::string & label )
    : m_data( n0, n1, n2, n3, label )
    { m_data.m_dimension[4] = 1 ; }

  MDArrayView( size_type n0 , size_type n1 ,
               size_type n2 ,
               const std::string & label )
    : m_data( n0, n1, n2, label )
    { m_data.m_dimension[3] = 1 ; }

  MDArrayView( size_type n0 , size_type n1 ,
               const std::string & label )
    : m_data( n0, n1, label )
    { m_data.m_dimension[2] = 1 ; }

  MDArrayView( size_type n0 ,
               const std::string & label )
    : m_data( n0, label )
    { m_data.m_dimension[1] = 1 ; }

  MDArrayView( const std::string & label )
    : m_data( label )
    { m_data.m_dimension[0] = 1 ; }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Kokkos

#endif /* KOKKOS_CUDAMDARRAYVIEW_HPP */


