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

#ifndef KOKKOS_CUDAMAPPEDARRAY_HPP
#define KOKKOS_CUDAMAPPEDARRAY_HPP

#include <Kokkos_BaseMappedArray.hpp>

#define KOKKOS_DEVICE_FUNCTION inline __device__ __host__

namespace Kokkos {

template< typename ValueType , class DeviceType > class MDArray ;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class CudaMap : public BaseMapInterface {
public:

  ~CudaMap();

  CudaMap( size_type parallel_work_count );

  size_type parallel_work_count() const { return m_parallel_work_count ; }

  void deallocate( BaseMappedArray & );

  void allocate( BaseMappedArray & ,
                 size_type sizeof_value ,
                 size_type rank ,
                 const size_type * const dimension );

private:

  std::list< BaseMappedArray > m_allocated_arrays ;
  size_type                    m_parallel_work_count ;

  CudaMap();
  CudaMap( const CudaMap & );
  CudaMap & operator = ( const CudaMap & );
};

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
class MDArray< ValueType , CudaMap > {
public:
  typedef ValueType                   value_type ;
  typedef CudaMap                     device_map_type ;
  typedef BaseMappedArray::size_type  size_type ;

  /*------------------------------------------------------------------*/
  /** \brief  Query dimension of the given ordinate of the array */
  template < typename iType >
  KOKKOS_DEVICE_FUNCTION
  size_type dimension( const iType & ordinal ) const
    {
      KOKKOS_MDARRAY_CHECK( m_base.require_in_rank(ordinal) );
      return m_base.m_dimension[ordinal] ;
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
      KOKKOS_MDARRAY_CHECK( m_base.require_in_bounds(i0,i1,i2,i3,i4,i5,i6,iP) );

      return ( (value_type*) m_base.m_ptr_on_device )
             [ ( iP + m_base.m_dimension[7] *
               ( i0 + m_base.m_dimension[0] * ( i1 + m_base.m_dimension[1] *
               ( i2 + m_base.m_dimension[2] * ( i3 + m_base.m_dimension[3] *
               ( i4 + m_base.m_dimension[4] * ( i5 + m_base.m_dimension[5] *
               ( i6 )))))))) ];
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
      KOKKOS_MDARRAY_CHECK( m_base.require_in_bounds(i0,i1,i2,i3,i4,i5,iP) );

      return ( (value_type*) m_base.m_ptr_on_device )
             [ ( iP + m_base.m_dimension[6] *
               ( i0 + m_base.m_dimension[0] * ( i1 + m_base.m_dimension[1] *
               ( i2 + m_base.m_dimension[2] * ( i3 + m_base.m_dimension[3] *
               ( i4 + m_base.m_dimension[4] * ( i5 ))))))) ];
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
      KOKKOS_MDARRAY_CHECK( m_base.require_in_bounds(i0,i1,i2,i3,i4,iP) );

      return ( (value_type*) m_base.m_ptr_on_device )
             [ ( iP + m_base.m_dimension[5] *
               ( i0 + m_base.m_dimension[0] * ( i1 + m_base.m_dimension[1] *
               ( i2 + m_base.m_dimension[2] * ( i3 + m_base.m_dimension[3] *
               ( i4 )))))) ];
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
      KOKKOS_MDARRAY_CHECK( m_base.require_in_bounds(i0,i1,i2,i3,iP) );

      return ( (value_type*) m_base.m_ptr_on_device )
             [ ( iP + m_base.m_dimension[4] *
               ( i0 + m_base.m_dimension[0] * ( i1 + m_base.m_dimension[1] *
               ( i2 + m_base.m_dimension[2] * ( i3 ))))) ];
    }

  /** \brief  Query value of a rank 4 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iTypeP >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK( m_base.require_in_bounds(i0,i1,i2,iP) );

      return ( (value_type*) m_base.m_ptr_on_device )
             [ ( iP + m_base.m_dimension[3] *
               ( i0 + m_base.m_dimension[0] * ( i1 + m_base.m_dimension[1] *
               ( i2 )))) ];
    }

  /** \brief  Query value of a rank 3 array */
  template< typename iType0 , typename iType1 ,
            typename iTypeP >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK( m_base.require_in_bounds(i0,i1,iP) );

      return ( (value_type*) m_base.m_ptr_on_device )
             [ ( iP + m_base.m_dimension[2] *
               ( i0 + m_base.m_dimension[0] * ( i1 ))) ];
    }

  /** \brief  Query value of a rank 2 array */
  template< typename iType0 , typename iTypeP >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iType0 & i0 , const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK( m_base.require_in_bounds(i0,iP) );

      return ( (value_type*) m_base.m_ptr_on_device )
             [ ( iP + m_base.m_dimension[1] * ( i0 )) ];
    }

  /** \brief  Query value of a rank 1 array */
  template< typename iTypeP >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK( m_base.require_in_bounds(iP) );

      return ( (value_type*) m_base.m_ptr_on_device )
             [ ( iP ) ];
    }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  inline
  MDArray() : m_base() {}

  /** \brief  Construct a view of the array */
  inline
  MDArray( const MDArray & rhs ) : m_base( rhs.m_base ) {}

  /** \brief  Assign a view of the array, the old view is discarded. */
  inline
  MDArray & operator = ( const MDArray & rhs )
    { m_base = rhs.m_base ; return *this ; }
  
  /**  Destroy this view of the array, the memory is not deallocated. */
  inline
  ~MDArray() {}

  /*------------------------------------------------------------------*/
  /** \brief  Allocate a rank 8 array on device using the device array map */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 >
  inline
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                 const iType2 & n2 , const iType3 & n3 ,
                 const iType4 & n4 , const iType5 & n5 ,
                 const iType6 & n6 , device_map_type & map )
    {
      size_type dimension[7] ;
      dimension[0] = n0 ;
      dimension[1] = n1 ;
      dimension[2] = n2 ;
      dimension[3] = n3 ;
      dimension[4] = n4 ;
      dimension[5] = n5 ;
      dimension[6] = n6 ;
      map.allocate( m_base , sizeof(value_type) , 8 , dimension );
    }

  /** \brief  Allocate a rank 7 array on device using the device array map */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  inline
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                 const iType2 & n2 , const iType3 & n3 ,
                 const iType4 & n4 , const iType5 & n5 ,
                 device_map_type & map )
    {
      size_type dimension[6] ;
      dimension[0] = n0 ;
      dimension[1] = n1 ;
      dimension[2] = n2 ;
      dimension[3] = n3 ;
      dimension[4] = n4 ;
      dimension[5] = n5 ;
      map.allocate( *this , sizeof(value_type) , 7 , dimension );
    }

  /** \brief  Allocate a rank 6 array on device using the device array map */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 >
  inline
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                 const iType2 & n2 , const iType3 & n3 ,
                 const iType4 & n4 , device_map_type & map )
    {
      size_type dimension[5] ;
      dimension[0] = n0 ;
      dimension[1] = n1 ;
      dimension[2] = n2 ;
      dimension[3] = n3 ;
      dimension[4] = n4 ;
      map.allocate( m_base , sizeof(value_type) , 6 , dimension );
    }

  /** \brief  Allocate a rank 5 array on device using the device array map */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 >
  inline
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                 const iType2 & n2 , const iType3 & n3 ,
                 device_map_type & map )
    {
      size_type dimension[4] ;
      dimension[0] = n0 ;
      dimension[1] = n1 ;
      dimension[2] = n2 ;
      dimension[3] = n3 ;
      map.allocate( m_base , sizeof(value_type) , 5 , dimension );
    }

  /** \brief  Allocate a rank 4 array on device using the device array map */
  template< typename iType0 , typename iType1 ,
            typename iType2 >
  inline
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                 const iType2 & n2 , device_map_type & map )
    {
      size_type dimension[3] ;
      dimension[0] = n0 ;
      dimension[1] = n1 ;
      dimension[2] = n2 ;
      map.allocate( m_base , sizeof(value_type) , 4 , dimension );
    }

  /** \brief  Allocate a rank 3 array on device using the device array map */
  template< typename iType0 , typename iType1 >
  inline
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                    device_map_type & map )
    {
      size_type dimension[2] ;
      dimension[0] = n0 ;
      dimension[1] = n1 ;
      map.allocate( m_base , sizeof(value_type) , 3 , dimension );
    }

  /** \brief  Allocate a rank 2 array on device using the device array map */
  template< typename iType0 >
  inline
  void allocate( const iType0 & n0 , device_map_type & map )
    {
      size_type dimension[1] ;
      dimension[0] = n0 ;
      map.allocate( m_base , sizeof(value_type) , 2 , dimension );
    }

  /** \brief  Allocate a rank 1 array on device using the device array map */
  inline
  void allocate( device_map_type & map )
    {
      map.allocate( m_base, sizeof(value_type) , 1 , (size_type *) NULL );
    }

  /** \brief  Deallocate an array, all views are set to NULL */
  inline
  void deallocate()
    {
      device_map_type * const map = m_base.map<device_map_type>();
      if ( map ) { map->deallocate( m_base ); }
    }

private:
  /** \brief  A single data member to manage multiple views. */
  BaseMappedArray m_base ;
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


} // namespace Kokkos

#endif /* KOKKOS_CUDAMAPPEDARRAY_HPP */


