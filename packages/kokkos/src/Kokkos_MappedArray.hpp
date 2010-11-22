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

#ifndef KOKKOS_MAPPEDARRAY_HPP
#define KOKKOS_MAPPEDARRAY_HPP

namespace Kokkos {

/** \brief  Multidimensional array allocated and mapped
 *          onto a compute device.
 */
template< typename ValueType , class DeviceMapType >
class MDArray {
private:
  typedef typename DeviceMapType::mapped_array_type storage_type ;
public:
  typedef ValueType                             value_type ;
  typedef DeviceMapType                         device_map_type ;
  typedef typename device_map_type::device_type device_type ;
  typedef typename device_map_type::size_type   size_type ;

  /*------------------------------------------------------------------*/

  template < typename iType >
  size_type dimension( const iType & ordinal ) const
    {
      KOKKOS_MDARRAY_CHECK( require_ordinal_in_bounds( m_rank , ordinal ) );
      return ordinal < m_rank - 1
             ? m_dimension[ordinal] : m_storage.work_size();
    }

  /*------------------------------------------------------------------*/

  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 , typename iTypeP >
  inline
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 , const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK(
        require_multi_index_in_bounds( m_rank, m_dimension,
                                       m_storage.parallel_work_size(),
                                       i0,i1,i2,i3,i4,i5,i6,iP) );

      // Map all but the parallel work index to a single offset
      const size_type offset = ( i0 + m_dimension[0] * 
                               ( i1 + m_dimension[1] *
                               ( i2 + m_dimension[2] *
                               ( i3 + m_dimension[3] *
                               ( i4 + m_dimension[4] *
                               ( i5 + m_dimension[5] *
                               ( i6 )))))));

      return m_storage.template value< value_type >( offset , iP );
    }

  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iTypeP >
  inline
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK(
        require_multi_index_in_bounds( m_rank, m_dimension,
                                       m_storage.parallel_work_size(),
                                       i0,i1,i2,i3,i4,i5,iP) );

      // Map all but the parallel work index to a single offset
      const size_type offset = ( i0 + m_dimension[0] * 
                               ( i1 + m_dimension[1] *
                               ( i2 + m_dimension[2] *
                               ( i3 + m_dimension[3] *
                               ( i4 + m_dimension[4] *
                               ( i5 ))))));

      return m_storage.template value< value_type >( offset , iP );
    }

  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iTypeP >
  inline
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK(
        require_multi_index_in_bounds( m_rank, m_dimension,
                                       m_storage.parallel_work_size(),
                                       i0,i1,i2,i3,i4,iP) );

      // Map all but the parallel work index to a single offset
      const size_type offset = ( i0 + m_dimension[0] * 
                               ( i1 + m_dimension[1] *
                               ( i2 + m_dimension[2] *
                               ( i3 + m_dimension[3] *
                               ( i4 )))));

      return m_storage.template value< value_type >( offset , iP );
    }

  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iTypeP >
  inline
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK(
        require_multi_index_in_bounds( m_rank, m_dimension,
                                       m_storage.parallel_work_size(),
                                       i0,i1,i2,i3,iP) );

      // Map all but the parallel work index to a single offset
      const size_type offset = ( i0 + m_dimension[0] * 
                               ( i1 + m_dimension[1] *
                               ( i2 + m_dimension[2] *
                               ( i3 ))));

      return m_storage.template value< value_type >( offset , iP );
    }

  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iTypeP >
  inline
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK(
        require_multi_index_in_bounds( m_rank, m_dimension,
                                       m_storage.parallel_work_size(),
                                       i0,i1,i2,iP) );

      // Map all but the parallel work index to a single offset
      const size_type offset = ( i0 + m_dimension[0] * 
                               ( i1 + m_dimension[1] *
                               ( i2 )));

      return m_storage.template value< value_type >( offset , iP );
    }

  template< typename iType0 , typename iType1 ,
            typename iTypeP >
  inline
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK(
        require_multi_index_in_bounds( m_rank, m_dimension,
                                       m_storage.parallel_work_size(),
                                       i0,i1,iP) );

      // Map all but the parallel work index to a single offset
      const size_type offset = ( i0 + m_dimension[0] * 
                               ( i1 ));

      return m_storage.template value< value_type >( offset , iP );
    }

  template< typename iType0 , typename iTypeP >
  inline
  value_type & operator()( const iType0 & i0 , const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK(
        require_multi_index_in_bounds( m_rank, m_dimension,
                                       m_storage.parallel_work_size(),
                                       i0,iP) );
      return m_storage.template value< value_type >( i0 , iP );
    }

  template< typename iTypeP >
  inline
  value_type & operator()( const iTypeP & iP ) const
    {
      KOKKOS_MDARRAY_CHECK(
        require_multi_index_in_bounds( m_rank, m_dimension,
                                       m_storage.parallel_work_size(),
                                       iP) );
      return m_storage.template value< value_type >( 0 , iP );
    }

  /*------------------------------------------------------------------*/

  MDArray()
    : m_storage()
    , m_rank(0) 
    {
      m_dimension[0] = 0 ;
      m_dimension[1] = 0 ;
      m_dimension[2] = 0 ;
      m_dimension[3] = 0 ;
      m_dimension[4] = 0 ;
      m_dimension[5] = 0 ;
      m_dimension[6] = 0 ;
    }

  /** \brief  Construct a view of the array */
  MDArray( const MDArray & rhs )
    : m_storage( rhs.m_storage )
    , m_rank(    rhs.m_rank )
    {
      m_dimension[0] = rhs.m_dimension[0] ;
      m_dimension[1] = rhs.m_dimension[1] ;
      m_dimension[2] = rhs.m_dimension[2] ;
      m_dimension[3] = rhs.m_dimension[3] ;
      m_dimension[4] = rhs.m_dimension[4] ;
      m_dimension[5] = rhs.m_dimension[5] ;
      m_dimension[6] = rhs.m_dimension[6] ;
    }

  /** \brief  Assign a view of the array */
  MDArray & operator = ( const MDArray & rhs )
    {
      m_storage      = rhs.m_storage ;
      m_rank         = rhs.m_rank ;
      m_dimension[0] = rhs.m_dimension[0] ;
      m_dimension[1] = rhs.m_dimension[1] ;
      m_dimension[2] = rhs.m_dimension[2] ;
      m_dimension[3] = rhs.m_dimension[3] ;
      m_dimension[4] = rhs.m_dimension[4] ;
      m_dimension[5] = rhs.m_dimension[5] ;
      m_dimension[6] = rhs.m_dimension[6] ;
      return *this ;
    }
  
  /**  Destroy this view of the array */
  inline
  ~MDArray() {}

  /*------------------------------------------------------------------*/
  /** \brief  Allocate array on device using the device array map */
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
      const size_type n =
       ( m_dimension[0] = n0 ) *
       ( m_dimension[1] = n1 ) *
       ( m_dimension[2] = n2 ) *
       ( m_dimension[3] = n3 ) *
       ( m_dimension[4] = n4 ) *
       ( m_dimension[5] = n5 ) *
       ( m_dimension[6] = n6 );
      m_rank = 8 ;
      map.allocate( m_storage , n , sizeof(ValueType) );
    }

  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  inline
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                 const iType2 & n2 , const iType3 & n3 ,
                 const iType4 & n4 , const iType5 & n5 ,
                 device_map_type & map )
    {
      const size_type n =
       ( m_dimension[0] = n0 ) *
       ( m_dimension[1] = n1 ) *
       ( m_dimension[2] = n2 ) *
       ( m_dimension[3] = n3 ) *
       ( m_dimension[4] = n4 ) *
       ( m_dimension[5] = n5 );
      m_rank = 7 ;
      map.allocate( m_storage , n , sizeof(ValueType) );
    }

  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 >
  inline
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                 const iType2 & n2 , const iType3 & n3 ,
                 const iType4 & n4 , device_map_type & map )
    {
      const size_type n =
       ( m_dimension[0] = n0 ) *
       ( m_dimension[1] = n1 ) *
       ( m_dimension[2] = n2 ) *
       ( m_dimension[3] = n3 ) *
       ( m_dimension[4] = n4 );
      m_rank = 6 ;
      map.allocate( m_storage , n , sizeof(ValueType) );
    }

  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 >
  inline
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                 const iType2 & n2 , const iType3 & n3 ,
                 device_map_type & map )
    {
      const size_type n =
       ( m_dimension[0] = n0 ) *
       ( m_dimension[1] = n1 ) *
       ( m_dimension[2] = n2 ) *
       ( m_dimension[3] = n3 );
      m_rank = 5 ;
      map.allocate( m_storage , n , sizeof(ValueType) );
    }

  template< typename iType0 , typename iType1 ,
            typename iType2 >
  inline
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                 const iType2 & n2 , device_map_type & map )
    {
      const size_type n =
       ( m_dimension[0] = n0 ) *
       ( m_dimension[1] = n1 ) *
       ( m_dimension[2] = n2 );
      m_rank = 4 ;
      map.allocate( m_storage , n , sizeof(ValueType) );
    }

  template< typename iType0 , typename iType1 >
  inline
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                    device_map_type & map )
    {
      const size_type n =
       ( m_dimension[0] = n0 ) *
       ( m_dimension[1] = n1 );
      m_rank = 3 ;
      map.allocate( m_storage , n , sizeof(ValueType) );
    }

  template< typename iType0 >
  inline
  void allocate( const iType0 & n0 , device_map_type & map )
    {
      const size_type n = ( m_dimension[0] = n0 );
      m_rank = 2 ;
      map.allocate( m_storage , n , sizeof(ValueType) );
    }

  inline
  void allocate( device_map_type & map )
    {
      m_rank = 1 ;
      map.allocate( m_storage , 1 , sizeof(ValueType) );
    }

  inline
  void deallocate()
    {
      device_map_type * const map = m_storage.map();
      if ( map ) { map->deallocate( m_storage ); }
    }

  /*------------------------------------------------------------------*/

private:
  storage_type m_storage ;
  size_type    m_rank ;
  size_type    m_dimension[7] ;
};


} // namespace Kokkos

#endif /* KOKKOS_MAPPEDARRAY_HPP */


