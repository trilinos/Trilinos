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
template< typename ValueType , class DeviceMapType >
class MDArray {
public:
  typedef ValueType                             value_type ;
  typedef DeviceMapType                         device_map_type ;
  typedef typename device_map_type::size_type   size_type ;

  /*------------------------------------------------------------------*/
  /** \brief  Query dimension of the given ordinate of the array */
  template < typename iType >
  size_type dimension( const iType & ordinal ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Query value of a rank 8 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 , typename iTypeP >
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 , const iTypeP & iP ) const ;

  /** \brief  Query value of a rank 7 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iTypeP >
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iTypeP & iP ) const ;

  /** \brief  Query value of a rank 6 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iTypeP >
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iTypeP & iP ) const ;

  /** \brief  Query value of a rank 5 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iTypeP >
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iTypeP & iP ) const ;

  /** \brief  Query value of a rank 4 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iTypeP >
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iType2 & i2 , const iTypeP & iP ) const ;

  /** \brief  Query value of a rank 3 array */
  template< typename iType0 , typename iType1 ,
            typename iTypeP >
  value_type & operator()( const iType0 & i0 , const iType1 & i1 ,
                           const iTypeP & iP ) const ;

  /** \brief  Query value of a rank 2 array */
  template< typename iType0 , typename iTypeP >
  value_type & operator()( const iType0 & i0 , const iTypeP & iP ) const ;

  /** \brief  Query value of a rank 1 array */
  template< typename iTypeP >
  value_type & operator()( const iTypeP & iP ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  MDArray();

  /** \brief  Construct a view of the array */
  MDArray( const MDArray & rhs );

  /** \brief  Assign a view of the array, the old view is discarded. */
  MDArray & operator = ( const MDArray & rhs );
  
  /**  Destroy this view of the array, the memory is not deallocated. */
  ~MDArray();

  /*------------------------------------------------------------------*/
  /** \brief  Allocate a rank 8 array on device using the device array map */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 >
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                 const iType2 & n2 , const iType3 & n3 ,
                 const iType4 & n4 , const iType5 & n5 ,
                 const iType6 & n6 , device_map_type & map );

  /** \brief  Allocate a rank 7 array on device using the device array map */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                 const iType2 & n2 , const iType3 & n3 ,
                 const iType4 & n4 , const iType5 & n5 ,
                 device_map_type & map );

  /** \brief  Allocate a rank 6 array on device using the device array map */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 >
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                 const iType2 & n2 , const iType3 & n3 ,
                 const iType4 & n4 , device_map_type & map );

  /** \brief  Allocate a rank 5 array on device using the device array map */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 >
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                 const iType2 & n2 , const iType3 & n3 ,
                 device_map_type & map );

  /** \brief  Allocate a rank 4 array on device using the device array map */
  template< typename iType0 , typename iType1 ,
            typename iType2 >
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                 const iType2 & n2 , device_map_type & map );

  /** \brief  Allocate a rank 3 array on device using the device array map */
  template< typename iType0 , typename iType1 >
  void allocate( const iType0 & n0 , const iType1 & n1 ,
                 device_map_type & map );

  /** \brief  Allocate a rank 2 array on device using the device array map */
  template< typename iType0 >
  void allocate( const iType0 & n0 , device_map_type & map );

  /** \brief  Allocate a rank 1 array on device using the device array map */
  void allocate( device_map_type & map );

  /** \brief  Deallocate an array, all views are set to NULL */
  void deallocate();

};

} // namespace Kokkos

#endif /* KOKKOS_MAPPEDARRAY_HPP */


