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

#ifndef KOKKOS_MDARRAY_HPP
#define KOKKOS_MDARRAY_HPP

#include <cstddef>
#include <string>
#include <impl/KokkosArray_forward.hpp>
#include <impl/KokkosArray_StaticAssert.hpp>
#include <impl/KokkosArray_ArrayBounds.hpp>

namespace KokkosArray {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  Multidimensional array allocated and mapped
 *          onto the memory space of a compute device.
 *
 *  The array is a simple rank-N container of simple scalar values
 *  where 1 <= N <= 8.
 *
 *  The first rank is the parallel work index.
 *
 *  No assumptions should be made as to the mapping, contiguity, or strides
 *  of the storage of these arrays.  The mapping will vary according to the
 *  underlying device.  The usage model is for algorithms to be parameterized
 *  with respect to the type of the mapped array and thus achieve portability
 *  across compute devices.
 *
 *  Several APIs for MDArray creation functions are available.
 *  The "labeled" group creates MDArray's with string labels that
 *  will appear in error messages:
 *  \code
 *  create_mdarray< MDArray<...> >( label , nP , ... );
 *  \endcode
 *  The "unlabeled" group creates MDArray's with NULL string labels:
 *  \code
 *  create_mdarray< MDArray<...> >( nP , ... );
 *  \endcode
 */
template< typename ValueType , class DeviceType >
class MDArray {
public:
  typedef ValueType                                       value_type ;
  typedef DeviceType                                      device_type ;
  typedef typename DeviceType::size_type                  size_type ;
  typedef typename DeviceType::template IndexMap<>::type  index_map ;

  typedef MDArray< value_type ,
                   typename HostMapped<device_type>::type > HostMirror ;

  /*------------------------------------------------------------------*/
  /** \brief  Query rank of the array */
  size_type rank() const ;

  /** \brief  Query dimension of the given ordinate of the array */
  template < typename iType >
  size_type dimension( const iType & rank_ordinate ) const ;

  /** \brief  Query all dimensions */
  template< typename iType >
  void dimensions( iType * const dims ) const ;

  /** \brief  Query total number of members */
  size_type size() const ;

  /*------------------------------------------------------------------*/
  /** \brief  Query value of a rank 8 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 , typename iType7 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 , const iType7 & i7 ) const ;

  /** \brief  Query value of a rank 7 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 ) const ;

  /** \brief  Query value of a rank 6 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ) const ;

  /** \brief  Query value of a rank 5 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 ) const ;

  /** \brief  Query value of a rank 4 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ) const ;

  /** \brief  Query value of a rank 3 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 ) const ;

  /** \brief  Query value of a rank 2 array */
  template< typename iTypeP , typename iType1 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ) const ;

  /** \brief  Query value of a rank 1 array */
  template< typename iTypeP >
  value_type & operator()( const iTypeP & iP ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Memory is contiguous, OK to return pointer */
  value_type * ptr_on_device() const ;

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  MDArray();

  /** \brief  Construct another view of the 'rhs' array */
  MDArray( const MDArray & rhs );

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  MDArray & operator = ( const MDArray & rhs );

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  ~MDArray();

  /*------------------------------------------------------------------*/
  /** \brief  Query if non-NULL view */
  operator bool () const ;

  bool operator == ( const MDArray & rhs ) const ;

  bool operator != ( const MDArray & rhs ) const ;
};

//----------------------------------------------------------------------------

/** \brief  Create an MDArray with a label */
template< typename MDArrayType >
inline
MDArray< typename MDArrayType::value_type ,
         typename MDArrayType::device_type >
create_mdarray( const std::string & label ,
                size_t nP ,     size_t n1 = 0 ,
                size_t n2 = 0 , size_t n3 = 0 ,
                size_t n4 = 0 , size_t n5 = 0 ,
                size_t n6 = 0 , size_t n7 = 0 )
{
  return Impl::Factory< MDArrayType , void >
           ::create( label, nP, n1, n2, n3, n4, n5, n6, n7 );
}

/** \brief  Create an MDArray without a label */
template< typename MDArrayType >
inline
MDArray< typename MDArrayType::value_type ,
         typename MDArrayType::device_type >
create_mdarray( size_t nP ,     size_t n1 = 0 ,
                size_t n2 = 0 , size_t n3 = 0 ,
                size_t n4 = 0 , size_t n5 = 0 ,
                size_t n6 = 0 , size_t n7 = 0 )
{
  return Impl::Factory< MDArrayType , void >
           ::create( std::string(), nP, n1, n2, n3, n4, n5, n6, n7 );
}

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceDst , class DeviceSrc >
inline
void deep_copy( const MDArray<ValueType,DeviceDst> & dst ,
                const MDArray<ValueType,DeviceSrc> & src )
{
  typedef MDArray<ValueType,DeviceDst>  dst_type ;
  typedef MDArray<ValueType,DeviceSrc>  src_type ;

  if ( dst.operator!=( src ) ) {
    Impl::mdarray_require_equal_dimension( dst , src );

    Impl::Factory< dst_type, src_type >::deep_copy( dst , src );
  }
}

//----------------------------------------------------------------------------

} // namespace KokkosArray

#include <impl/KokkosArray_MDArray_factory.hpp>

#endif /* KOKKOS_MDARRAY_HPP */


