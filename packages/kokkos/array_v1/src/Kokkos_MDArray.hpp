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
#include <impl/Kokkos_forward.hpp>
#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_ArrayBounds.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------

enum { MDArrayMaxRank = 8 };

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
 *  will appear in error messages.
 *
 *  create_labeled_mdarray< ValueType , DeviceType >( label , nP , ... );
 *  create_labeled_mdarray< MDArray<...> >( label , nP , ... );
 *
 *  The "unlabeled" group creates MDArray's with NULL string labels.
 *
 *  create_mdarray< ValueType , DeviceType >( nP , ... );
 *  create_mdarray< MDArray<...> >( nP , ... );
 */
template< typename ValueType , class DeviceType >
class MDArray {
public:
  typedef ValueType                          value_type ;
  typedef DeviceType                         device_type ;
  typedef typename DeviceType::mdarray_map   mdarray_map ;
  typedef typename DeviceType::size_type     size_type ;

  typedef MDArray< value_type , void /* HostMapped< mdarray_map > */ > HostMirror ;

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

template< typename ValueType , class DeviceType >
MDArray< ValueType , DeviceType >
create_labeled_mdarray( const std::string & label ,
                        size_t nP , size_t n1 = 0 ,
                        size_t n2 = 0 , size_t n3 = 0 ,
                        size_t n4 = 0 , size_t n5 = 0 ,
                        size_t n6 = 0 , size_t n7 = 0 );

template< typename MDArrayType >
inline
MDArray< typename MDArrayType::value_type ,
         typename MDArrayType::device_type >
create_labeled_mdarray( const std::string & label ,
                        size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                                    size_t n3 = 0 , size_t n4 = 0 ,
                                    size_t n5 = 0 , size_t n6 = 0 ,
                                    size_t n7 = 0 );

template< typename ValueType , class DeviceType >
inline
MDArray< ValueType , DeviceType >
create_mdarray( size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                            size_t n3 = 0 , size_t n4 = 0 ,
                            size_t n5 = 0 , size_t n6 = 0 ,
                            size_t n7 = 0 );

template< typename MDArrayType >
inline
MDArray< typename MDArrayType::value_type ,
         typename MDArrayType::device_type >
create_mdarray( size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                            size_t n3 = 0 , size_t n4 = 0 ,
                            size_t n5 = 0 , size_t n6 = 0 ,
                            size_t n7 = 0 );

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
typename MDArray< ValueType , DeviceType >::HostMirror
create_mirror( const MDArray< ValueType , DeviceType > & );

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceDst , class DeviceSrc >
inline
void deep_copy( const MDArray<ValueType,DeviceDst> & dst ,
                const MDArray<ValueType,DeviceSrc> & src );

//----------------------------------------------------------------------------
/** \brief  This is THE creation function.
 *          All other versions call this function.
 */
template< typename ValueType , class DeviceType >
inline
MDArray< ValueType , DeviceType >
create_labeled_mdarray( const std::string & label ,
                        size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                        size_t n4 , size_t n5 , size_t n6 , size_t n7 )
{
  return MDArray< ValueType , DeviceType >
           ( label , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

template< typename MDArrayType >
inline
MDArray< typename MDArrayType::value_type , typename MDArrayType::device_type >
create_labeled_mdarray( const std::string & label ,
                        size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                        size_t n4 , size_t n5 , size_t n6 , size_t n7 )
{
  return create_labeled_mdarray< typename MDArrayType::value_type ,
                                 typename MDArrayType::device_type >
           ( label , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

template< typename ValueType , class DeviceType >
inline
MDArray< ValueType , DeviceType >
create_mdarray( size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                size_t n4 , size_t n5 , size_t n6 , size_t n7 )
{
  return create_labeled_mdarray< ValueType , DeviceType >
           ( std::string() , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

template< typename MDArrayType >
inline
MDArray< typename MDArrayType::value_type , typename MDArrayType::device_type >
create_mdarray( size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                size_t n4 , size_t n5 , size_t n6 , size_t n7 )
{
  return create_labeled_mdarray< typename MDArrayType::value_type ,
                                 typename MDArrayType::device_type >
           ( std::string() , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Impl {

template< unsigned N >
class Rank { public: enum { value = N }; };

template< typename ValueType , class MDArrayMap >
class MDArrayHostMirror ;

template< typename ValueType , class Device >
class CreateMirror< MDArray< ValueType , Device > , true /* view */ >
{
public:
  typedef  MDArray< ValueType , Device >            View ;
  typedef  typename MDArray< ValueType , Device >::HostMirror  HostMirror ;

  static
  HostMirror create( const View & v ) { return HostMirror( v ); }
};

template< typename ValueType , class Device >
class CreateMirror< MDArray< ValueType , Device > , false /* copy */ >
{
public:
  typedef  MDArray< ValueType , Device >                     View ;
  typedef  typename MDArray< ValueType , Device >::HostMirror  HostMirror ;
  typedef  typename HostMirror::device_type                    HostDevice ;

  static
  HostMirror create( const View & a )
    {
      return create_labeled_mdarray< ValueType , HostDevice >(
               std::string() ,
               a.dimension(0) , a.dimension(1) ,
               a.dimension(2) , a.dimension(3) ,
               a.dimension(4) , a.dimension(5) ,
               a.dimension(6) , a.dimension(7) );
    }
};

}

template< typename ValueType , class DeviceType >
inline
typename MDArray< ValueType , DeviceType >::HostMirror
create_mirror( const MDArray< ValueType , DeviceType > & a )
{
  typedef typename MDArray< ValueType , DeviceType >::HostMirror  host_view ;
  typedef typename host_view::device_type                     host_device ;
  typedef typename host_device::memory_space                  host_memory ;
  typedef typename host_device::mdarray_map                   host_mdarray_map ;
  typedef typename DeviceType::memory_space                   memory ;
  typedef typename DeviceType::mdarray_map                    mdarray_map ;

  enum { OPTIMIZE = Impl::SameType< memory , host_memory >::value &&
                    Impl::SameType< mdarray_map , host_mdarray_map >::value &&
#if defined( KOKKOS_MIRROR_VIEW_OPTIMIZE )
                    KOKKOS_MIRROR_VIEW_OPTIMIZE
#else
                    false
#endif
       };

  return Impl::CreateMirror< MDArray< ValueType , DeviceType > , OPTIMIZE >
    ::create( a );
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

    Impl::DeepCopy< dst_type, src_type >::run( dst , src );
  }
}

//----------------------------------------------------------------------------

} // namespace Kokkos

#endif /* KOKKOS_MDARRAY_HPP */


