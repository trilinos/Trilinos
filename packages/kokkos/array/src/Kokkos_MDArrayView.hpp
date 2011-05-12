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

#ifndef KOKKOS_MDARRAYVIEW_HPP
#define KOKKOS_MDARRAYVIEW_HPP

#include <cstddef>
#include <string>

namespace Kokkos {

//----------------------------------------------------------------------------

template< typename ValueType ,
          class DeviceType ,
          class MapOption = typename DeviceType::default_mdarray_map >
class MDArrayView ;

template< typename ValueType , class DeviceType , class MapOption >
MDArrayView< ValueType , DeviceType , MapOption >
create_labeled_mdarray( const std::string & label ,
                        size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                                    size_t n3 = 0 , size_t n4 = 0 ,
                                    size_t n5 = 0 , size_t n6 = 0 ,
                                    size_t n7 = 0 );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Multidimensional array allocated and mapped
 *          onto a compute device.
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
 *  create_labeled_mdarray< ValueType , DeviceType , MapOption >( nP , ... , label );
 *  create_labeled_mdarray< ValueType , DeviceType >( nP , ... , label );
 *  create_labeled_mdarray< MDArrayView<...> >( nP , ... , label );
 *
 *  The "unlabeled" group creates MDArray's with NULL string labels.
 *
 *  create_mdarray< ValueType , DeviceType , MapOption >( nP , ... );
 *  create_mdarray< ValueType , DeviceType >( nP , ... );
 *  create_mdarray< MDArrayView<...> >( nP , ... );
 */

template< typename ValueType , class DeviceType , class MapOption >
class MDArrayView {
public:
  typedef ValueType                       value_type ;
  typedef DeviceType                      device_type ;
  typedef MapOption                       map_option ;
  typedef typename DeviceType::size_type  size_type ;

  /*------------------------------------------------------------------*/
  /** \brief  Query rank of the array */
  size_type rank() const ;

  /** \brief  Query dimension of the given ordinate of the array */
  template < typename iType >
  size_type dimension( const iType & rank_ordinate ) const ;

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
  /** \brief  Construct a NULL view */
  inline
  MDArrayView();

  /** \brief  Construct another view of the 'rhs' array */
  inline
  MDArrayView( const MDArrayView & rhs );

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  inline
  MDArrayView & operator = ( const MDArrayView & rhs );
  
  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  inline
  ~MDArrayView();

private:

  MDArrayView( const std::string & label ,
               size_t nP , size_t n1 , size_t n2 , size_t n3 ,
               size_t n4 , size_t n5 , size_t n6 , size_t n7 );

  template< typename V , class D , class M >
  friend
  MDArrayView< V , D , M >
  create_labeled_mdarray( const std::string & label ,
                          size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                          size_t n4 , size_t n5 , size_t n6 , size_t n7 );
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  This is THE creation function.
 *          All other versions call this function.
 */
template< typename ValueType , class DeviceType , class MapType >
inline
MDArrayView< ValueType , DeviceType , MapType >
create_labeled_mdarray( const std::string & label ,
                        size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                        size_t n4 , size_t n5 , size_t n6 , size_t n7 )
{
  return MDArrayView< ValueType , DeviceType , MapType >
           ( label , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
inline
MDArrayView< ValueType , DeviceType , typename DeviceType::default_map >
create_labeled_mdarray( const std::string & label ,
                        size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                                    size_t n3 = 0 , size_t n4 = 0 ,
                                    size_t n5 = 0 , size_t n6 = 0 ,
                                    size_t n7 = 0 )
{
  return create_labeled_mdarray< ValueType , DeviceType , typename DeviceType::default_map >
           ( label , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------

template< typename MDArrayType >
inline
MDArrayView< typename MDArrayType::value_type ,
             typename MDArrayType::device_type ,
             typename MDArrayType::map_option >
create_labeled_mdarray( const std::string & label ,
                        size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                                    size_t n3 = 0 , size_t n4 = 0 ,
                                    size_t n5 = 0 , size_t n6 = 0 ,
                                    size_t n7 = 0 )
{
  return create_labeled_mdarray< typename MDArrayType::value_type ,
                                 typename MDArrayType::device_type ,
                                 typename MDArrayType::map_option >
           ( label , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType , class MapType >
inline
MDArrayView< ValueType , DeviceType , MapType >
create_mdarray( size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                            size_t n3 = 0 , size_t n4 = 0 ,
                            size_t n5 = 0 , size_t n6 = 0 ,
                            size_t n7 = 0 )
{
  return create_labeled_mdarray< ValueType , DeviceType , MapType >
           ( std::string() , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
inline
MDArrayView< ValueType , DeviceType , typename DeviceType::default_map >
create_mdarray( size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                            size_t n3 = 0 , size_t n4 = 0 ,
                            size_t n5 = 0 , size_t n6 = 0 ,
                            size_t n7 = 0 )
{
  return create_labeled_mdarray< ValueType , DeviceType , typename DeviceType::default_map >
           ( std::string() , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------

template< typename MDArrayType >
inline
MDArrayView< typename MDArrayType::value_type ,
             typename MDArrayType::device_type ,
             typename MDArrayType::map_option >
create_mdarray( size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                            size_t n3 = 0 , size_t n4 = 0 ,
                            size_t n5 = 0 , size_t n6 = 0 ,
                            size_t n7 = 0 )
{
  return create_labeled_mdarray< typename MDArrayType::value_type ,
                                 typename MDArrayType::device_type ,
                                 typename MDArrayType::map_option >
           ( std::string() , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <impl/Kokkos_ArrayBounds.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Partial specializations for known devices and default index maps

#if defined( KOKKOS_DEVICE_HOST )
#include <macros/Kokkos_DeviceHost_macros.hpp>
#include <macros/Kokkos_MDArrayIndexMapLeft_macros.hpp>
#include <macros/Kokkos_MDArrayView_macros.hpp>
#include <macros/Kokkos_DeviceClear_macros.hpp>
#endif

#if defined( KOKKOS_DEVICE_TPI )
#include <macros/Kokkos_DeviceTPI_macros.hpp>
#include <macros/Kokkos_MDArrayIndexMapLeft_macros.hpp>
#include <macros/Kokkos_MDArrayView_macros.hpp>
#include <macros/Kokkos_DeviceClear_macros.hpp>
#endif

#if defined( KOKKOS_DEVICE_CUDA )
#include <macros/Kokkos_DeviceCuda_macros.hpp>
#include <macros/Kokkos_MDArrayIndexMapLeft_macros.hpp>
#include <macros/Kokkos_MDArrayView_macros.hpp>
#include <macros/Kokkos_DeviceClear_macros.hpp>
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_MDARRAYVIEW_HPP */


