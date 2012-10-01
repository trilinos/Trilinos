/*
//@HEADER
// ************************************************************************
//
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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

#ifndef KOKKOSARRAY_VIEW_HPP
#define KOKKOSARRAY_VIEW_HPP

#include <string>
#include <impl/KokkosArray_StaticAssert.hpp>
#include <impl/KokkosArray_ArrayTraits.hpp>
#include <impl/KokkosArray_Shape.hpp>
#include <impl/KokkosArray_AnalyzeShape.hpp>

#include <KokkosArray_Macros.hpp>

namespace KokkosArray {

class Host ;

//----------------------------------------------------------------------------

template < class >
struct View_requires_a_device_specific_specialization_which_is_not_found ;

template < class >
struct View_create_requires_non_const_data_type ;

//----------------------------------------------------------------------------

template< class DataType , class LayoutType , class DeviceType = LayoutType >
class View {
public:
  typedef DataType    data_type ;
  typedef LayoutType  layout_type ;
  typedef DeviceType  device_type ;

  typedef View< data_type , layout_type , Host >  HostMirror ;

  typedef typename Impl::AnalyzeShape<DataType>::value_type  value_type ;
  typedef typename LayoutType::array_layout                  array_layout ;
  typedef typename device_type::memory_space                 memory_space ;
  typedef typename device_type::size_type                    size_type ;

private:

  typedef typename Impl::AnalyzeShape<DataType,LayoutType>::shape shape_type ;

public:

  typedef Impl::unsigned_< shape_type::rank > Rank ;

  size_type rank() const ;

  template< typename iType >
  size_type dimension( const iType & rank ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Value for rank-1 array */
  template< typename iType0 >
  value_type & operator()( const iType0 & i0 ) const ;

  /** \brief  Value for rank-2 array */
  template< typename iType0 , typename iType1 >
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ) const ;

  /** \brief  Value for rank-3 array */
  template< typename iType0 , typename iType1 , typename iType2 >
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ) const ;

  /** \brief  Value for rank-4 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 >
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ) const ;

  /** \brief  Value for rank-5 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 >
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ) const ;

  /** \brief  Value for rank-6 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ,
                           const iType5 & i5 ) const ;

  /** \brief  Value for rank-7 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 >
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ,
                           const iType5 & i5 ,
                           const iType6 & i6 ) const ;

  /** \brief  Value for rank-8 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 , typename iType7 >
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ,
                           const iType5 & i5 ,
                           const iType6 & i6 ,
                           const iType7 & i7 ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  View()
  {
    View_requires_a_device_specific_specialization_which_is_not_found<device_type>::fail();
  }

  /** \brief  Construct a view of the array */
  View( const View & rhs )
  {
    View_requires_a_device_specific_specialization_which_is_not_found<device_type>::fail();
  }

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  View & operator = ( const View & rhs );

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  ~View();

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  operator bool () const ;

  /** \brief  Query if view to same memory */
  bool operator == ( const View & ) const ;

  /** \brief  Query if not view to same memory */
  bool operator != ( const View & ) const ;

  /*------------------------------------------------------------------*/

  explicit View( const std::string & )
  {
    View_requires_a_device_specific_specialization_which_is_not_found<device_type>::fail();
  }

  View( const std::string & , const unsigned N0 )
  {
    View_requires_a_device_specific_specialization_which_is_not_found<device_type>::fail();
  }

  View( const std::string & , const unsigned N0 ,
                              const unsigned N1 )
  {
    View_requires_a_device_specific_specialization_which_is_not_found<device_type>::fail();
  }

  View( const std::string & , const unsigned N0 ,
                              const unsigned N1 ,
                              const unsigned N2 )
  {
    View_requires_a_device_specific_specialization_which_is_not_found<device_type>::fail();
  }

  View( const std::string & , const unsigned N0 ,
                              const unsigned N1 ,
                              const unsigned N2 ,
                              const unsigned N3 )
  {
    View_requires_a_device_specific_specialization_which_is_not_found<device_type>::fail();
  }

  View( const std::string & , const unsigned N0 ,
                              const unsigned N1 ,
                              const unsigned N2 ,
                              const unsigned N3 ,
                              const unsigned N4 )
  {
    View_requires_a_device_specific_specialization_which_is_not_found<device_type>::fail();
  }

  View( const std::string & , const unsigned N0 ,
                              const unsigned N1 ,
                              const unsigned N2 ,
                              const unsigned N3 ,
                              const unsigned N4 ,
                              const unsigned N5 )
  {
    View_requires_a_device_specific_specialization_which_is_not_found<device_type>::fail();
  }

  View( const std::string & , const unsigned N0 ,
                              const unsigned N1 ,
                              const unsigned N2 ,
                              const unsigned N3 ,
                              const unsigned N4 ,
                              const unsigned N5 ,
                              const unsigned N6 )
  {
    View_requires_a_device_specific_specialization_which_is_not_found<device_type>::fail();
  }

  View( const std::string & , const unsigned N0 ,
                              const unsigned N1 ,
                              const unsigned N2 ,
                              const unsigned N3 ,
                              const unsigned N4 ,
                              const unsigned N5 ,
                              const unsigned N6 ,
                              const unsigned N7 )
  {
    View_requires_a_device_specific_specialization_which_is_not_found<device_type>::fail();
  }
};

//----------------------------------------------------------------------------

template< class DataType , class LayoutType , class DeviceType >
typename View< DataType , LayoutType , DeviceType >::HostMirror
create_mirror( const View<DataType,LayoutType,DeviceType> & );

//----------------------------------------------------------------------------

template < class T >
struct View_create_requires_non_const_data_type { static void success(){} };

template < class T >
struct View_create_requires_non_const_data_type< const T > {};

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class MemorySpace , typename ValueType , class Shape >
class ViewOper ;

template< typename ValueType , class DeviceDst , class DeviceSrc >
struct DeepCopy ;

}
}

#include <impl/KokkosArray_View_factory.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_VIEW_HPP */

