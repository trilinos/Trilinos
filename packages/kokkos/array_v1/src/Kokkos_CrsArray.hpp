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

#ifndef KOKKOS_CRSARRAY_HPP
#define KOKKOS_CRSARRAY_HPP

#include <cstddef>
#include <string>
#include <utility>
#include <iterator>
#include <limits>
#include <impl/Kokkos_forward.hpp>
#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_ArrayBounds.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------
/** \brief  Compressed row storage array.
 *
 *  map ( I , J ) -> value
 *  where I in [ 0 .. N ) and J in [ 0 .. M_I )
 */
template< typename ValueType , class DeviceType >
class CrsArray {
public:
  typedef DeviceType                       device_type ;
  typedef ValueType                        value_type ;
  typedef typename device_type::size_type  size_type ;

  typedef CrsArray< value_type , void /* Host */ > HostView ;

  /*------------------------------------------------------------------*/
  /** \brief  Number of rows */
  size_type row_dimension() const ;

  /** \brief  Number of columns in a given row,
   *          only available on the device.
   */
  template< typename iType >
  size_type column_dimension( const iType & row ) const ;

  /** \brief  Total number of values,
   *          only available on the device.
   */
  size_type size() const ;

  /*------------------------------------------------------------------*/
  /** \brief  Query value */
  template< typename iType , typename jType >
  value_type & operator()( const iType & row , const jType & j ) const ;
  
  typedef std::pair<value_type*,value_type*> value_range_type ;
  
  /** \brief  Query a row of values */
  template< typename iType >
  value_range_type operator()( const iType & row ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  CrsArray();

  /** \brief  Construct a view of the array */
  CrsArray( const CrsArray & rhs );

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  CrsArray & operator = ( const CrsArray & rhs );

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  ~CrsArray();

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  operator bool () const ;

  /** \brief  Query if view to same memory */
  bool operator == ( const CrsArray & ) const ;

  /** \brief  Query if not view to same memory */
  bool operator != ( const CrsArray & ) const ;
};

//----------------------------------------------------------------------------
/** \brief  Create a CrsArray with the input row sizes.
 *
 *  Constructed array has the following properties:
 *
 *  row_dimension()       == std::distance( row_size_begin , row_size_end );
 *  column_dimension(row) == *row_iter ;
 *    where  row_iter = row_size_begin ; std::advance( row_iter , row );
 */
template< typename ValueType , class DeviceType , typename IteratorType >
CrsArray< ValueType , DeviceType >
create_labeled_crsarray( const std::string & label ,
                         const IteratorType row_size_begin ,
                         const IteratorType row_size_end );

template< typename ValueType , class DeviceType , typename IteratorType >
inline
CrsArray< ValueType , DeviceType >
create_crsarray( const std::string & label ,
                 const IteratorType row_size_begin ,
                 const IteratorType row_size_end )
{ return create_crsarray( std::string() , row_size_begin , row_size_end ); }

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType , typename IteratorType >
typename CrsArray< ValueType , DeviceType >::HostView
create_mirror( const CrsArray< ValueType , DeviceType > & );

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceDst , class DeviceSrc >
void deep_copy( const CrsArray< ValueType , DeviceDst > & dst ,
                const CrsArray< ValueType , DeviceSrc > & src );

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Impl {

template< typename ValueType , class Device >
class CreateCrsArray {
public:
  template< typename IteratorType >
  static CrsArray< ValueType , Device > create( const std::string & ,
                                                const IteratorType ,
                                                const IteratorType );
};

template< typename ValueType , class Device >
class CreateMirror< CrsArray< ValueType , Device > , true /* view */ >
{
public:
  typedef  CrsArray< ValueType , Device >  View ;
  typedef  typename View::HostView         HostView ;

  static
  HostView create( const View & v ) { return HostView( v ); }
};

template< typename ValueType , class Device >
class CreateMirror< CrsArray< ValueType , Device > , false /* copy */ >
{
public:
  typedef  CrsArray< ValueType , Device >  View ;
  typedef  typename View::HostView         HostView ;
  typedef  typename HostView::device_type  HostDevice ;

  static
  HostView create( const View & v );
};

} // namespace Impl

//----------------------------------------------------------------------------
/** \brief  Create a CrsArray with the input row sizes.
 *
 *  Constructed array has the following properties:
 *
 *  row_dimension()       == std::distance( row_size_begin , row_size_end );
 *  column_dimension(row) == *row_iter ;
 *    where  row_iter = row_size_begin ; std::advance( row_iter , row );
 */
template< typename ValueType , class DeviceType , typename IteratorType >
CrsArray< ValueType , DeviceType >
create_labeled_crsarray( const std::string & label ,
                         const IteratorType row_size_begin ,
                         const IteratorType row_size_end )
{
  // Try to verify the iterators dereference to integers

  typedef std::iterator_traits<IteratorType> traits ;
  typedef typename traits::value_type value_type ;
  enum { iterator_to_integer = std::numeric_limits<value_type>::is_integer };
  enum { OK = Impl::StaticAssert< iterator_to_integer >::value };

  return Impl::CreateCrsArray< ValueType , DeviceType >
    ::create( label , row_size_begin , row_size_end );
}

template< typename ValueType , class DeviceType >
inline
typename CrsArray< ValueType , DeviceType >::HostView
create_mirror( const CrsArray< ValueType , DeviceType > & v )
{
  typedef CrsArray< ValueType , DeviceType >     view_type ;
  typedef typename view_type::HostView           host_view ;
  typedef typename host_view::device_type        host_device ;
  typedef typename host_device::memory_space     host_memory ;
  typedef typename DeviceType::memory_space      memory ;

  enum { optimize = Impl::SameType< memory , host_memory >::value &&
#if defined( KOKKOS_MIRROR_VIEW_OPTIMIZE )
                    KOKKOS_MIRROR_VIEW_OPTIMIZE
#else
                    false
#endif
       };

  return Impl::CreateMirror< view_type , optimize >::create( v );
}


template< typename ValueType , class DeviceDst , class DeviceSrc >
inline
void deep_copy( const CrsArray< ValueType , DeviceDst > & dst ,
                const CrsArray< ValueType , DeviceSrc > & src )
{
  typedef CrsArray< ValueType , DeviceDst > dst_type ;
  typedef CrsArray< ValueType , DeviceSrc > src_type ;

  if ( dst.operator!=( src ) ) {
    Impl::crsarray_require_equal_dimension( dst.row_dimension() , dst.size() ,
                                            src.row_dimension() , src.size() );

    Impl::DeepCopy<dst_type,src_type>::run( dst , src );
  }
}


//----------------------------------------------------------------------------

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_CRSARRAY_HPP */

