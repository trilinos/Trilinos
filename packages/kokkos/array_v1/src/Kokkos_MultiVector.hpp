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

#ifndef KOKKOS_MULTIVECTOR_HPP
#define KOKKOS_MULTIVECTOR_HPP

#include <cstddef>
#include <string>
#include <impl/Kokkos_forward.hpp>
#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_ArrayBounds.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------
/** \brief  Multivector allocated and mapped
 *          onto a compute device.
 *
 *  The first rank corresponds to the parallel work index.
 *  The second rank selects one vector of the multivector.
 *
 *  No assumptions should be made as to the mapping, contiguity, or strides
 *  of the storage of these multivectors.  The mapping will vary according to the
 *  underlying device.  The usage model is for algorithms to be parameterized
 *  with respect to the type of the mapped multivector and thus achieve portability
 *  across compute devices.
 */
template< typename ValueType , class DeviceType >
class MultiVector {
public:
  typedef ValueType                       value_type ;
  typedef DeviceType                      device_type ;
  typedef typename DeviceType::size_type  size_type ;

  typedef MultiVector< value_type , void /* Host */ > HostMirror ;

  /*------------------------------------------------------------------*/
  /** \brief  Query length of vectors */
  size_type length() const ;

  /** \brief  Query count of vectors */
  size_type count() const ;

  /*------------------------------------------------------------------*/
  /** \brief  Query value */
  template< typename iTypeP , typename iTypeV >
  value_type & operator()( const iTypeP & iP , const iTypeV & iV ) const ;

  template< typename iTypeP >
  value_type & operator()( const iTypeP & iP ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  MultiVector();

  /** \brief  Construct a view of the array */
  MultiVector( const MultiVector & rhs );

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  MultiVector & operator = ( const MultiVector & rhs );

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  ~MultiVector();

  /*------------------------------------------------------------------*/
  /** \brief View to a single vector */

  MultiVector( const MultiVector & rhs , size_type iV );

  MultiVector( const MultiVector & rhs , size_type iVbeg , size_type iVend );

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  operator bool () const ;

  /** \brief  Query if view to same memory */
  bool operator == ( const MultiVector & ) const ;

  /** \brief  Query if not view to same memory */
  bool operator != ( const MultiVector & ) const ;
};

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
MultiVector< ValueType , DeviceType >
create_labeled_multivector( const std::string & label ,
                            size_t length , size_t count = 1 );

template< class MultiVectorType >
inline
MultiVector< typename MultiVectorType::value_type ,
             typename MultiVectorType::device_type >
create_labeled_multivector( const std::string & label ,
                            size_t length , size_t count = 1 );

template< typename ValueType , class DeviceType >
MultiVector< ValueType , DeviceType >
create_multivector( size_t length , size_t count = 1 );

template< class MultiVectorType >
inline
MultiVector< typename MultiVectorType::value_type ,
             typename MultiVectorType::device_type >
create_multivector( size_t length , size_t count = 1 );

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
typename MultiVector< ValueType , DeviceType >::HostMirror
create_mirror( const MultiVector< ValueType , DeviceType > & );

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceDst , class DeviceSrc >
void deep_copy( const MultiVector< ValueType , DeviceDst > & dst ,
                const MultiVector< ValueType , DeviceSrc > & src );

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  The master creation method that all other versions call */
template< typename ValueType , class DeviceType >
inline
MultiVector< ValueType , DeviceType >
create_labeled_multivector( const std::string & label ,
                            size_t length , size_t count )
{
  return MultiVector< ValueType , DeviceType >( label , length , count );
}

//----------------------------------------------------------------------------

template< class MultiVectorType >
inline
MultiVector< typename MultiVectorType::value_type ,
             typename MultiVectorType::device_type >
create_labeled_multivector( const std::string & label ,
                            size_t length , size_t count )
{
  return create_labeled_multivector< typename MultiVectorType::value_type ,
                                     typename MultiVectorType::device_type >
         ( label , length , count );
}

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
inline
MultiVector< ValueType , DeviceType >
create_multivector( size_t length , size_t count )
{
  return create_labeled_multivector< ValueType , DeviceType >
    ( std::string() , length , count );
}

template< class MultiVectorType >
inline
MultiVector< typename MultiVectorType::value_type ,
             typename MultiVectorType::device_type >
create_multivector( size_t length , size_t count )
{
  return create_labeled_multivector< typename MultiVectorType::value_type ,
                                     typename MultiVectorType::device_type >
         ( std::string() , length , count );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Impl {

template< typename ValueType , class Device >
class CreateMirror< MultiVector< ValueType , Device > , true /* view */ >
{
public:
  typedef  MultiVector< ValueType , Device >            View ;
  typedef  typename MultiVector< ValueType , Device >::HostMirror  HostMirror ;

  static
  HostMirror create( const View & v ) { return HostMirror( v ); }
};

template< typename ValueType , class Device >
class CreateMirror< MultiVector< ValueType , Device > , false /* copy */ >
{
public:
  typedef  MultiVector< ValueType , Device >            View ;
  typedef  typename MultiVector< ValueType , Device >::HostMirror  HostMirror ;
  typedef  typename HostMirror::device_type  HostDevice ;

  static
  HostMirror create( const View & v )
    {
      const size_t length = v.length();
      const size_t count  = v.count();
      return create_labeled_multivector< ValueType , HostDevice  >( std::string() , length , count );
    }
};

}

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
inline
typename MultiVector< ValueType , DeviceType >::HostMirror
create_mirror( const MultiVector< ValueType , DeviceType > & v )
{
  typedef MultiVector< ValueType , DeviceType >  view_type ;
  typedef typename view_type::HostMirror           host_view ;
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

  return Impl::CreateMirror< MultiVector< ValueType , DeviceType > , optimize >
    ::create( v );
}

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceDst , class DeviceSrc >
inline
void deep_copy( const MultiVector< ValueType , DeviceDst > & dst ,
                const MultiVector< ValueType , DeviceSrc > & src )
{
  typedef MultiVector< ValueType , DeviceDst > dst_type ;
  typedef MultiVector< ValueType , DeviceSrc > src_type ;

  if ( dst.operator!=( src ) ) {
    Impl::multivector_require_equal_dimension( dst.length() , dst.count() ,
                                               src.length() , src.count() );

    Impl::DeepCopy<dst_type,src_type>::run( dst , src );
  }
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_MULTIVECTOR_HPP */

