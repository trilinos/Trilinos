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

#ifndef KOKKOS_MULTIVECTORVIEW_HPP
#define KOKKOS_MULTIVECTORVIEW_HPP

#include <cstddef>
#include <string>
#include <Kokkos_ArrayForwardDeclarations.hpp>
#include <impl/Kokkos_ArrayBounds.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
MultiVectorView< ValueType , DeviceType >
create_labeled_multivector( const std::string & label ,
                            size_t length , size_t count = 1 );

//----------------------------------------------------------------------------
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
class MultiVectorView {
public:
  typedef ValueType                       value_type ;
  typedef DeviceType                      device_type ;
  typedef typename DeviceType::size_type  size_type ;

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
  MultiVectorView();

  /** \brief  Construct a view of the array */
  MultiVectorView( const MultiVectorView & rhs );

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  MultiVectorView & operator = ( const MultiVectorView & rhs );
  
  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  ~MultiVectorView();

  /*------------------------------------------------------------------*/
  /** \brief View to a single vector */

  MultiVectorView( const MultiVectorView & rhs , size_type iV );

  MultiVectorView( const MultiVectorView & rhs , size_type iVbeg , size_type iVend );

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  operator bool () const ;

  /** \brief  Query if view to same memory */
  bool operator == ( const MultiVectorView & ) const ;

  /** \brief  Query if not view to same memory */
  bool operator != ( const MultiVectorView & ) const ;
private:

  MultiVectorView( const std::string & label ,
                   size_type length , size_type count );

  template< typename V , class M >
  friend
  MultiVectorView< V , M >
  create_labeled_multivector( const std::string & label ,
                              size_t length , size_t count );
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  The master creation method that all other versions call */
template< typename ValueType , class DeviceType >
inline
MultiVectorView< ValueType , DeviceType >
create_labeled_multivector( const std::string & label ,
                            size_t length , size_t count )
{
  return MultiVectorView< ValueType , DeviceType >( label , length , count );
}

//----------------------------------------------------------------------------

template< class MultiVectorType >
inline
MultiVectorView< typename MultiVectorType::value_type ,
                 typename MultiVectorType::device_type >
create_labeled_multivector( const std::string & label ,
                            size_t length , size_t count = 1 )
{
  return create_labeled_multivector< typename MultiVectorType::value_type ,
                                     typename MultiVectorType::device_type >
         ( label , length , count );
}

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
inline
MultiVectorView< ValueType , DeviceType >
create_multivector( size_t length , size_t count = 1 )
{
  return create_labeled_multivector< ValueType , DeviceType >
    ( std::string() , length , count );
}

template< class MultiVectorType >
inline
MultiVectorView< typename MultiVectorType::value_type ,
                 typename MultiVectorType::device_type >
create_multivector( size_t length , size_t count = 1 )
{
  return create_labeled_multivector< typename MultiVectorType::value_type ,
                                     typename MultiVectorType::device_type >
         ( std::string() , length , count );
}

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceDst , class DeviceSrc >
class MultiVectorDeepCopy ;

template< typename ValueType , class DeviceDst , class DeviceSrc >
inline
void deep_copy( const MultiVectorView< ValueType , DeviceDst > & dst ,
                const MultiVectorView< ValueType , DeviceSrc > & src )
{ MultiVectorDeepCopy<ValueType,DeviceDst,DeviceSrc>::run( dst , src ); }

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Partial specializations for known devices
  
#if defined( KOKKOS_DEVICE_HOST )
#include <Kokkos_DeviceHost_macros.hpp>
#include <impl/Kokkos_MultiVectorView_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>
#endif
  
#if defined( KOKKOS_DEVICE_TPI )
#include <Kokkos_DeviceTPI_macros.hpp>
#include <impl/Kokkos_MultiVectorView_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>
#endif

#if defined( KOKKOS_DEVICE_CUDA )
#include <Kokkos_DeviceCuda_macros.hpp>
#include <impl/Kokkos_MultiVectorView_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>
#include <DeviceCuda/Kokkos_DeviceCuda_MultiVectorView.hpp>
#endif

#endif /* KOKKOS_MULTIVECTORVIEW_HPP */


