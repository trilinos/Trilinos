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
#include <impl/KokkosArray_forward.hpp>
#include <impl/KokkosArray_MemoryView.hpp>
#include <impl/KokkosArray_StaticAssert.hpp>
#include <impl/KokkosArray_ArrayTraits.hpp>
#include <impl/KokkosArray_ArrayBounds.hpp>

namespace KokkosArray {

//----------------------------------------------------------------------------
/// \class MultiVector
/// \brief View of a collection of one or more vectors, allocated and mapped onto a compute device.
///
/// The first rank corresponds to the parallel work index.
/// The second rank selects one vector of the multivector.
/// All vectors in the multivector have the same length (number of entries).
///
/// No assumptions should be made as to the mapping, contiguity, or
/// strides of the storage of these multivectors.  The mapping will
/// vary according to the underlying device.  The usage model is for
/// algorithms to be parameterized with respect to the type of the
/// mapped multivector and thus achieve portability across compute
/// devices.
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

  /// \brief Get the value at work index iP and vector index iV.
  ///
  /// \warning For a DeviceType which is not a host device, this
  ///   method should only be called in a compute kernel.
  template< typename iTypeP , typename iTypeV >
  value_type & operator()( const iTypeP & iP , const iTypeV & iV ) const ;

  /// \brief Fetch the value at work index iP (if there is only one vector).
  ///
  /// \warning For a DeviceType which is not a host device, this
  ///   method should only be called in a compute kernel.
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

  /// \brief Constructor producing a view of a single vector in rhs.
  ///
  /// \param iV [in] Index of the vector of rhs to view.
  MultiVector( const MultiVector & rhs , size_type iV );

  /// \brief Constructor producing a view of a contiguous range of vectors in rhs.
  ///
  /// \param iVbeg [in] Index of the first vector in the range.
  /// \param iVend [in] Exclusive index of the last vector in the range.
  ///   (For example, for the range 0, 1, ... N-1, iVbeg=0 and iVend=N.)
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

template< class MultiVectorType >
MultiVector< typename MultiVectorType::value_type ,
             typename MultiVectorType::device_type >
create_multivector(
  const std::string & label , size_t length , size_t count = 1 )
{
  typedef MultiVector< typename MultiVectorType::value_type ,
                       typename MultiVectorType::device_type > type ;
  return Impl::Factory< type , void >::create( label , length , count );
}

template< class MultiVectorType >
MultiVector< typename MultiVectorType::value_type ,
             typename MultiVectorType::device_type >
create_multivector( const size_t length , size_t count = 1 )
{
  typedef MultiVector< typename MultiVectorType::value_type ,
                       typename MultiVectorType::device_type > type ;
  return Impl::Factory< type , void >::create( std::string() , length , count );
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

    Impl::Factory<dst_type,src_type>::deep_copy( dst , src );
  }
}

template< typename ValueType , class DeviceDst , class DeviceSrc >
inline
void deep_copy( const MultiVector< ValueType , DeviceDst > & dst ,
                const MultiVector< ValueType , DeviceSrc > & src ,
                const size_t count )
{
  typedef MultiVector< ValueType , DeviceDst > dst_type ;
  typedef MultiVector< ValueType , DeviceSrc > src_type ;

  if ( dst.operator!=( src ) ) {
    Impl::Factory<dst_type,src_type>::deep_copy( dst , src , count );
  }
}

//----------------------------------------------------------------------------

template< typename ValueType , class Device >
inline
void update( const ValueType & alpha ,
             const MultiVector< ValueType , Device > & x ,
             const ValueType & beta ,
             const MultiVector< ValueType , Device > & y )
{
  typedef MultiVector< ValueType , Device > multivec_type ;

  if ( x.operator!=( y ) ) {
    Impl::multivector_require_equal_dimension( x.length() , x.count() ,
                                               y.length() , y.count() );
  }

  Impl::Update<multivec_type>::run( alpha , x , beta, y );
}


} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <impl/KokkosArray_MultiVector_factory.hpp>

#endif /* KOKKOS_MULTIVECTOR_HPP */

