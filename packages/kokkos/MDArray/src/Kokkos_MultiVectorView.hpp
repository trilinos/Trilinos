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

#ifndef KOKKOS_BOUNDS_CHECK
#define KOKKOS_BOUNDS_CHECK( EXPR )  EXPR
#else
#define KOKKOS_BOUNDS_CHECK( EXPR )  /* EXPR */
#endif

namespace Kokkos {

//----------------------------------------------------------------------------


template< typename ValueType , class DeviceMapType > class MultiVectorView ;

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceMapType >
inline
MultiVectorView< ValueType , DeviceMapType >
create_multivector( typename DeviceMapType::size_type length ,
                    typename DeviceMapType::size_type count )
{
  return MultiVectorView< ValueType , DeviceMapType >( length , count , std::string() );
}

template< typename ValueType , class DeviceMapType >
inline
MultiVectorView< ValueType , DeviceMapType >
create_labeled_multivector( typename DeviceMapType::size_type length ,
                            typename DeviceMapType::size_type count ,
                            const std::string & label )
{
  return MultiVectorView< ValueType , DeviceMapType >( length , count , label );
}

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceMapDest , class DeviceMapSrc >
struct MultiVectorDeepCopy ;

template< typename ValueType , class DeviceMapDest , class DeviceMapSrc >
void deep_copy( const MultiVectorView<ValueType,DeviceMapDest> & dest ,
                const MultiVectorView<ValueType,DeviceMapSrc>  & src )
{
  MultiVectorDeepCopy<ValueType,DeviceMapDest,DeviceMapSrc>::run( dest , src );
}

//----------------------------------------------------------------------------

void require_less(  size_t i , size_t j );
void require_equal( size_t i , size_t j );

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
template< typename ValueType , class DeviceMapType >
class MultiVectorView {
public:
  typedef DeviceMapType                          device_map_type ;
  typedef typename device_map_type::device_type  device_type ;
  typedef typename device_map_type::size_type    size_type ;
  typedef ValueType                              value_type ;

  /*------------------------------------------------------------------*/
  /** \brief  Query length of vectors */
  size_type length() const ;

  /** \brief  Query count of vectors */
  size_type count()  const ;

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

private:

  template< typename V , class M >
  friend
  MultiVectorView< V , M >
  create_multivector( typename M::size_type length ,
                      typename M::size_type count );

  template< typename V , class M >
  friend
  MultiVectorView< V , M >
  create_labeled_multivector( typename M::size_type length ,
                              typename M::size_type count ,
                              const std::string & label );

  MultiVectorView( size_type length , size_type count , const std::string & label );
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Kokkos

#endif /* KOKKOS_MULTIVECTORVIEW_HPP */


