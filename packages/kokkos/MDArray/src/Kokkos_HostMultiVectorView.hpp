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

#ifndef KOKKOS_HOSTMVECTORVIEW_HPP
#define KOKKOS_HOSTMVECTORVIEW_HPP

#include <Kokkos_ViewTracker.hpp>
#include <Kokkos_MVectorView.hpp>
#include <Kokkos_HostDevice.hpp>
#include <Kokkos_HostMap.hpp>

namespace Kokkos {

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
template< typename ValueType >
class MVectorView< ValueType , HostMap > {
public:
  typedef HostMap                        device_map_type ;
  typedef typename HostMap::device_type  device_type ;
  typedef typename HostMap::size_type    size_type ;
  typedef ValueType                      value_type ;

  /*------------------------------------------------------------------*/
  /** \brief  Query length of vectors */
  size_type length() const { return m_length ; }

  /** \brief  Query count of vectors */
  size_type count()  const { return m_count ; }

  /*------------------------------------------------------------------*/
  /** \brief  Query value */
  template< typename iTypeP , typename iTypeV >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP , const iTypeV & iV ) const
    {
      KOKKOS_BOUNDS_CHECK( require_less( iP , m_length ) );
      KOKKOS_BOUNDS_CHECK( require_less( iV , m_count ) );

      return m_ptr_on_device[ iP + m_length * iV ];
    }

  template< typename iTypeP >
  KOKKOS_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP ) const
    {
      KOKKOS_BOUNDS_CHECK( require_less( iP , m_length ) );

      return m_ptr_on_device[ iP ];
    }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  inline
  MVectorView()
    : m_tracker()
    , m_alloc_ptr_on_device( NULL )
    , m_ptr_on_device( NULL )
    , m_length( 0 )
    , m_count( 0 )
    {}

  /** \brief  Construct view to an existing multivector */
  MVectorView( const MVectorView & rhs );

  /** \brief  Clear the old view and assign this view
   *          to the 'rhs' multivector.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  MVectorView & operator = ( const MVectorView & rhs );
  
  /**  \brief  Destroy this view of the multivector.
   *           If the last view then allocated memory is deallocated.
   */
  ~MVectorView();

  /*------------------------------------------------------------------*/
  /** \brief Construct view to a single vector */
  MVectorView( const MVectorView & rhs , size_type iV );

  /** \brief Construct view to a multivector range */
  MVectorView( const MVectorView & rhs , size_type iVbeg ,
                                         size_type iVend );

  /*------------------------------------------------------------------*/

  value_type * address_on_device() const { return m_ptr_on_device ; }

  /*------------------------------------------------------------------*/

private:

  /*------------------------------------------------------------------*/

  ViewTracker  m_tracker ;             ///< Track equivalent views.
  value_type * m_alloc_ptr_on_device ; ///< Allocated memory
  value_type * m_ptr_on_device ;       ///< Referenced memory for subviews
  size_type    m_length ;              ///< Length of vectors
  size_type    m_count ;               ///< Count of vectors
  
  /*------------------------------------------------------------------*/
  /** \brief  Another view to the same memory. */
  inline
  void insert_view( const MVectorView & rhs ,
                    size_type iVbeg , size_type iVend )
    {
      KOKKOS_BOUNDS_CHECK( require_less( iVbeg , iVend ) );
      KOKKOS_BOUNDS_CHECK( require_less( iVend , rhs.m_count + 1 ) );

      m_tracker.insert( rhs.m_tracker );

      m_alloc_ptr_on_device = rhs.m_alloc_ptr_on_device ;
      m_ptr_on_device       = rhs.m_ptr_on_device + rhs.m_length * iVbeg ;
      m_length              = rhs.m_length ;
      m_count               = iVend - iVbeg ;
    }

  /** \brief  Remove this view from the shared views.
   *          If the last view then deallocate memory.
   */ 
  inline
  void clear_view()
    {
      if ( m_tracker.remove_and_query_is_last() ) {
        // If the last view then destroy the memory
        device_type::deallocate_memory( m_alloc_ptr_on_device );
      }
      m_alloc_ptr_on_device = NULL ;
      m_ptr_on_device       = NULL ;
      m_length              = 0 ;
      m_count               = 0 ;
    }

  /*------------------------------------------------------------------*/

  template< typename V , class M >
  friend
  MVectorView< V , M >
  create_mvector( typename M::size_type length ,
                  typename M::size_type count );

  template< typename V , class M >
  friend
  MVectorView< V , M >
  create_labeled_mvector( typename M::size_type length ,
                          typename M::size_type count ,
                          const std::string & label );

  /** \brief  Constructor that allocates */
  MVectorView( size_type arg_length , size_type arg_count ,
               const std::string & label )
    : m_tracker()
    , m_alloc_ptr_on_device( NULL )
    , m_ptr_on_device( NULL )
    , m_length( 0 )
    , m_count( 0 )
    {
      // 'allocate_memory' throws an exception if allocation fails.
      const size_t total = arg_length * arg_count ;
      value_type * const ptr = (value_type *)
        device_type::allocate_memory( sizeof(ValueType) , total , label );

      m_tracker.insert( m_tracker ); // First view
      m_alloc_ptr_on_device = ptr ;
      m_ptr_on_device       = ptr ;
      m_length              = arg_length ;
      m_count               = arg_count ;
    }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename ValueType >
inline
MVectorView< ValueType , HostMap >::
MVectorView( const MVectorView< ValueType , HostMap > & rhs )
  : m_tracker() { insert_view( rhs , 0 , rhs.m_count ); }

template< typename ValueType >
inline
MVectorView< ValueType , HostMap > &
MVectorView< ValueType , HostMap >::
operator = ( const MVectorView< ValueType , HostMap > & rhs )
  { clear_view(); insert_view( rhs , 0 , rhs.m_count ); return *this ; }

template< typename ValueType >
inline
MVectorView< ValueType , HostMap >::~MVectorView()
  { clear_view(); }

template< typename ValueType >
inline
MVectorView< ValueType , HostMap >::
MVectorView( const MVectorView< ValueType , HostMap > & rhs , 
             MVectorView< ValueType , HostMap >::size_type iV )
 : m_tracker()
 { insert_view( rhs , iV , iV + 1 ); }

template< typename ValueType >
inline
MVectorView< ValueType , HostMap >::
MVectorView( const MVectorView< ValueType , HostMap > & rhs , 
             MVectorView< ValueType , HostMap >::size_type iVbeg ,
             MVectorView< ValueType , HostMap >::size_type iVend )
 : m_tracker()
 { insert_view( rhs , iVbeg , iVend ); }

//----------------------------------------------------------------------------

template< typename ValueType >
struct MVectorDeepCopy< ValueType , HostMap , HostMap > {

  typedef MVectorView< ValueType , HostMap > MVector ;
  typedef HostMap::size_type                 size_type ;

  static 
  void run( const MVector & dest , const MVector & src )
  {
    require_equal( dest.length() , src.length() );
    require_equal( dest.count() ,  src.count() );

    const size_type n = dest.length() * dest.count();

    ValueType * d = dest.address_on_device();
    const ValueType * s = src.address_on_device();
    const ValueType * sEnd = s + n ;

    while ( s < sEnd ) { *d++ = *s++ ; }
  }
};

//----------------------------------------------------------------------------


} // namespace Kokkos

#endif /* KOKKOS_HOSTMVECTORVIEW_HPP */


