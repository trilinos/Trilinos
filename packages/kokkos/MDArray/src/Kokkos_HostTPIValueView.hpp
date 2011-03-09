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

#ifndef KOKKOS_HOSTTPIVALUEVIEW_HPP
#define KOKKOS_HOSTTPIVALUEVIEW_HPP

#include <Kokkos_ValueView.hpp>
#include <Kokkos_ViewTracker.hpp>
#include <Kokkos_HostTPI.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------
/** \brief  Plain-old-data value allocated on a compute device.
 */
template< typename ValueType >
class ValueView<ValueType,HostTPI> {
public:
  typedef HostTPI   device_type ;
  typedef ValueType value_type ;

  /*------------------------------------------------------------------*/
  /** \brief  Query value */
  KOKKOS_DEVICE_FUNCTION
  value_type & operator * () const { return *m_ptr_on_device ; }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  ValueView() : m_tracker(), m_ptr_on_device( NULL ) {}

  /** \brief  Construct a view of the array */
  ValueView( const ValueView & rhs )
    : m_tracker() { insert_view( rhs ); }

  /** \brief  Assign to a view of the rhs.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  ValueView & operator = ( const ValueView & rhs )
    { clear_view(); insert_view( rhs ); return *this ; }
  
  /**  \brief  Destroy this view of the value.
   *           If the last view then allocated memory is deallocated.
   */
  ~ValueView() { clear_view(); }

  /*------------------------------------------------------------------*/

  KOKKOS_DEVICE_AND_HOST_FUNCTION
  ValueType * address_on_device() const { return m_ptr_on_device ; }

private:

  ViewTracker  m_tracker ;
  value_type * m_ptr_on_device ;

  /*------------------------------------------------------------------*/

  inline
  void insert_view( const ValueView & rhs )
    {
      m_tracker.insert( rhs.m_tracker );
      m_ptr_on_device = rhs.m_ptr_on_device ;
    }

  inline
  void clear_view()
    {
      if ( m_tracker.remove_and_query_is_last() ) {
        // If the last view then destroy the memory
        device_type::deallocate_memory( m_ptr_on_device );
      }
      m_ptr_on_device = NULL ; 
    }

  /*------------------------------------------------------------------*/

private:

  template< typename V , class D >
  friend
  ValueView< V , D > create_value();

  template< typename V , class D >
  friend
  ValueView< V , D > create_labeled_value( const std::string & label );

  ValueView( const std::string & label )
  : m_tracker(), m_ptr_on_device( NULL )
  {
    m_ptr_on_device = (ValueType *) device_type::allocate_memory( sizeof(ValueType) , 1 , label );
    m_tracker.insert( m_tracker ); // First view
  }
};

//----------------------------------------------------------------------------

template< typename ValueType >
struct ValueDeepCopy<ValueType,HostTPI> {

  static void run( const ValueView<ValueType,HostTPI> & dest ,
                   const ValueType & src )
  { *dest = src ; }

  static void run( ValueType & dest ,
                   const ValueView<ValueType,HostTPI> & src )
  { dest = *src ; }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Kokkos

#endif /* KOKKOS_HOSTTPIVALUEVIEW_HPP */


