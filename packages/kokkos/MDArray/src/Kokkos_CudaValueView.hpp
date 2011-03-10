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

#ifndef KOKKOS_CUDAVALUEVIEW_HPP
#define KOKKOS_CUDAVALUEVIEW_HPP

#include <Kokkos_ValueView.hpp>
#include <Kokkos_ViewTracker.hpp>
#include <Kokkos_CudaDevice.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------
/** \brief  Plain-old-data value allocated on a compute device.
 */
template< typename ValueType >
class ValueView<ValueType,CudaDevice> {
public:
  typedef CudaDevice device_type ;
  typedef ValueType  value_type ;

  /*------------------------------------------------------------------*/
  /** \brief  Access value */
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

  /** \brief  Allow the ValueView to be a parallel reduce
   *          'finalize functor' that assigns the reduced value
   *          on the device.
   */
  KOKKOS_DEVICE_FUNCTION
  void operator()( const value_type & rhs )
  { *m_ptr_on_device = rhs ; }

  KOKKOS_DEVICE_AND_HOST_FUNCTION
  ValueType * address_on_device() const { return m_ptr_on_device ; }

private:

  ViewTracker  m_tracker ;
  value_type * m_ptr_on_device ;

  /*------------------------------------------------------------------*/

  inline
  void insert_view( const ValueView & rhs )
    {
#ifndef __CUDA_ARCH__
      m_tracker.insert( rhs.m_tracker );
#endif
      m_ptr_on_device = rhs.m_ptr_on_device ;
    }
  inline

  void clear_view()
    {
#ifndef __CUDA_ARCH__
      if ( m_tracker.remove_and_query_is_last() ) {
        // If the last view then destroy the memory
        device_type::deallocate_memory( m_ptr_on_device );
      }
#endif
      m_ptr_on_device = NULL ; 
    }

  /*------------------------------------------------------------------*/

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
struct ValueDeepCopy<ValueType,CudaDevice> {

  static void run( const ValueView<ValueType,CudaDevice> & dest ,
                   const ValueType & src )
  {
    // ( device , host , size , code )
    cudaMemcpy( dest.address_on_device() , & src , sizeof(ValueType) , cudaMemcpyHostToDevice );
  }

  static void run( ValueType & dest ,
                   const ValueView<ValueType,CudaDevice> & src )
  {
    cudaThreadSynchronize(); // Wait for device
    cudaMemcpy( & dest , src.address_on_device() , sizeof(ValueType) , cudaMemcpyDeviceToHost );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Kokkos

#endif /* KOKKOS_CUDAVALUEVIEW_HPP */


