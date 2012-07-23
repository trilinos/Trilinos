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

#if ! defined(KOKKOS_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOS_MACRO_DEVICE)                  || \
    ! defined(KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <KokkosArray_PrefixSum_macros.hpp> without macros defined"

#else

namespace KokkosArray {

//----------------------------------------------------------------------------

template< typename IntType , class LayoutType >
class PrefixSum< IntType , LayoutType , KOKKOS_MACRO_DEVICE >
{
public:
  typedef KOKKOS_MACRO_DEVICE    device_type ;
  typedef IntType                size_type ;

  typedef PrefixSum< size_type , LayoutType , Host > HostMirror ;

  /*------------------------------------------------------------------*/

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type length() const { return m_data.dimension_0() - 1 ; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type sum() const { return m_sum ; }

  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  size_type operator[]( const iType & i ) const
    { return m_data[i]; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  const size_type * ptr_on_device() const
    { return m_data.ptr_on_device(); }
  
  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  PrefixSum()
    : m_data()
    , m_sum(0)
    {}

  /** \brief  Construct a view of the array */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  PrefixSum( const PrefixSum & rhs )
    : m_data( rhs.m_data )
    , m_sum(  rhs.m_sum )
    {}

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  PrefixSum & operator = ( const PrefixSum & rhs )
  {
    m_data = rhs.m_data ;
    m_sum  = rhs.m_sum ;
    return *this ;
  }

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ~PrefixSum() { m_sum = 0 ; }

  /*------------------------------------------------------------------*/

  operator bool () const
  { return m_data.ptr_on_device() != 0 ; }

  /** \brief  Query if view to same memory */
  bool operator == ( const PrefixSum & rhs ) const
  { return m_data.ptr_on_device() == rhs.m_data.ptr_on_device() ; }

  /** \brief  Query if not view to same memory */
  bool operator != ( const PrefixSum & rhs ) const
  { return m_data.ptr_on_device() != rhs.m_data.ptr_on_device() ; }

  template < class DeviceRHS >
  bool operator == ( const PrefixSum<size_type,DeviceRHS> & rhs ) const
  { return m_data.ptr_on_device() == rhs.m_data.ptr_on_device() ; }
  
  template < class DeviceRHS >
  bool operator != ( const PrefixSum<size_type,DeviceRHS> & rhs ) const
  { return m_data.ptr_on_device() != rhs.m_data.ptr_on_device() ; }
  

private:

  typedef View< size_type[] , LayoutType , device_type > view_type ;

  view_type  m_data ;
  size_type  m_sum ;

  template< typename , class , class > friend class PrefixSum ;
  template< class Dst , class Src >  friend class Impl::Factory ;
};

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

