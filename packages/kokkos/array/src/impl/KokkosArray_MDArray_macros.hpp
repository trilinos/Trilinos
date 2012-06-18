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

#error "Including <impl/KokkosArray_MDArray_macros.hpp> without macros defined"

#else

namespace KokkosArray {

template< typename ValueType >
class MDArray< ValueType , KOKKOS_MACRO_DEVICE >
{
public:
  typedef ValueType                               value_type ;
  typedef KOKKOS_MACRO_DEVICE                     device_type ;
  typedef typename device_type::size_type         size_type ;
  typedef typename device_type::IndexMap<>::type  index_map ;

  typedef MDArray< value_type , HostMapped< device_type >::type > HostMirror ;

  /*------------------------------------------------------------------*/
  /** \brief  Query rank of the array */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type rank() const { return m_map.rank(); }

  /** \brief  Query dimension of the given ordinate of the array */
  template < typename iType >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension( const iType & rank_ordinate ) const
  { return m_map.dimension( rank_ordinate ); }

  template < typename iType >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  void dimensions( iType * const dims ) const
  { m_map.dimensions( dims ); }

  /** \brief  Query rank of the array */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type size() const { return m_map.size(); }

  /*------------------------------------------------------------------*/

#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )

  /** \brief  Because memory is contiguous this is exposed */
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type * ptr_on_device() const
  {
    // TBD: If memory is not contiguous and can throw then throw !
    return m_memory.ptr_on_device();
  }

  /*------------------------------------------------------------------*/
  /** \brief  Query value of a rank 8 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 , typename iType7 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 , const iType7 & i7 ) const
  { return m_memory.ptr_on_device()[ m_map.offset(iP,i1,i2,i3,i4,i5,i6,i7) ]; }

  /** \brief  Query value of a rank 7 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 ) const
  { return m_memory.ptr_on_device()[ m_map.offset(iP,i1,i2,i3,i4,i5,i6) ]; }

  /** \brief  Query value of a rank 6 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ) const
  { return m_memory.ptr_on_device()[ m_map.offset(iP,i1,i2,i3,i4,i5) ]; }

  /** \brief  Query value of a rank 5 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 ) const
  { return m_memory.ptr_on_device()[ m_map.offset(iP,i1,i2,i3,i4) ]; }

  /** \brief  Query value of a rank 4 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ) const
  { return m_memory.ptr_on_device()[ m_map.offset(iP,i1,i2,i3) ]; }

  /** \brief  Query value of a rank 3 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 ) const
  { return m_memory.ptr_on_device()[ m_map.offset(iP,i1,i2) ]; }

  /** \brief  Query value of a rank 2 array */
  template< typename iTypeP , typename iType1 >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ) const
  { return m_memory.ptr_on_device()[ m_map.offset(iP,i1) ]; }

  /** \brief  Query value of a rank 1 array */
  template< typename iTypeP >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP ) const
  { return m_memory.ptr_on_device()[ m_map.offset(iP) ]; }

  /*------------------------------------------------------------------*/

#endif /* defined( KOKKOS_MACRO_DEVICE_FUNCTION ) */

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MDArray() : m_memory(), m_map() {}

  /** \brief  Construct a view of the array */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MDArray( const MDArray & rhs )
    : m_memory( rhs.m_memory ), m_map( rhs.m_map ) {}

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MDArray & operator = ( const MDArray & rhs )
    {
      m_memory.operator=( rhs.m_memory );
      m_map   .operator=( rhs.m_map );
      return *this ;
    }

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ~MDArray() {}

  /*------------------------------------------------------------------*/
  /** \brief  Assign view if memory space and map are compatible.
   *          Intended use is to optimize mirroring by eliminating
   *          unnecessary deep copies.
   */
  template< class DeviceRHS >
  inline
  explicit
  MDArray( const MDArray< value_type , DeviceRHS > & rhs )
    : m_memory( rhs.m_memory )
    , m_map( rhs.m_map.dimension(0), rhs.m_map.dimension(1),
             rhs.m_map.dimension(2), rhs.m_map.dimension(3),
             rhs.m_map.dimension(4), rhs.m_map.dimension(5),
             rhs.m_map.dimension(6), rhs.m_map.dimension(7) )
  {
    typedef typename DeviceRHS::index_map rhs_index_map ;

    enum { SameMap = Impl::SameType< index_map , rhs_index_map >::value };

    Impl::StaticAssert< SameMap >::ok();
  }
  /*------------------------------------------------------------------*/
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  operator bool () const
  { return m_memory.operator bool(); }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator == ( const MDArray & rhs ) const
  { return m_memory.operator == ( rhs.m_memory ); }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator != ( const MDArray & rhs ) const
  { return m_memory.operator != ( rhs.m_memory ); }

  /*------------------------------------------------------------------*/
  template< class DeviceRHS >
  inline
  bool operator == ( const MDArray< value_type , DeviceRHS > & rhs ) const
  { return m_memory.operator == ( rhs.m_memory ); }

  template< class DeviceRHS >
  inline
  bool operator != ( const MDArray< value_type , DeviceRHS > & rhs ) const
  { return m_memory.operator != ( rhs.m_memory ); }

  /*------------------------------------------------------------------*/

private:

  typedef typename device_type::memory_space  memory_space ;
  typedef Impl::MemoryView< value_type , memory_space >  memory_view ;

  memory_view  m_memory ;
  index_map    m_map ;

  template< typename V , class D >
  friend
  MDArray< V , D >
  create_labeled_mdarray( const std::string & label ,
                          size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                          size_t n4 , size_t n5 , size_t n6 , size_t n7 );

  template< typename V , class D >  friend class MDArray ;
  template< class , class >  friend class Impl::Factory ;
};

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif


