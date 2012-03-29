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

#error "Including <impl/Kokkos_MDArray_macros.hpp> without macros defined"

#else

namespace Kokkos {

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

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< typename ValueType >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > , ValueType , 8 >
{
  typedef KOKKOS_MACRO_DEVICE::memory_space           device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type              size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE >  dst_type ;
  typedef ValueType                                   src_type ;

  const dst_type dst ;
  const src_type src ;
  const size_type N1 , N2 , N3 , N4 , N5 , N6 , N7 ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    , N1( arg_dst.dimension(1) )
    , N2( arg_dst.dimension(2) )
    , N3( arg_dst.dimension(3) )
    , N4( arg_dst.dimension(4) )
    , N5( arg_dst.dimension(5) )
    , N6( arg_dst.dimension(6) )
    , N7( arg_dst.dimension(7) )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < N5 ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < N6 ; ++i6 ) {
    for ( size_type i7 = 0 ; i7 < N7 ; ++i7 ) {
      dst(i0,i1,i2,i3,i4,i5,i6,i7) = src ;
    }}}}}}}
  }
};

template< typename ValueType >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > , ValueType , 7 >
{
  typedef KOKKOS_MACRO_DEVICE::memory_space           device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type              size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE >  dst_type ;
  typedef ValueType                                   src_type ;

  const dst_type dst ;
  const src_type src ;
  const size_type N1 , N2 , N3 , N4 , N5 , N6 ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    , N1( arg_dst.dimension(1) )
    , N2( arg_dst.dimension(2) )
    , N3( arg_dst.dimension(3) )
    , N4( arg_dst.dimension(4) )
    , N5( arg_dst.dimension(5) )
    , N6( arg_dst.dimension(6) )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < N5 ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < N6 ; ++i6 ) {
      dst(i0,i1,i2,i3,i4,i5,i6) = src ;
    }}}}}}
  }
};

template< typename ValueType >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > , ValueType , 6 >
{
  typedef KOKKOS_MACRO_DEVICE::memory_space           device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type              size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE >  dst_type ;
  typedef ValueType                                   src_type ;

  const dst_type dst ;
  const src_type src ;
  const size_type N1 , N2 , N3 , N4 , N5 ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    , N1( arg_dst.dimension(1) )
    , N2( arg_dst.dimension(2) )
    , N3( arg_dst.dimension(3) )
    , N4( arg_dst.dimension(4) )
    , N5( arg_dst.dimension(5) )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < N5 ; ++i5 ) {
      dst(i0,i1,i2,i3,i4,i5) = src ;
    }}}}}
  }
};

template< typename ValueType >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > , ValueType , 5 >
{
  typedef KOKKOS_MACRO_DEVICE::memory_space           device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type              size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE >  dst_type ;
  typedef ValueType                                   src_type ;

  const dst_type dst ;
  const src_type src ;
  const size_type N1 , N2 , N3 , N4 ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    , N1( arg_dst.dimension(1) )
    , N2( arg_dst.dimension(2) )
    , N3( arg_dst.dimension(3) )
    , N4( arg_dst.dimension(4) )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
      dst(i0,i1,i2,i3,i4) = src ;
    }}}}
  }
};

template< typename ValueType >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > , ValueType , 4 >
{
  typedef KOKKOS_MACRO_DEVICE::memory_space           device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type              size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE >  dst_type ;
  typedef ValueType                                   src_type ;

  const dst_type dst ;
  const src_type src ;
  const size_type N1 , N2 , N3 ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    , N1( arg_dst.dimension(1) )
    , N2( arg_dst.dimension(2) )
    , N3( arg_dst.dimension(3) )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
      dst(i0,i1,i2,i3) = src ;
    }}}
  }
};

template< typename ValueType >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > , ValueType , 3 >
{
  typedef KOKKOS_MACRO_DEVICE::memory_space           device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type              size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE >  dst_type ;
  typedef ValueType                                   src_type ;

  const dst_type dst ;
  const src_type src ;
  const size_type N1 , N2 ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    , N1( arg_dst.dimension(1) )
    , N2( arg_dst.dimension(2) )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
      dst(i0,i1,i2) = src ;
    }}
  }
};

template< typename ValueType >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > , ValueType , 2 >
{
  typedef KOKKOS_MACRO_DEVICE::memory_space           device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type              size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE >  dst_type ;
  typedef ValueType                                   src_type ;

  const dst_type dst ;
  const src_type src ;
  const size_type N1 ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    , N1( arg_dst.dimension(1) )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
      dst(i0,i1) = src ;
    }
  }
};

template< typename ValueType >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > , ValueType , 1 >
{
  typedef KOKKOS_MACRO_DEVICE::memory_space           device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type              size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE >  dst_type ;
  typedef ValueType                                   src_type ;

  const dst_type dst ;
  const src_type src ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    dst(i0) = src ;
  }
};

//----------------------------------------------------------------------------

template< typename ValueType , class SrcDevice >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > ,
  MDArray< ValueType , SrcDevice > , 8 >
{
  enum { OK = StaticAssert<
                SameType< KOKKOS_MACRO_DEVICE::memory_space ,
                          typename  SrcDevice::memory_space >::value >::value };

  typedef KOKKOS_MACRO_DEVICE::memory_space          device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type             size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE > dst_type ;
  typedef MDArray< ValueType , SrcDevice >           src_type ;

  const dst_type dst ;
  const src_type src ;
  const size_type N1 , N2 , N3 , N4 , N5 , N6 , N7 ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    , N1( arg_dst.dimension(1) )
    , N2( arg_dst.dimension(2) )
    , N3( arg_dst.dimension(3) )
    , N4( arg_dst.dimension(4) )
    , N5( arg_dst.dimension(5) )
    , N6( arg_dst.dimension(6) )
    , N7( arg_dst.dimension(7) )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < N5 ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < N6 ; ++i6 ) {
    for ( size_type i7 = 0 ; i7 < N7 ; ++i7 ) {
      dst(i0,i1,i2,i3,i4,i5,i6,i7) = src(i0,i1,i2,i3,i4,i5,i6,i7);
    }}}}}}}
  }
};

template< typename ValueType , class SrcDevice >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > ,
  MDArray< ValueType , SrcDevice > , 7 >
{
  enum { OK = StaticAssert<
                SameType< KOKKOS_MACRO_DEVICE::memory_space ,
                          typename  SrcDevice::memory_space >::value >::value };

  typedef KOKKOS_MACRO_DEVICE::memory_space          device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type             size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE > dst_type ;
  typedef MDArray< ValueType , SrcDevice >           src_type ;

  const dst_type dst ;
  const src_type src ;
  const size_type N1 , N2 , N3 , N4 , N5 , N6 ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    , N1( arg_dst.dimension(1) )
    , N2( arg_dst.dimension(2) )
    , N3( arg_dst.dimension(3) )
    , N4( arg_dst.dimension(4) )
    , N5( arg_dst.dimension(5) )
    , N6( arg_dst.dimension(6) )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < N5 ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < N6 ; ++i6 ) {
      dst(i0,i1,i2,i3,i4,i5,i6) = src(i0,i1,i2,i3,i4,i5,i6);
    }}}}}}
  }
};

template< typename ValueType , class SrcDevice >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > ,
  MDArray< ValueType , SrcDevice > , 6 >
{
  enum { OK = StaticAssert<
                SameType< KOKKOS_MACRO_DEVICE::memory_space ,
                          typename  SrcDevice::memory_space >::value >::value };

  typedef KOKKOS_MACRO_DEVICE::memory_space          device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type             size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE > dst_type ;
  typedef MDArray< ValueType , SrcDevice >           src_type ;

  const dst_type dst ;
  const src_type src ;
  const size_type N1 , N2 , N3 , N4 , N5 ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    , N1( arg_dst.dimension(1) )
    , N2( arg_dst.dimension(2) )
    , N3( arg_dst.dimension(3) )
    , N4( arg_dst.dimension(4) )
    , N5( arg_dst.dimension(5) )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < N5 ; ++i5 ) {
      dst(i0,i1,i2,i3,i4,i5) = src(i0,i1,i2,i3,i4,i5);
    }}}}}
  }
};

template< typename ValueType , class SrcDevice >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > ,
  MDArray< ValueType , SrcDevice > , 5 >
{
  enum { OK = StaticAssert<
                SameType< KOKKOS_MACRO_DEVICE::memory_space ,
                          typename  SrcDevice::memory_space >::value >::value };

  typedef KOKKOS_MACRO_DEVICE::memory_space          device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type             size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE > dst_type ;
  typedef MDArray< ValueType , SrcDevice >           src_type ;

  const dst_type dst ;
  const src_type src ;
  const size_type N1 , N2 , N3 , N4 ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    , N1( arg_dst.dimension(1) )
    , N2( arg_dst.dimension(2) )
    , N3( arg_dst.dimension(3) )
    , N4( arg_dst.dimension(4) )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
      dst(i0,i1,i2,i3,i4) = src(i0,i1,i2,i3,i4);
    }}}}
  }
};

template< typename ValueType , class SrcDevice >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > ,
  MDArray< ValueType , SrcDevice > , 4 >
{
  enum { OK = StaticAssert<
                SameType< KOKKOS_MACRO_DEVICE::memory_space ,
                          typename  SrcDevice::memory_space >::value >::value };

  typedef KOKKOS_MACRO_DEVICE::memory_space          device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type             size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE > dst_type ;
  typedef MDArray< ValueType , SrcDevice >           src_type ;

  const dst_type dst ;
  const src_type src ;
  const size_type N1 , N2 , N3 ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    , N1( arg_dst.dimension(1) )
    , N2( arg_dst.dimension(2) )
    , N3( arg_dst.dimension(3) )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
      dst(i0,i1,i2,i3) = src(i0,i1,i2,i3);
    }}}
  }
};

template< typename ValueType , class SrcDevice >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > ,
  MDArray< ValueType , SrcDevice > , 3 >
{
  enum { OK = StaticAssert<
                SameType< KOKKOS_MACRO_DEVICE::memory_space ,
                          typename  SrcDevice::memory_space >::value >::value };

  typedef KOKKOS_MACRO_DEVICE::memory_space          device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type             size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE > dst_type ;
  typedef MDArray< ValueType , SrcDevice >           src_type ;

  const dst_type dst ;
  const src_type src ;
  const size_type N1 , N2 ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    , N1( arg_dst.dimension(1) )
    , N2( arg_dst.dimension(2) )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
      dst(i0,i1,i2) = src(i0,i1,i2);
    }}
  }
};

template< typename ValueType , class SrcDevice >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > ,
  MDArray< ValueType , SrcDevice > , 2 >
{
  enum { OK = StaticAssert<
                SameType< KOKKOS_MACRO_DEVICE::memory_space ,
                          typename  SrcDevice::memory_space >::value >::value };

  typedef KOKKOS_MACRO_DEVICE::memory_space          device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type             size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE > dst_type ;
  typedef MDArray< ValueType , SrcDevice >           src_type ;

  const dst_type dst ;
  const src_type src ;
  const size_type N1 ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    , N1( arg_dst.dimension(1) )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
      dst(i0,i1) = src(i0,i1);
    }
  }
};

template< typename ValueType , class SrcDevice >
struct DeepCopyKernelMDArray<
  MDArray< ValueType , KOKKOS_MACRO_DEVICE > ,
  MDArray< ValueType , SrcDevice > , 1 >
{
  enum { OK = StaticAssert<
                SameType< KOKKOS_MACRO_DEVICE::memory_space ,
                          typename  SrcDevice::memory_space >::value >::value };

  typedef KOKKOS_MACRO_DEVICE::memory_space          device_type ;
  typedef KOKKOS_MACRO_DEVICE::size_type             size_type ;
  typedef MDArray< ValueType , KOKKOS_MACRO_DEVICE > dst_type ;
  typedef MDArray< ValueType , SrcDevice >           src_type ;

  const dst_type dst ;
  const src_type src ;

  DeepCopyKernelMDArray( const dst_type & arg_dst , 
                         const src_type & arg_src )
    : dst( arg_dst ), src( arg_src )
    {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    dst(i0) = src(i0);
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif


