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

#error "Including <impl/Kokkos_IndexMapLeft_macros.hpp> without macros defined"

#else

#include <impl/Kokkos_IndexMap.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** 
 *  Compile time dimension, except for leading dimension.
 */
template< unsigned Rank , unsigned N1 , unsigned N2 , unsigned N3 ,
          unsigned N4 ,   unsigned N5 , unsigned N6 , unsigned N7 >
struct IndexMapLeft< KOKKOS_MACRO_DEVICE::memory_space , 
                     Rank, N1, N2, N3, N4, N5, N6, N7 >
{
public:
  typedef KOKKOS_MACRO_DEVICE::size_type size_type ;

private:
  static const unsigned NSize =
    1 < Rank ? 1 : N1 * (
    2 < Rank ? 1 : N2 * (
    3 < Rank ? 1 : N3 * (
    4 < Rank ? 1 : N4 * (
    5 < Rank ? 1 : N5 * (
    6 < Rank ? 1 : N6 * (
    7 < Rank ? 1 : N7 ))))));


  static const unsigned S2 = N1 ;
  static const unsigned S3 = N2 * S2 ;
  static const unsigned S4 = N3 * S3 ;
  static const unsigned S5 = N4 * S4 ;
  static const unsigned S6 = N5 * S5 ;
  static const unsigned S7 = N6 * S6 ;

  size_type N0 ;
  size_type NAlloc ;

public:

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type rank() const { return Rank ; }

  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension( const iType & ordinal ) const
  {
    KOKKOS_MACRO_CHECK( require_less( ordinal , Rank ) );
    const size_type d[] = { N0, N1, N2, N3, N4, N5, N6, N7 };
    return d[ordinal];
  }

  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  void dimensions( iType * const dims ) const
  {
    const size_type d[] = { N0, N1, N2, N3, N4, N5, N6, N7 };
    for ( size_type i = 0 ; i < Rank ; ++i ) { dims[i] = d[i] ; }
  }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type size() const { return N0 * NSize ; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type allocation_size() const { return NAlloc * NSize ; }

  //------------------------------

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 , typename iType5 ,
             typename iType6 , typename iType7 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 , const iType5 & i5 ,
                    const iType6 & i6 , const iType7 & i7 ) const
  {
    enum { InputRank = 8 };

    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( Rank, N0, N1, N2, N3, N4, N5, N6, N7,
                            InputRank, i0, i1, i2, i3, i4, i5, i6, i7 ) );

    return i0 + NAlloc * ( i1 + S2 * i2 + S3 * i3 + S4 * i4 +
                                S5 * i5 + S6 * i6 + S7 * i7 );
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 , typename iType5 ,
             typename iType6 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 , const iType5 & i5 ,
                    const iType6 & i6 ) const
  {
    enum { InputRank = 7 };

    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( Rank, N0, N1, N2, N3, N4, N5, N6, 0,
                            InputRank, i0, i1, i2, i3, i4, i5, i6, 0 ) );

    return i0 + NAlloc * ( i1 + S2 * i2 + S3 * i3 + S4 * i4 +
                                S5 * i5 + S6 * i6 );
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 , typename iType5 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 , const iType5 & i5 ) const
  {
    enum { InputRank = 6 };

    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( Rank, N0, N1, N2, N3, N4, N5, 0, 0,
                            InputRank, i0, i1, i2, i3, i4, i5, 0, 0 ) );

    return i0 + NAlloc * ( i1 + S2 * i2 + S3 * i3 + S4 * i4 + S5 * i5 );
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 ) const
  {
    enum { InputRank = 5 };

    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( Rank, N0, N1, N2, N3, N4, 0, 0, 0,
                            InputRank, i0, i1, i2, i3, i4, 0, 0, 0 ) );

    return i0 + NAlloc * ( i1 + S2 * i2 + S3 * i3 + S4 * i4 );
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ) const
  {
    enum { InputRank = 4 };

    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( Rank, N0, N1, N2, N3, 0, 0, 0, 0,
                            InputRank, i0, i1, i2, i3, 0, 0, 0, 0 ) );

    return i0 + NAlloc * ( i1 + S2 * i2 + S3 * i3 );
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 ) const
  {
    enum { InputRank = 3 };

    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( Rank, N0, N1, N2, 0, 0, 0, 0, 0,
                            InputRank, i0, i1, i2, 0, 0, 0, 0, 0 ) );

    return i0 + NAlloc * ( i1 + S2 * i2 );
  }

  template < typename iType0 , typename iType1 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 , const iType1 & i1 ) const
  {
    enum { InputRank = 2 };

    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( Rank, N0, N1, 0, 0, 0, 0, 0, 0,
                            InputRank, i0, i1, 0, 0, 0, 0, 0, 0 ) );

    return i0 + NAlloc * ( i1 );
  }

  template < typename iType0 , typename iType1 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 ) const
  {
    enum { InputRank = 1 };

    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( Rank, N0, 0, 0, 0, 0, 0, 0, 0,
                            InputRank, i0, 0, 0, 0, 0, 0, 0, 0 ) );

    return i0 ;
  }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  IndexMapLeft() : N0(0), NAlloc(0) {}

  template< typename ValueType >
  inline
  void assign( size_t n0 )
    {
      typedef KOKKOS_MACRO_DEVICE::memory_space   memory_space ;
      typedef Impl::MemoryManager< memory_space > memory_manager ;

      N0 = n0 ;
      NAlloc = memory_manager::preferred_alignment< ValueType >( n0 );
    }
};

//----------------------------------------------------------------------------

/** \brief  Runtime dimensioned */
template<>
struct IndexMapLeft< KOKKOS_MACRO_DEVICE::memory_space, 0, 0,0,0,0,0,0,0 >
{
public:
  typedef KOKKOS_MACRO_DEVICE::size_type size_type ;

private:


  // Stride for the rank #0 index,
  // If contiguous then the rank #0 dimension.
  // If padded to the alignment then rank #0 dimension + padding.
  static const unsigned S0 = IndexMapMaxRank ;

  size_type m_rank ;
  size_type m_dims[ IndexMapMaxRank + 1 ];

public:

  /** \brief  Rank of the index space */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type rank() const { return m_rank ; }

  /** \brief  Dimension of a given rank of the index space */
  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension( const iType & ordinal ) const
  {
    KOKKOS_MACRO_CHECK( require_less( ordinal , IndexMapMaxRank ) );
    return m_dims[ordinal];
  }

  /** \brief  Dimension of all ranks of the index space.
   *          The 'dims' output array must be at least of length 'rank()'.
   */
  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  void dimensions( iType * const dims ) const
  {
    for ( size_type i = 0 ; i < m_rank ; ++i ) {
      dims[i] = m_dims[i] ;
    }
  }

  /** \brief  The size of the index space; i.e., produce of dimensions. */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type size() const
    {
      size_type n = m_dims[0]; //  The actual rank #0 dimension
      for ( size_type i = 1 ; i < m_rank ; ++i ) { n *= m_dims[i] ; }
      return n ;
    }

  /** \brief  Amount of memory to allocate for the array
   *          including the necessary padding.
   */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type allocation_size() const
    {
      size_type n = m_dims[ S0 ]; // The padded rank #0 dimension
      for ( size_type i = 1 ; i < m_rank ; ++i ) { n *= m_dims[i] ; }
      return n ;
    }

  //------------------------------------

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 , typename iType5 ,
             typename iType6 , typename iType7 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 , const iType5 & i5 ,
                    const iType6 & i6 , const iType7 & i7 ) const
  {
    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( m_rank,
                                 m_dims[0], m_dims[1], m_dims[2], m_dims[3],
                                 m_dims[4], m_dims[5], m_dims[6], m_dims[7],
                                 8, i0, i1, i2, i3, i4, i5, i6, i7 ) );

    return ( i0 + m_dims[S0] * ( i1 + m_dims[1] *
           ( i2 + m_dims[2]  * ( i3 + m_dims[3] *
           ( i4 + m_dims[4]  * ( i5 + m_dims[5] *
           ( i6 + m_dims[6]  * ( i7 ))))))));
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 , typename iType5 ,
             typename iType6 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 , const iType5 & i5 ,
                    const iType6 & i6 ) const
  {
    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( m_rank,
                                 m_dims[0], m_dims[1], m_dims[2], m_dims[3],
                                 m_dims[4], m_dims[5], m_dims[6], m_dims[7],
                                 7, i0, i1, i2, i3, i4, i5, i6, 0 ) );

    return ( i0 + m_dims[S0] * ( i1 + m_dims[1] *
           ( i2 + m_dims[2]  * ( i3 + m_dims[3] *
           ( i4 + m_dims[4]  * ( i5 + m_dims[5] *
           ( i6 )))))));
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 , typename iType5 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 , const iType5 & i5 ) const
  {
    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( m_rank,
                                 m_dims[0], m_dims[1], m_dims[2], m_dims[3],
                                 m_dims[4], m_dims[5], m_dims[6], m_dims[7],
                                 6, i0, i1, i2, i3, i4, i5, 0, 0 ) );

    return ( i0 + m_dims[S0] * ( i1 + m_dims[1] *
           ( i2 + m_dims[2]  * ( i3 + m_dims[3] *
           ( i4 + m_dims[4]  * ( i5 ))))));
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 ,
             typename iType4 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ,
                    const iType4 & i4 ) const
  {
    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( m_rank,
                                 m_dims[0], m_dims[1], m_dims[2], m_dims[3],
                                 m_dims[4], m_dims[5], m_dims[6], m_dims[7],
                                 5, i0, i1, i2, i3, i4, 0, 0, 0 ) );

    return ( i0 + m_dims[S0] * ( i1 + m_dims[1] *
           ( i2 + m_dims[2]  * ( i3 + m_dims[3] *
           ( i4 )))));
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 , typename iType3 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 , const iType3 & i3 ) const
  {
    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( m_rank,
                                 m_dims[0], m_dims[1], m_dims[2], m_dims[3],
                                 m_dims[4], m_dims[5], m_dims[6], m_dims[7],
                                 4, i0, i1, i2, i3, 0, 0, 0, 0 ) );

    return ( i0 + m_dims[S0] * ( i1 + m_dims[1] *
           ( i2 + m_dims[2]  * ( i3 ))));
  }

  template < typename iType0 , typename iType1 ,
             typename iType2 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 , const iType1 & i1 ,
                    const iType2 & i2 ) const
  {
    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( m_rank,
                                 m_dims[0], m_dims[1], m_dims[2], m_dims[3],
                                 m_dims[4], m_dims[5], m_dims[6], m_dims[7],
                                 3, i0, i1, i2, 0, 0, 0, 0, 0 ) );

    return ( i0 + m_dims[S0] * ( i1 + m_dims[1] * ( i2 )));
  }

  template < typename iType0 , typename iType1 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 , const iType1 & i1 ) const
  {
    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( m_rank,
                                 m_dims[0], m_dims[1], m_dims[2], m_dims[3],
                                 m_dims[4], m_dims[5], m_dims[6], m_dims[7],
                                 2, i0, i1, 0, 0, 0, 0, 0, 0 ) );

    return ( i0 + m_dims[S0] * ( i1 ));
  }

  template < typename iType0 >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const iType0 & i0 ) const
  {
    KOKKOS_MACRO_CHECK(
      mdarray_require_dimension( m_rank,
                                 m_dims[0], m_dims[1], m_dims[2], m_dims[3],
                                 m_dims[4], m_dims[5], m_dims[6], m_dims[7],
                                 1, i0, 0, 0, 0, 0, 0, 0, 0 ) );

    return ( i0 );
  }

  //------------------------------------

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  void reverse_index( size_type offset, size_type * const indices ) const
  {
    for ( size_type i = 0 ; i < m_rank ; ++i ) {
      indices[i] = offset % m_dims[i] ;
      offset /= m_dims[i] ;
    }
  }

  //--------------------------------------

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  IndexMapLeft() : m_rank(0)
  {
    m_dims[0] = 0 ; m_dims[1] = 0 ;
    m_dims[2] = 0 ; m_dims[3] = 0 ;
    m_dims[4] = 0 ; m_dims[5] = 0 ;
    m_dims[6] = 0 ; m_dims[7] = 0 ;
    m_dims[S0] = 0 ;
  }

  template< typename ValueType >
  inline
  void assign( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
               size_t n4 , size_t n5 , size_t n6 , size_t n7 )
  {
    typedef KOKKOS_MACRO_DEVICE::memory_space   memory_space ;
    typedef Impl::MemoryManager< memory_space > memory_manager ;

    m_rank = mdarray_deduce_rank( n0, n1, n2, n3, n4, n5, n6, n7 );
    m_dims[0] = n0 ; m_dims[1] = n1 ; m_dims[2] = n2 ; m_dims[3] = n3 ;
    m_dims[4] = n4 ; m_dims[5] = n5 ; m_dims[6] = n6 ; m_dims[7] = n7 ;
    m_dims[S0] = memory_manager::preferred_alignment< ValueType >( n0 );
  }
};

} // Impl namespace
} // Kokkos namespace

#endif

