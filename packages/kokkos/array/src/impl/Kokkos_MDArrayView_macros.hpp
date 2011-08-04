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

#if ! defined(KOKKOS_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOS_MACRO_DEVICE)                  || \
    ! defined(KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <impl/Kokkos_MDArrayView_macros.hpp> without macros defined"

#else

namespace Kokkos {

template< typename ValueType >
class MDArrayView< ValueType , KOKKOS_MACRO_DEVICE >
{
public:
  typedef ValueType                           value_type ;
  typedef KOKKOS_MACRO_DEVICE                 device_type ;
  typedef typename device_type::mdarray_map   mdarray_map ;
  typedef typename device_type::memory_space  memory_space ;
  typedef typename device_type::size_type     size_type ;

  typedef MDArrayView< value_type , Serial< HostMemory , mdarray_map > > HostView ;

  /*------------------------------------------------------------------*/
  enum { Contiguous = true };

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
  value_type * ptr_on_device() const { return m_memory.ptr_on_device(); }

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

#endif /* defined( KOKKOS_MACRO_DEVICE_FUNCTION ) */

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MDArrayView() : m_memory(), m_map() {}

  /** \brief  Construct a view of the array */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MDArrayView( const MDArrayView & rhs )
    : m_memory(), m_map( rhs.m_map )
    { memory_space::assign_memory_view( m_memory , rhs.m_memory); }

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MDArrayView & operator = ( const MDArrayView & rhs )
    {
      memory_space::assign_memory_view( m_memory , rhs.m_memory );
      m_map = rhs.m_map ;
      return *this ;
    }
  
  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ~MDArrayView() {}

  /*------------------------------------------------------------------*/
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  operator bool () const
  { return m_memory.operator bool(); }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator == ( const MDArrayView & rhs ) const
  { return m_memory.operator == ( rhs.m_memory ); }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator != ( const MDArrayView & rhs ) const
  { return m_memory.operator != ( rhs.m_memory ); }

  /*------------------------------------------------------------------*/

private:

  MemoryView< value_type , memory_space >            m_memory ;
  Impl::MDArrayIndexMap< memory_space , mdarray_map > m_map ;

  inline
  MDArrayView( const std::string & label ,
               size_t nP , size_t n1 , size_t n2 , size_t n3 ,
               size_t n4 , size_t n5 , size_t n6 , size_t n7 )
    : m_memory()
    , m_map( nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 )
    { memory_space::allocate_memory_view( m_memory , m_map.size() , label ); }

  template< typename V , class D >
  friend
  MDArrayView< V , D >
  create_labeled_mdarray( const std::string & label ,
                          size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                          size_t n4 , size_t n5 , size_t n6 , size_t n7 );

  template< typename V , class DeviceDst , class DeviceSrc ,
            bool , bool , bool >
  friend
  class Impl::MDArrayDeepCopy ;
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )

namespace Kokkos {
namespace Impl {

template< class MDArrayDst , class MDArraySrc >
class MDArrayDeepCopyFunctor< KOKKOS_MACRO_DEVICE , MDArrayDst , MDArraySrc , 8 >
{
public:
  typedef KOKKOS_MACRO_DEVICE    device_type ;
  typedef device_type::size_type size_type ;

  MDArrayDst dst ; 
  MDArraySrc src ; 

  MDArrayDeepCopyFunctor( const MDArrayDst & arg_dst ,
                          const MDArraySrc & arg_src )
    : dst( arg_dst ), src( arg_src ) {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < dst.dimension(1) ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < dst.dimension(2) ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < dst.dimension(3) ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < dst.dimension(4) ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < dst.dimension(5) ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < dst.dimension(6) ; ++i6 ) {
    for ( size_type i7 = 0 ; i7 < dst.dimension(7) ; ++i7 ) {
      dst(i0,i1,i2,i3,i4,i5,i6,i7) = src(i0,i1,i2,i3,i4,i5,i6,i7);
    }}}}}}}
  }
};

template< class MDArrayDst , class MDArraySrc >
class MDArrayDeepCopyFunctor< KOKKOS_MACRO_DEVICE , MDArrayDst , MDArraySrc , 7 >
{
public:
  typedef KOKKOS_MACRO_DEVICE    device_type ;
  typedef device_type::size_type size_type ;

  MDArrayDst dst ; 
  MDArraySrc src ; 

  MDArrayDeepCopyFunctor( const MDArrayDst & arg_dst ,
                          const MDArraySrc & arg_src )
    : dst( arg_dst ), src( arg_src ) {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < dst.dimension(1) ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < dst.dimension(2) ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < dst.dimension(3) ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < dst.dimension(4) ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < dst.dimension(5) ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < dst.dimension(6) ; ++i6 ) {
      dst(i0,i1,i2,i3,i4,i5,i6) = src(i0,i1,i2,i3,i4,i5,i6);
    }}}}}}
  }
};

template< class MDArrayDst , class MDArraySrc >
class MDArrayDeepCopyFunctor< KOKKOS_MACRO_DEVICE , MDArrayDst , MDArraySrc , 6 >
{
public:
  typedef KOKKOS_MACRO_DEVICE    device_type ;
  typedef device_type::size_type size_type ;

  MDArrayDst dst ; 
  MDArraySrc src ; 

  MDArrayDeepCopyFunctor( const MDArrayDst & arg_dst ,
                          const MDArraySrc & arg_src )
    : dst( arg_dst ), src( arg_src ) {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < dst.dimension(1) ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < dst.dimension(2) ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < dst.dimension(3) ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < dst.dimension(4) ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < dst.dimension(5) ; ++i5 ) {
      dst(i0,i1,i2,i3,i4,i5) = src(i0,i1,i2,i3,i4,i5);
    }}}}}
  }
};

template< class MDArrayDst , class MDArraySrc >
class MDArrayDeepCopyFunctor< KOKKOS_MACRO_DEVICE , MDArrayDst , MDArraySrc , 5 >
{
public:
  typedef KOKKOS_MACRO_DEVICE    device_type ;
  typedef device_type::size_type size_type ;

  MDArrayDst dst ; 
  MDArraySrc src ; 

  MDArrayDeepCopyFunctor( const MDArrayDst & arg_dst ,
                          const MDArraySrc & arg_src )
    : dst( arg_dst ), src( arg_src ) {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < dst.dimension(1) ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < dst.dimension(2) ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < dst.dimension(3) ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < dst.dimension(4) ; ++i4 ) {
      dst(i0,i1,i2,i3,i4) = src(i0,i1,i2,i3,i4);
    }}}}
  }
};

template< class MDArrayDst , class MDArraySrc >
class MDArrayDeepCopyFunctor< KOKKOS_MACRO_DEVICE , MDArrayDst , MDArraySrc , 4 >
{
public:
  typedef KOKKOS_MACRO_DEVICE    device_type ;
  typedef device_type::size_type size_type ;

  MDArrayDst dst ; 
  MDArraySrc src ; 

  MDArrayDeepCopyFunctor( const MDArrayDst & arg_dst ,
                          const MDArraySrc & arg_src )
    : dst( arg_dst ), src( arg_src ) {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < dst.dimension(1) ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < dst.dimension(2) ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < dst.dimension(3) ; ++i3 ) {
      dst(i0,i1,i2,i3) = src(i0,i1,i2,i3);
    }}}
  }
};

template< class MDArrayDst , class MDArraySrc >
class MDArrayDeepCopyFunctor< KOKKOS_MACRO_DEVICE , MDArrayDst , MDArraySrc , 3 >
{
public:
  typedef KOKKOS_MACRO_DEVICE    device_type ;
  typedef device_type::size_type size_type ;

  MDArrayDst dst ; 
  MDArraySrc src ; 

  MDArrayDeepCopyFunctor( const MDArrayDst & arg_dst ,
                          const MDArraySrc & arg_src )
    : dst( arg_dst ), src( arg_src ) {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < dst.dimension(1) ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < dst.dimension(2) ; ++i2 ) {
      dst(i0,i1,i2) = src(i0,i1,i2);
    }}
  }
};

template< class MDArrayDst , class MDArraySrc >
class MDArrayDeepCopyFunctor< KOKKOS_MACRO_DEVICE , MDArrayDst , MDArraySrc , 2 >
{
public:
  typedef KOKKOS_MACRO_DEVICE    device_type ;
  typedef device_type::size_type size_type ;

  MDArrayDst dst ; 
  MDArraySrc src ; 

  MDArrayDeepCopyFunctor( const MDArrayDst & arg_dst ,
                          const MDArraySrc & arg_src )
    : dst( arg_dst ), src( arg_src ) {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < dst.dimension(1) ; ++i1 ) {
      dst(i0,i1) = src(i0,i1);
    }
  }
};

} // namespace Impl
} // namespace Kokkos

#endif /* defined( KOKKOS_MACRO_DEVICE_FUNCTION ) */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  Different maps on the same memory space */
template< typename ValueType , class DeviceSrc , bool Contig >
class MDArrayDeepCopy< ValueType , KOKKOS_MACRO_DEVICE , DeviceSrc ,
                       true   /* Same Memory Space */ ,
                       false  /* Different MDArray Maps */ ,
                       Contig /* Don't care */ >
{
public:
  typedef KOKKOS_MACRO_DEVICE    DeviceDst ;

  typedef MDArrayView<ValueType, DeviceDst > dst_type ;
  typedef MDArrayView<ValueType, DeviceSrc > src_type ;
  
  static void run( const dst_type & dst , const src_type & src )
  {
    switch( dst.rank() ) {
    case 8 :
      parallel_for( dst.dimension(0), Impl::MDArrayDeepCopyFunctor< DeviceDst , dst_type , src_type , 8 >( dst , src ) );
      break ;
    
    case 7 :
      parallel_for( dst.dimension(0), Impl::MDArrayDeepCopyFunctor< DeviceDst , dst_type , src_type , 7 >( dst , src ) );
      break ;

    case 6 :
      parallel_for( dst.dimension(0), Impl::MDArrayDeepCopyFunctor< DeviceDst , dst_type , src_type , 6 >( dst , src ) );
      break ;
    
    case 5 :
      parallel_for( dst.dimension(0), Impl::MDArrayDeepCopyFunctor< DeviceDst , dst_type , src_type , 5 >( dst , src ) );
      break ;

    case 4 :
      parallel_for( dst.dimension(0), Impl::MDArrayDeepCopyFunctor< DeviceDst , dst_type , src_type , 4 >( dst , src ) );
      break ;

    case 3 :
      parallel_for( dst.dimension(0), Impl::MDArrayDeepCopyFunctor< DeviceDst , dst_type , src_type , 3 >( dst , src ) );
      break ;

    case 2 :
      parallel_for( dst.dimension(0), Impl::MDArrayDeepCopyFunctor< DeviceDst , dst_type , src_type , 2 >( dst , src ) );

    case 1 :
      parallel_for( dst.dimension(0), Impl::DeepCopyContiguous< ValueType , DeviceDst >( dst.m_memory.ptr_on_device() , src.m_memory.ptr_on_device() ) );
      break ;

    default: break ;
    }
  }
};

/** \brief  Different maps on the same memory space */
template< typename ValueType , class DeviceDst , bool Contig >
class MDArrayDeepCopy< ValueType , DeviceDst , KOKKOS_MACRO_DEVICE ,
                       true   /* Same Memory Space */ ,
                       false  /* Different MDArray Maps */ ,
                       Contig /* Don't care */ >
{
public:
  typedef KOKKOS_MACRO_DEVICE    DeviceSrc ;

  typedef MDArrayView<ValueType, DeviceDst > dst_type ;
  typedef MDArrayView<ValueType, DeviceSrc > src_type ;
  
  static void run( const dst_type & dst , const src_type & src )
  {
    switch( dst.rank() ) {
    case 8 :
      parallel_for( dst.dimension(0), Impl::MDArrayDeepCopyFunctor< DeviceSrc , dst_type , src_type , 8 >( dst , src ) );
      break ;
    
    case 7 :
      parallel_for( dst.dimension(0), Impl::MDArrayDeepCopyFunctor< DeviceSrc , dst_type , src_type , 7 >( dst , src ) );
      break ;

    case 6 :
      parallel_for( dst.dimension(0), Impl::MDArrayDeepCopyFunctor< DeviceSrc , dst_type , src_type , 6 >( dst , src ) );
      break ;
    
    case 5 :
      parallel_for( dst.dimension(0), Impl::MDArrayDeepCopyFunctor< DeviceSrc , dst_type , src_type , 5 >( dst , src ) );
      break ;

    case 4 :
      parallel_for( dst.dimension(0), Impl::MDArrayDeepCopyFunctor< DeviceSrc , dst_type , src_type , 4 >( dst , src ) );
      break ;

    case 3 :
      parallel_for( dst.dimension(0), Impl::MDArrayDeepCopyFunctor< DeviceSrc , dst_type , src_type , 3 >( dst , src ) );
      break ;

    case 2 :
      parallel_for( dst.dimension(0), Impl::MDArrayDeepCopyFunctor< DeviceSrc , dst_type , src_type , 2 >( dst , src ) );

    case 1 :
      parallel_for( dst.dimension(0), Impl::DeepCopyContiguous< ValueType , DeviceSrc >( dst.m_memory.ptr_on_device() , src.m_memory.ptr_on_device() ) );
      break ;

    default: break ;
    }
  }
};

} // namespace Impl
} // namespace Kokkos


#endif


