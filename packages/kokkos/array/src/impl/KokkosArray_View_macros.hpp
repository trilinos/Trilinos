/*
//@HEADER
// ************************************************************************
//
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <iostream>
#include <KokkosArray_Macros.hpp>

#if ! defined(KOKKOSARRAY_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE)                  || \
    ! defined(KOKKOSARRAY_INLINE_FUNCTION)

#error "Including <KokkosArray_View_macros.hpp> without macros defined"

#else

#include <impl/KokkosArray_AnalyzeShape.hpp>

namespace KokkosArray {

//-----------------------------------------------------------------------------
//View -- MemoryManaged
//-----------------------------------------------------------------------------
template< class DataType , class LayoutType >
class View< DataType , LayoutType , KOKKOSARRAY_MACRO_DEVICE, MemoryManaged >
  : public Impl::ViewOper<
      KOKKOSARRAY_MACRO_DEVICE::memory_space ,
      typename Impl::AnalyzeShape< DataType , typename LayoutType::array_layout >::value_type ,
      typename Impl::AnalyzeShape< DataType , typename LayoutType::array_layout >::shape >
{
private:

  typedef Impl::AnalyzeShape<DataType,
                             typename LayoutType::array_layout >  analysis ;
  typedef typename analysis::const_type      const_data_type ;
  typedef typename analysis::non_const_type  non_const_data_type ;

public:

  typedef DataType                  data_type ;
  typedef LayoutType                layout_type ;
  typedef KOKKOSARRAY_MACRO_DEVICE  device_type ;
  typedef MemoryManaged             memory_management_type;

  typedef View<           data_type, layout_type, device_type, memory_management_type > type ;
  typedef View<     const_data_type, layout_type, device_type, memory_management_type > const_type ;
//  typedef View< non_const_data_type, layout_type, device_type > non_const_type ;

  typedef View< data_type , layout_type , Host, memory_management_type >  HostMirror ;

  typedef typename analysis::scalar_type      scalar_type ;
  typedef typename analysis::value_type       value_type ;
  typedef typename LayoutType::array_layout   array_layout ;
  typedef typename device_type::memory_space  memory_space ;
  typedef typename device_type::size_type     size_type ;
  typedef typename analysis::shape            shape_type ;

private:

  typedef Impl::ViewOper< KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                          typename analysis::value_type ,
                          shape_type > oper_type ;

  template< class , class , class , class> friend class View ;

public:

  /*------------------------------------------------------------------*/

  static const unsigned Rank = shape_type::rank ;

  KOKKOSARRAY_INLINE_FUNCTION
  size_type rank() const { return shape_type::rank ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_0() const { return oper_type::m_shape.N0 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_1() const { return oper_type::m_shape.N1 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_2() const { return oper_type::m_shape.N2 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_3() const { return oper_type::m_shape.N3 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_4() const { return oper_type::m_shape.N4 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_5() const { return oper_type::m_shape.N5 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_6() const { return oper_type::m_shape.N6 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_7() const { return oper_type::m_shape.N7 ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension( const iType & r ) const
    {
      return iType(0) == r ? oper_type::m_shape.N0 : (
             iType(1) == r ? oper_type::m_shape.N1 : (
             iType(2) == r ? oper_type::m_shape.N2 : (
             iType(3) == r ? oper_type::m_shape.N3 : (
             iType(4) == r ? oper_type::m_shape.N4 : (
             iType(5) == r ? oper_type::m_shape.N5 : (
             iType(6) == r ? oper_type::m_shape.N6 : (
             iType(7) == r ? oper_type::m_shape.N7 : 0 )))))));
    }

  /*------------------------------------------------------------------*/
  KOKKOSARRAY_INLINE_FUNCTION
  scalar_type * ptr_on_device() const { return oper_type::m_ptr_on_device ; }

  /*------------------------------------------------------------------*/
  /** \brief Query shape */
  KOKKOSARRAY_INLINE_FUNCTION
  shape_type shape() const { return oper_type::m_shape ; }

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  KOKKOSARRAY_INLINE_FUNCTION
  operator bool () const
  { return 0 != oper_type::m_ptr_on_device ; }

  /** \brief  Query if view to same memory */
  KOKKOSARRAY_INLINE_FUNCTION
  bool operator == ( const View & rhs ) const
  {
    return oper_type::m_ptr_on_device == rhs.m_ptr_on_device &&
           oper_type::m_shape == rhs.m_shape ;
  }

  /** \brief  Query if not view to same memory */
  KOKKOSARRAY_INLINE_FUNCTION
  bool operator != ( const View & rhs ) const
  {
    return oper_type::m_ptr_on_device != rhs.m_ptr_on_device ||
           oper_type::m_shape != rhs.m_shape ;
  }

  /** \brief  Query if view to same memory */
  template< class rhsDataType , class rhsLayout , class rhsMemory, class rhsMemManagement >
  KOKKOSARRAY_INLINE_FUNCTION
  bool operator == ( const View<rhsDataType,rhsLayout,rhsMemory,rhsMemManagement> & rhs ) const
  {
    return oper_type::m_ptr_on_device == rhs.m_ptr_on_device &&
           oper_type::m_shape == rhs.m_shape ;
  }

  /** \brief  Query if not view to same memory */
  template< class rhsDataType , class rhsLayout , class rhsMemory, class rhsMemManagement >
  KOKKOSARRAY_INLINE_FUNCTION
  bool operator != ( const View<rhsDataType,rhsLayout,rhsMemory,rhsMemManagement> & rhs ) const
  {
    return oper_type::m_ptr_on_device != rhs.m_ptr_on_device ||
           oper_type::m_shape != rhs.m_shape ;
  }

  /*------------------------------------------------------------------*/

private:

  KOKKOSARRAY_INLINE_FUNCTION
  void internal_private_assign( const shape_type & shape , scalar_type * ptr )
  {
    oper_type::m_shape          = shape ;
    oper_type::m_ptr_on_device  = ptr ;
    memory_space::increment( oper_type::m_ptr_on_device );
  }

  KOKKOSARRAY_INLINE_FUNCTION
  void internal_private_clear()
  {
    memory_space::decrement( oper_type::m_ptr_on_device );
    oper_type::m_ptr_on_device = 0 ;
  }

  /** \brief  This device-specialized method may only be called
   *          within a view constructor.
   */
  void internal_private_create( const std::string & label ,
                                const shape_type shape );

public:

  /** \brief  Construct a NULL view */
  KOKKOSARRAY_INLINE_FUNCTION
  View()
    {
      oper_type::m_ptr_on_device = 0 ;
      oper_type::m_shape = shape_type();
    }

  /** \brief  Construct a view of the array */
  KOKKOSARRAY_INLINE_FUNCTION
  View( const View & rhs )
    {
      internal_private_assign( rhs.m_shape , rhs.m_ptr_on_device );
    }

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  KOKKOSARRAY_INLINE_FUNCTION
  ~View()
  {
    internal_private_clear();
  }

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  KOKKOSARRAY_INLINE_FUNCTION
  View & operator = ( const View & rhs )
  {
    internal_private_clear();
    internal_private_assign( rhs.m_shape , rhs.m_ptr_on_device );
    return *this ;
  }

  /*------------------------------------------------------------------*/
  // Construct a compatible view:

  template< class rhsType , class rhsLayout , class rhsMemory, class rhsMemManagement >
  View( const View< rhsType , rhsLayout , rhsMemory, rhsMemManagement > & rhs )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type  rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type data( rhs.shape() );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  // Assign a compatible view:

  template< class rhsType , class rhsLayout , class rhsMemory, class rhsMemManagement >
  View & operator = ( const View< rhsType , rhsLayout , rhsMemory, rhsMemManagement > & rhs )
    {
      typedef View< rhsType , rhsLayout , rhsMemory, rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type data( rhs.shape() );

      internal_private_clear();
      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );

      return *this ;
    }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a subview */

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement, class ArgType0 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > & rhs ,
        const ArgType0 & arg0 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type data( rhs.shape() , arg0 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement ,
            class ArgType0 , class ArgType1 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement> & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement> rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type data( rhs.shape() , arg0 , arg1 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement ,
            class ArgType0 , class ArgType1 , class ArgType2 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type data( rhs.shape() , arg0 , arg1 , arg2 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement ,
            class ArgType0 , class ArgType1 , class ArgType2 ,
            class ArgType3 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement> & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 ,
        const ArgType3 & arg3 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type data( rhs.shape() , arg0 , arg1 , arg2 , arg3 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement ,
            class ArgType0 , class ArgType1 , class ArgType2 ,
            class ArgType3 , class ArgType4 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement> & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 ,
        const ArgType3 & arg3 ,
        const ArgType4 & arg4 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type data( rhs.shape(), arg0, arg1, arg2, arg3, arg4 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement ,
            class ArgType0 , class ArgType1 , class ArgType2 ,
            class ArgType3 , class ArgType4 , class ArgType5 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement> & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 ,
        const ArgType3 & arg3 ,
        const ArgType4 & arg4 ,
        const ArgType5 & arg5 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type
         data( rhs.shape() , arg0 , arg1 , arg2 , arg3 , arg4 , arg5 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement ,
            class ArgType0 , class ArgType1 , class ArgType2 ,
            class ArgType3 , class ArgType4 , class ArgType5 ,
            class ArgType6 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement> & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 ,
        const ArgType3 & arg3 ,
        const ArgType4 & arg4 ,
        const ArgType5 & arg5 ,
        const ArgType6 & arg6 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type
         data( rhs.shape() , arg0 , arg1 , arg2 , arg3 , arg4 , arg5 , arg6 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement ,
            class ArgType0 , class ArgType1 , class ArgType2 ,
            class ArgType3 , class ArgType4 , class ArgType5 ,
            class ArgType6 , class ArgType7 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement> & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 ,
        const ArgType3 & arg3 ,
        const ArgType4 & arg4 ,
        const ArgType5 & arg5 ,
        const ArgType6 & arg6 ,
        const ArgType7 & arg7 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type
         data( rhs.shape(), arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/
  /* Creation with allocation of memory on the device */

  View( const std::string & label , const shape_type shape )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( label , shape );
  }

  explicit View( const std::string & label )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>() );
  }

  View( const std::string & label ,
        const size_t n0 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0) );
  }

  View( const std::string & label ,
        const size_t n0 ,
        const size_t n1 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0,n1) );
  }

  View( const std::string & label ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0,n1,n2) );
  }

  View( const std::string & label ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 ,
        const size_t n3 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0,n1,n2,n3) );
  }

  View( const std::string & label ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 ,
        const size_t n3 ,
        const size_t n4 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0,n1,n2,n3,n4) );
  }

  View( const std::string & label ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 ,
        const size_t n3 ,
        const size_t n4 ,
        const size_t n5 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0,n1,n2,n3,n4,n5) );
  }

  View( const std::string & label ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 ,
        const size_t n3 ,
        const size_t n4 ,
        const size_t n5 ,
        const size_t n6 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0,n1,n2,n3,n4,n5,n6) );
  }

  View( const std::string & label ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 ,
        const size_t n3 ,
        const size_t n4 ,
        const size_t n5 ,
        const size_t n6 ,
        const size_t n7 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0,n1,n2,n3,n4,n5,n6,n7) );
  }

  /*------------------------------------------------------------------*/
  /** \brief  For unit testing shape mapping. */
  explicit View( const shape_type & shape )
    { internal_private_assign( shape , 0 ); }
};

//----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//View -- MemoryUnmanaged
//-----------------------------------------------------------------------------
template< class DataType , class LayoutType >
class View< DataType , LayoutType , KOKKOSARRAY_MACRO_DEVICE, MemoryUnmanaged >
  : public Impl::ViewOper<
      KOKKOSARRAY_MACRO_DEVICE::memory_space ,
      typename Impl::AnalyzeShape< DataType , typename LayoutType::array_layout >::value_type ,
      typename Impl::AnalyzeShape< DataType , typename LayoutType::array_layout >::shape >
{
private:

  typedef Impl::AnalyzeShape<DataType,
                             typename LayoutType::array_layout >  analysis ;
  typedef typename analysis::const_type      const_data_type ;
  typedef typename analysis::non_const_type  non_const_data_type ;

public:

  typedef DataType                  data_type ;
  typedef LayoutType                layout_type ;
  typedef KOKKOSARRAY_MACRO_DEVICE  device_type ;
  typedef MemoryUnmanaged           memory_management_type;

  typedef View<           data_type, layout_type, device_type, memory_management_type > type ;
  typedef View<     const_data_type, layout_type, device_type, memory_management_type > const_type ;
//  typedef View< non_const_data_type, layout_type, device_type > non_const_type ;

  typedef View< data_type , layout_type , Host, memory_management_type >  HostMirror ;

  typedef typename analysis::scalar_type      scalar_type ;
  typedef typename analysis::value_type       value_type ;
  typedef typename LayoutType::array_layout   array_layout ;
  typedef typename device_type::memory_space  memory_space ;
  typedef typename device_type::size_type     size_type ;
  typedef typename analysis::shape            shape_type ;

private:

  typedef Impl::ViewOper< KOKKOSARRAY_MACRO_DEVICE::memory_space ,
                          typename analysis::value_type ,
                          shape_type > oper_type ;

  template< class , class , class , class> friend class View ;

public:

  /*------------------------------------------------------------------*/

  static const unsigned Rank = shape_type::rank ;

  KOKKOSARRAY_INLINE_FUNCTION
  size_type rank() const { return shape_type::rank ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_0() const { return oper_type::m_shape.N0 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_1() const { return oper_type::m_shape.N1 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_2() const { return oper_type::m_shape.N2 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_3() const { return oper_type::m_shape.N3 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_4() const { return oper_type::m_shape.N4 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_5() const { return oper_type::m_shape.N5 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_6() const { return oper_type::m_shape.N6 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension_7() const { return oper_type::m_shape.N7 ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension( const iType & r ) const
    {
      return iType(0) == r ? oper_type::m_shape.N0 : (
             iType(1) == r ? oper_type::m_shape.N1 : (
             iType(2) == r ? oper_type::m_shape.N2 : (
             iType(3) == r ? oper_type::m_shape.N3 : (
             iType(4) == r ? oper_type::m_shape.N4 : (
             iType(5) == r ? oper_type::m_shape.N5 : (
             iType(6) == r ? oper_type::m_shape.N6 : (
             iType(7) == r ? oper_type::m_shape.N7 : 0 )))))));
    }

  /*------------------------------------------------------------------*/
  KOKKOSARRAY_INLINE_FUNCTION
  scalar_type * ptr_on_device() const { return oper_type::m_ptr_on_device ; }

  /*------------------------------------------------------------------*/
  /** \brief Query shape */
  KOKKOSARRAY_INLINE_FUNCTION
  shape_type shape() const { return oper_type::m_shape ; }

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  KOKKOSARRAY_INLINE_FUNCTION
  operator bool () const
  { return 0 != oper_type::m_ptr_on_device ; }

  /** \brief  Query if view to same memory */
  KOKKOSARRAY_INLINE_FUNCTION
  bool operator == ( const View & rhs ) const
  {
    return oper_type::m_ptr_on_device == rhs.m_ptr_on_device &&
           oper_type::m_shape == rhs.m_shape ;
  }

  /** \brief  Query if not view to same memory */
  KOKKOSARRAY_INLINE_FUNCTION
  bool operator != ( const View & rhs ) const
  {
    return oper_type::m_ptr_on_device != rhs.m_ptr_on_device ||
           oper_type::m_shape != rhs.m_shape ;
  }

  /** \brief  Query if view to same memory */
  template< class rhsDataType , class rhsLayout , class rhsMemory, class rhsMemManagement >
  KOKKOSARRAY_INLINE_FUNCTION
  bool operator == ( const View<rhsDataType,rhsLayout,rhsMemory,rhsMemManagement> & rhs ) const
  {
    return oper_type::m_ptr_on_device == rhs.m_ptr_on_device &&
           oper_type::m_shape == rhs.m_shape ;
  }

  /** \brief  Query if not view to same memory */
  template< class rhsDataType , class rhsLayout , class rhsMemory, class rhsMemManagement >
  KOKKOSARRAY_INLINE_FUNCTION
  bool operator != ( const View<rhsDataType,rhsLayout,rhsMemory,rhsMemManagement> & rhs ) const
  {
    return oper_type::m_ptr_on_device != rhs.m_ptr_on_device ||
           oper_type::m_shape != rhs.m_shape ;
  }

  /*------------------------------------------------------------------*/

private:

  KOKKOSARRAY_INLINE_FUNCTION
  void internal_private_assign( const shape_type & shape , scalar_type * ptr )
  {
    oper_type::m_shape          = shape ;
    oper_type::m_ptr_on_device  = ptr ;
  }

  KOKKOSARRAY_INLINE_FUNCTION
  void internal_private_clear()
  {
    oper_type::m_ptr_on_device = 0 ;
  }

  KOKKOSARRAY_INLINE_FUNCTION
  void internal_private_create( scalar_type * data_ptr ,
                                const shape_type shape )
  {
    oper_type::m_ptr_on_device = data_ptr;
    oper_type::m_shape = shape ;
  }

public:

  /** \brief  Construct a NULL view */
  KOKKOSARRAY_INLINE_FUNCTION
  View()
    {
      oper_type::m_ptr_on_device = 0 ;
      oper_type::m_shape = shape_type();
    }

  /** \brief  Construct a view of the array */
  KOKKOSARRAY_INLINE_FUNCTION
  View( const View & rhs )
    {
      internal_private_assign( rhs.m_shape , rhs.m_ptr_on_device );
    }

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  KOKKOSARRAY_INLINE_FUNCTION
  ~View()
  {
    internal_private_clear();
  }

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  KOKKOSARRAY_INLINE_FUNCTION
  View & operator = ( const View & rhs )
  {
    internal_private_clear();
    internal_private_assign( rhs.m_shape , rhs.m_ptr_on_device );
    return *this ;
  }

  /*------------------------------------------------------------------*/
  // Construct a compatible view:

  template< class rhsType , class rhsLayout , class rhsMemory, class rhsMemManagement >
  View( const View< rhsType , rhsLayout , rhsMemory, rhsMemManagement > & rhs )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type data( rhs.shape() );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  // Assign a compatible view:

  template< class rhsType , class rhsLayout , class rhsMemory, class rhsMemManagement >
  View & operator = ( const View< rhsType , rhsLayout , rhsMemory, rhsMemManagement > & rhs )
    {
      typedef View< rhsType , rhsLayout , rhsMemory, rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type data( rhs.shape() );

      internal_private_clear();
      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );

      return *this ;
    }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a subview */

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement, class ArgType0 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > & rhs ,
        const ArgType0 & arg0 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type data( rhs.shape() , arg0 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement ,
            class ArgType0 , class ArgType1 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement> & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement> rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type data( rhs.shape() , arg0 , arg1 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement ,
            class ArgType0 , class ArgType1 , class ArgType2 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type data( rhs.shape() , arg0 , arg1 , arg2 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement ,
            class ArgType0 , class ArgType1 , class ArgType2 ,
            class ArgType3 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement> & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 ,
        const ArgType3 & arg3 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type data( rhs.shape() , arg0 , arg1 , arg2 , arg3 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement ,
            class ArgType0 , class ArgType1 , class ArgType2 ,
            class ArgType3 , class ArgType4 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement> & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 ,
        const ArgType3 & arg3 ,
        const ArgType4 & arg4 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type data( rhs.shape(), arg0, arg1, arg2, arg3, arg4 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement ,
            class ArgType0 , class ArgType1 , class ArgType2 ,
            class ArgType3 , class ArgType4 , class ArgType5 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement> & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 ,
        const ArgType3 & arg3 ,
        const ArgType4 & arg4 ,
        const ArgType5 & arg5 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type
         data( rhs.shape() , arg0 , arg1 , arg2 , arg3 , arg4 , arg5 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement ,
            class ArgType0 , class ArgType1 , class ArgType2 ,
            class ArgType3 , class ArgType4 , class ArgType5 ,
            class ArgType6 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement> & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 ,
        const ArgType3 & arg3 ,
        const ArgType4 & arg4 ,
        const ArgType5 & arg5 ,
        const ArgType6 & arg6 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type
         data( rhs.shape() , arg0 , arg1 , arg2 , arg3 , arg4 , arg5 , arg6 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory , class rhsMemManagement ,
            class ArgType0 , class ArgType1 , class ArgType2 ,
            class ArgType3 , class ArgType4 , class ArgType5 ,
            class ArgType6 , class ArgType7 >
  View( const View< rhsType , rhsLayout , rhsMemory , rhsMemManagement> & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 ,
        const ArgType3 & arg3 ,
        const ArgType4 & arg4 ,
        const ArgType5 & arg5 ,
        const ArgType6 & arg6 ,
        const ArgType7 & arg7 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory , rhsMemManagement > rhs_view ;
      typedef typename rhs_view::scalar_type   rhs_scalar_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< scalar_type , memory_space ,
                                 rhs_scalar_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef typename
        Impl::SubShape< shape_type , rhs_shape_type >::type
          subshape_type ;

      const subshape_type
         data( rhs.shape(), arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/
  /* Creation with allocation of memory on the device */

  View( scalar_type * data_ptr , const shape_type shape )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( data_ptr , shape );
  }

  explicit View( scalar_type * data_ptr )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( data_ptr , shape_type::template create<memory_space>() );
  }

  View( scalar_type * data_ptr ,
        const size_t n0 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( data_ptr , shape_type::template create<memory_space>(n0) );
  }

  View( scalar_type * data_ptr ,
        const size_t n0 ,
        const size_t n1 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( data_ptr , shape_type::template create<memory_space>(n0,n1) );
  }

  View( scalar_type * data_ptr ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( data_ptr , shape_type::template create<memory_space>(n0,n1,n2) );
  }

  View( scalar_type * data_ptr ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 ,
        const size_t n3 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( data_ptr , shape_type::template create<memory_space>(n0,n1,n2,n3) );
  }

  View( scalar_type * data_ptr ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 ,
        const size_t n3 ,
        const size_t n4 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( data_ptr , shape_type::template create<memory_space>(n0,n1,n2,n3,n4) );
  }

  View( scalar_type * data_ptr ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 ,
        const size_t n3 ,
        const size_t n4 ,
        const size_t n5 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( data_ptr , shape_type::template create<memory_space>(n0,n1,n2,n3,n4,n5) );
  }

  View( scalar_type * data_ptr ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 ,
        const size_t n3 ,
        const size_t n4 ,
        const size_t n5 ,
        const size_t n6 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( data_ptr , shape_type::template create<memory_space>(n0,n1,n2,n3,n4,n5,n6) );
  }

  View( scalar_type * data_ptr ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 ,
        const size_t n3 ,
        const size_t n4 ,
        const size_t n5 ,
        const size_t n6 ,
        const size_t n7 )
  {
    View_create_requires_non_const_data_type< scalar_type >::success();
    internal_private_create( data_ptr , shape_type::template create<memory_space>(n0,n1,n2,n3,n4,n5,n6,n7) );
  }

  /*------------------------------------------------------------------*/
  /** \brief  For unit testing shape mapping. */
  explicit View( const shape_type & shape )
    { internal_private_assign( shape , 0 ); }
};

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

