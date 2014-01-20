/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/


/* DualView: container class to manage data structures which exist both on Host and Device
 * member functions:
 * DualView()
 * DualView(label,dim0,dim1,dim2,...)
 * view<DeviceType>()
 * sync<DeviceType>()
 * modify<DeviceType>()
 * resize(dim0,dim1,dim2,...)
 */
#ifndef KOKKOS_DUALVIEW_HPP
#define KOKKOS_DUALVIEW_HPP

#include <Kokkos_View.hpp>
namespace Kokkos {

template< class T , class L , class D, class M = MemoryManaged>
class DualView {
public:

  typedef D device_type;
  typedef typename D::host_mirror_device_type host_mirror_device_type;

  /* Define base types for Device and Host */

  typedef Kokkos::View<T,L,D,M> t_dev ;
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_UVM)
  typedef t_dev t_host;
#else
  typedef typename t_dev::HostMirror t_host ;
#endif
  /* Define typedefs for different usage scenarios */

  // Define const view types
  typedef Kokkos::View<typename t_dev::const_data_type,L,D,M> t_dev_const;
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_UVM)
  typedef t_dev_const t_host_const;
#else
  typedef typename t_dev_const::HostMirror t_host_const;
#endif
  // Define const randomread view types
  typedef Kokkos::View<typename t_dev::const_data_type,L,D,Kokkos::MemoryRandomAccess> t_dev_const_randomread ;
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_UVM)
  typedef t_dev_const_randomread t_host_const_randomread;
#else
  typedef typename t_dev_const_randomread::HostMirror t_host_const_randomread;
#endif

  // Define unmanaged view types
  typedef Kokkos::View<T,L,D,Kokkos::MemoryUnmanaged> t_dev_um;
  typedef Kokkos::View<typename t_host::data_type,typename t_host::array_layout,
                       typename t_host::device_type,Kokkos::MemoryUnmanaged> t_host_um;

  // Define const unmanaged view types
  typedef Kokkos::View<typename t_dev::const_data_type,L,D,Kokkos::MemoryUnmanaged> t_dev_const_um;
  typedef Kokkos::View<typename t_host::const_data_type,typename t_host::array_layout,
                       typename t_host::device_type,Kokkos::MemoryUnmanaged> t_host_const_um;

  /* provide the same typedefs as a view for scalar, data and value types */

  typedef typename t_dev::value_type value_type;
  typedef typename t_dev::const_value_type const_value_type;
  typedef typename t_dev::scalar_type scalar_type;
  typedef typename t_dev::const_scalar_type const_scalar_type;
  typedef typename t_dev::non_const_scalar_type non_const_scalar_type;

  /* Instances of base types */

  t_dev d_view;
  t_host h_view;


  /* Counters to keep track of changes (dirty-flags) */

  View<unsigned int,LayoutLeft,host_mirror_device_type> modified_device;
  View<unsigned int,LayoutLeft,host_mirror_device_type> modified_host;

  /* Return view on specific device via view<Device>() */

  template< class Device >
  const typename Kokkos::Impl::if_c< Kokkos::Impl::is_same< typename t_dev::memory_space ,
                                      typename Device::memory_space >::value ,
                             t_dev , t_host >::type view() const
  {
    return Kokkos::Impl::if_c< Kokkos::Impl::is_same< typename t_dev::memory_space ,
                                typename Device::memory_space >::value ,
                       t_dev , t_host >::select( d_view , h_view );
  }


  /* Construct views */

  /* Empty Constructor */

  DualView() {
    modified_host = View<unsigned int,LayoutLeft,host_mirror_device_type>("DualView::modified_host");
    modified_device = View<unsigned int,LayoutLeft,host_mirror_device_type>("DualView::modified_device");
  }

  /* Create view with allocation on both host and device */

  DualView( const std::string & label ,
    const size_t n0 = 0 ,
    const size_t n1 = 0 ,
    const size_t n2 = 0 ,
    const size_t n3 = 0 ,
    const size_t n4 = 0 ,
    const size_t n5 = 0 ,
    const size_t n6 = 0 ,
    const size_t n7 = 0 )
    : d_view( label, n0, n1, n2, n3, n4, n5, n6, n7 )
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_UVM)
    , h_view( d_view )
#else
    , h_view( create_mirror_view( d_view ) )
#endif
  {

    modified_host = View<unsigned int,LayoutLeft,host_mirror_device_type>("DualView::modified_host");;
    modified_device = View<unsigned int,LayoutLeft,host_mirror_device_type>("DualView::modified_device");;
  }

  /* Update data on device or host only if other space is polluted */

  template<class SS, class LS, class DS, class MS>
  DualView(const DualView<SS,LS,DS,MS> src) {
    d_view = src.d_view;
    h_view = src.h_view;
    modified_host = src.modified_host;
    modified_device = src.modified_device;
  }

  template<class Device>
  void sync() {
    unsigned int dev = Kokkos::Impl::if_c< Kokkos::Impl::is_same< typename t_dev::memory_space ,
                                  typename Device::memory_space >::value ,
                                  unsigned int , unsigned int >::select( 1, 0 );

    if(dev) {
      if((modified_host() > 0) && (modified_host() >= modified_device())) {
      Kokkos::deep_copy(d_view,h_view);
      modified_host() = modified_device() = 0;
      }
    } else {
      if((modified_device() > 0) && (modified_device() >= modified_host())) {
      Kokkos::deep_copy(h_view,d_view);
      modified_host() = modified_device() = 0;
      }
    }
  }

  /* Mark data as dirty on a device */

  template<class Device>
  void modify() {
    unsigned int dev = Kokkos::Impl::if_c< Kokkos::Impl::is_same< typename t_dev::memory_space ,
                                  typename Device::memory_space >::value ,
                                  unsigned int , unsigned int >::select( 1, 0 );

    if(dev) {
      modified_device() = (modified_device() > modified_host() ? modified_device() : modified_host())  + 1;
    } else {
      modified_host() = (modified_device() > modified_host() ? modified_device() : modified_host())  + 1;
    }
  }

  /* Realloc both views, no deep copy */

  void realloc( const size_t n0 = 0 ,
           const size_t n1 = 0 ,
           const size_t n2 = 0 ,
           const size_t n3 = 0 ,
           const size_t n4 = 0 ,
           const size_t n5 = 0 ,
           const size_t n6 = 0 ,
           const size_t n7 = 0 ) {
     Kokkos::realloc(d_view,n0,n1,n2,n3,n4,n5,n6,n7);
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_UVM)
     h_view = d_view ;
#else
     h_view = create_mirror_view( d_view );
#endif
     /* Reset dirty flags */
     modified_device() = modified_host() = 0;
  }

  /* Resize both views, only do deep_copy in space which was last marked as dirty */

  void resize( const size_t n0 = 0 ,
           const size_t n1 = 0 ,
           const size_t n2 = 0 ,
           const size_t n3 = 0 ,
           const size_t n4 = 0 ,
           const size_t n5 = 0 ,
           const size_t n6 = 0 ,
           const size_t n7 = 0 ) {
   if(modified_device() >= modified_host()) {
     /* Resize on Device */
     Kokkos::resize(d_view,n0,n1,n2,n3,n4,n5,n6,n7);
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_UVM)
     h_view = d_view ;
#else
     h_view = create_mirror_view( d_view );
#endif

     /* Mark Device copy as modified */
     modified_device() = modified_device()+1;

   } else {
     /* Realloc on Device */

     Kokkos::realloc(d_view,n0,n1,n2,n3,n4,n5,n6,n7);
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_UVM)
     t_host temp_view = d_view ;
#else
     t_host temp_view = create_mirror_view( d_view );
#endif

     /* Remap on Host */
     Kokkos::Impl::ViewRemap< t_host , t_host >( temp_view , h_view );
     h_view = temp_view;

     /* Mark Host copy as modified */
     modified_host() = modified_host()+1;
   }
  }

  size_t capacity() const {
    return d_view.capacity();
  }

  template< typename iType>
  void stride(iType* stride_) {
    d_view.stride(stride_);
  }

  size_t dimension_0() const {return d_view.dimension_0();}
  size_t dimension_1() const {return d_view.dimension_1();}
  size_t dimension_2() const {return d_view.dimension_2();}
  size_t dimension_3() const {return d_view.dimension_3();}
  size_t dimension_4() const {return d_view.dimension_4();}
  size_t dimension_5() const {return d_view.dimension_5();}
  size_t dimension_6() const {return d_view.dimension_6();}
  size_t dimension_7() const {return d_view.dimension_7();}
};

template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}


template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 , class ArgType1 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0,arg1);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0,arg1);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}

template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 , class ArgType1 , class ArgType2 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0,arg1,arg2);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0,arg1,arg2);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}

template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0,arg1,arg2,arg3);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0,arg1,arg2,arg3);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}

template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 ,
         const ArgType4 & arg4 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0,arg1,arg2,arg3,arg4);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0,arg1,arg2,arg3,arg4);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}

template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 , class ArgType5 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 ,
         const ArgType4 & arg4 ,
         const ArgType5 & arg5 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0,arg1,arg2,arg3,arg4,arg5);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0,arg1,arg2,arg3,arg4,arg5);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}

template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 , class ArgType5 , class ArgType6 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 ,
         const ArgType4 & arg4 ,
         const ArgType5 & arg5 ,
         const ArgType6 & arg6 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0,arg1,arg2,arg3,arg4,arg5,arg6);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0,arg1,arg2,arg3,arg4,arg5,arg6);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}

template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 , class ArgType5 , class ArgType6 , class ArgType7 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 ,
         const ArgType4 & arg4 ,
         const ArgType5 & arg5 ,
         const ArgType6 & arg6 ,
         const ArgType7 & arg7 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}

template< class DT , class DL , class DD , class DM ,
          class ST , class SL , class SD , class SM >
void deep_copy( DualView<DT,DL,DD,DM> dst, const DualView<ST,SL,SD,SM> & src) {
  if(src.modified_device() >= src.modified_host()) {
    Kokkos::deep_copy(dst.d_view,src.d_view);
    dst.template modify<typename DualView<DT,DL,DD,DM>::device_type >();
  } else {
    Kokkos::deep_copy(dst.h_view,src.h_view);
    dst.template modify<typename DualView<DT,DL,DD,DM>::host_mirror_device_type >();
  }
}
}
#endif
