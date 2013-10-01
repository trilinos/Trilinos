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

template< class T , class L , class D>
class DualView {
public:

  /* Define base types for Device and Host */

  typedef Kokkos::View<T,L,D> t_dev ;
  typedef typename t_dev::HostMirror t_host ;

  /* Define typedefs for different usage scenarios */

  // Define const view types
  typedef Kokkos::View<typename t_dev::const_data_type,L,D> t_dev_const;
  typedef typename t_dev_const::HostMirror t_host_const;

  // Define const randomread view types
  typedef Kokkos::View<typename t_dev::const_data_type,L,D,Kokkos::MemoryRandomRead> t_dev_const_randomread ;
  typedef typename t_dev_const_randomread::HostMirror t_host_const_randomread;

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

  unsigned int modified_device;
  unsigned int modified_host;

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
    modified_host = 0;
    modified_device = 0;
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
    , h_view( create_mirror_view( d_view ) )
  {
    modified_host = 0;
    modified_device = 0;
  }

  /* Update data on device or host only if other space is polluted */

  template<class Device>
  void sync() {
    unsigned int dev = Kokkos::Impl::if_c< Kokkos::Impl::is_same< typename t_dev::memory_space ,
                                  typename Device::memory_space >::value ,
                                  unsigned int , unsigned int >::select( 1, 0 );

    if(dev) {
      if((modified_host > 0) && (modified_host >= modified_device)) {
      Kokkos::deep_copy(d_view,h_view);
      modified_host = modified_device = 0;
      }
    } else {
      if((modified_device > 0) && (modified_device >= modified_host)) {
      Kokkos::deep_copy(h_view,d_view);
      modified_host = modified_device = 0;
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
      modified_device = (modified_device > modified_host ? modified_device : modified_host)  + 1;
    } else {
      modified_host = (modified_device > modified_host ? modified_device : modified_host)  + 1;
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
     h_view = create_mirror_view( d_view );

     /* Reset dirty flags */
     modified_device = modified_host = 0;
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
   if(modified_device >= modified_host) {
     /* Resize on Device */
     Kokkos::resize(d_view,n0,n1,n2,n3,n4,n5,n6,n7);
     h_view = create_mirror_view( d_view );

     /* Mark Device copy as modified */
     modified_device++;

   } else {
     /* Realloc on Device */

     Kokkos::realloc(d_view,n0,n1,n2,n3,n4,n5,n6,n7);
     t_host temp_view = create_mirror_view( d_view );

     /* Remap on Host */
     Kokkos::Impl::ViewRemap< t_host , t_host >( temp_view , h_view );
     h_view = temp_view;

     /* Mark Host copy as modified */
     modified_host++;
   }
  }

  size_t capacity() const {
    return d_view.capacity();
  }
};
}
#endif
