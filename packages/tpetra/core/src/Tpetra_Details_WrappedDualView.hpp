// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_WRAPPEDDUALVIEW_HPP
#define TPETRA_DETAILS_WRAPPEDDUALVIEW_HPP

#include <Tpetra_Access.hpp>
#include <Kokkos_DualView.hpp>
#include "Teuchos_TestForException.hpp"
#include <sstream>

//#define DEBUG_UVM_REMOVAL  // Works only with gcc > 4.8

#ifdef DEBUG_UVM_REMOVAL

#define DEBUG_UVM_REMOVAL_ARGUMENT ,const char* callerstr = __builtin_FUNCTION()

#define DEBUG_UVM_REMOVAL_PRINT_CALLER(fn) \
  { \
  auto envVarSet = std::getenv("TPETRA_UVM_REMOVAL"); \
  if (envVarSet && (std::strcmp(envVarSet,"1") == 0)) \
    std::cout << (fn) << " called from " << callerstr \
              << " host cnt " << dualView.h_view.use_count()  \
              << " device cnt " << dualView.d_view.use_count()  \
              << std::endl; \
  }

#else

#define DEBUG_UVM_REMOVAL_ARGUMENT
#define DEBUG_UVM_REMOVAL_PRINT_CALLER(fn)

#endif

//! Namespace for Tpetra classes and methods
namespace Tpetra {

/// \brief Namespace for Tpetra implementation details.
/// \warning Do NOT rely on the contents of this namespace.
namespace Details {

namespace impl {

template <typename DualViewType>
struct hasConstData {
  using valueType = typename DualViewType::value_type;
  using constValueType = typename DualViewType::const_value_type;
  static constexpr bool value = std::is_same<valueType, constValueType>::value;
};

template <typename DualViewType>
using enableIfConstData = std::enable_if_t<hasConstData<DualViewType>::value>;

template <typename DualViewType>
using enableIfNonConstData = std::enable_if_t<!hasConstData<DualViewType>::value>;

template <typename DualViewType>
enableIfNonConstData<DualViewType>
sync_host(DualViewType dualView) {
  dualView.sync_host();
}

template <typename DualViewType>
enableIfConstData<DualViewType>
sync_host(DualViewType dualView) { }

template <typename DualViewType>
enableIfNonConstData<DualViewType>
sync_device(DualViewType dualView) {
  dualView.sync_device();
}

template <typename DualViewType>
enableIfConstData<DualViewType>
sync_device(DualViewType dualView) { }

}


template <typename DualViewType>
class WrappedDualView {
public:
  using HostViewType = typename DualViewType::t_host;
  using DeviceViewType = typename DualViewType::t_dev;

  using HostType   = typename HostViewType::device_type;
  using DeviceType = typename DeviceViewType::device_type;

  using DVT = DualViewType;
  using t_host = typename DualViewType::t_host;
  using t_dev  = typename DualViewType::t_dev;

private:
  static constexpr bool dualViewHasNonConstData = !impl::hasConstData<DualViewType>::value;
  static constexpr bool deviceMemoryIsHostAccessible =
    Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace, typename DeviceViewType::memory_space>::accessible;

public:
  WrappedDualView() {}

  WrappedDualView(DualViewType dualV)
    : originalDualView(dualV),
      dualView(originalDualView)
  { }

  // Should this be expert only?
  WrappedDualView(DualViewType dualV,DualViewType origDualV)
    : originalDualView(origDualV),
      dualView(dualV)
  { }


  WrappedDualView(const DeviceViewType deviceView) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        deviceView.data() != nullptr && deviceView.use_count() == 0,
        std::invalid_argument,
        "Tpetra::Details::WrappedDualView: cannot construct with a device view that\n"
        "does not own its memory (i.e. constructed with a raw pointer and dimensions)\n"
        "because the WrappedDualView needs to assume ownership of the memory.");
    //If the provided view is default-constructed (null, 0 extent, 0 use count),
    //leave the host mirror default-constructed as well in order to have a matching use count of 0.
    HostViewType hostView;
    if(deviceView.use_count() != 0)
    {
      hostView = Kokkos::create_mirror_view_and_copy(
          typename HostViewType::memory_space(),
          deviceView);
    }
    originalDualView = DualViewType(deviceView, hostView);
    dualView = originalDualView;
  }

  WrappedDualView(const WrappedDualView parent, int offset, int numEntries) {
    originalDualView = parent.originalDualView;
    dualView = getSubview(parent.dualView, offset, numEntries);
  }

  WrappedDualView(const WrappedDualView parent,const Kokkos::pair<size_t,size_t>& rowRng, const Kokkos::Impl::ALL_t& colRng) {
    originalDualView = parent.originalDualView;
    dualView = getSubview2D(parent.dualView,rowRng,colRng);
  }

  WrappedDualView(const WrappedDualView parent,const Kokkos::Impl::ALL_t &rowRng, const Kokkos::pair<size_t,size_t>& colRng) {
    originalDualView = parent.originalDualView;
    dualView = getSubview2D(parent.dualView,rowRng,colRng);
  }

  WrappedDualView(const WrappedDualView parent,const Kokkos::pair<size_t,size_t>& rowRng, const Kokkos::pair<size_t,size_t>& colRng) {
    originalDualView = parent.originalDualView;
    dualView = getSubview2D(parent.dualView,rowRng,colRng);
  }


  size_t extent(const int i) const {
    return dualView.h_view.extent(i);
  }
  

  void stride(size_t * stride_) const {
    dualView.stride(stride_);
  }


  typename HostViewType::const_type
  getHostView(Access::ReadOnlyStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  ) const 
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getHostViewReadOnly");
    if(needsSyncPath()) {
      throwIfDeviceViewAlive();
      impl::sync_host(originalDualView);
    }
    return dualView.view_host();
  }

  HostViewType
  getHostView(Access::ReadWriteStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  ) 
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getHostViewReadWrite");
    static_assert(dualViewHasNonConstData,
        "ReadWrite views are not available for DualView with const data");
    if(needsSyncPath()) {
      throwIfDeviceViewAlive();
      impl::sync_host(originalDualView);
      originalDualView.modify_host();
    }

    return dualView.view_host();
  }

  HostViewType
  getHostView(Access::OverwriteAllStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  ) 
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getHostViewOverwriteAll");
    static_assert(dualViewHasNonConstData,
        "OverwriteAll views are not available for DualView with const data");
    if (iAmASubview()) {
      return getHostView(Access::ReadWrite);
    }
    if(needsSyncPath()) {
      if (deviceMemoryIsHostAccessible) Kokkos::fence();      
      throwIfDeviceViewAlive();
      dualView.clear_sync_state();
      dualView.modify_host();
    }
    return dualView.view_host();
  }

  typename DeviceViewType::const_type
  getDeviceView(Access::ReadOnlyStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  ) const 
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getDeviceViewReadOnly");
    if(needsSyncPath()) {
      throwIfHostViewAlive();
      impl::sync_device(originalDualView);
    }
    return dualView.view_device();
  }

  DeviceViewType
  getDeviceView(Access::ReadWriteStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  ) 
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getDeviceViewReadWrite");
    static_assert(dualViewHasNonConstData,
        "ReadWrite views are not available for DualView with const data");
    if(needsSyncPath()) {
      throwIfHostViewAlive();
      impl::sync_device(originalDualView);
      originalDualView.modify_device();
    }
    return dualView.view_device();
  }

  DeviceViewType
  getDeviceView(Access::OverwriteAllStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  ) 
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getDeviceViewOverwriteAll");
    static_assert(dualViewHasNonConstData,
        "OverwriteAll views are not available for DualView with const data");
    if (iAmASubview()) {
      return getDeviceView(Access::ReadWrite);
    }
    if(needsSyncPath()) {
      if (deviceMemoryIsHostAccessible) Kokkos::fence();
      throwIfHostViewAlive();
      dualView.clear_sync_state();
      dualView.modify_device();
    }
    return dualView.view_device();
  }

  template<class TargetDeviceType>
  typename std::remove_reference<decltype(std::declval<DualViewType>().template view<TargetDeviceType>())>::type::const_type
  getView (Access::ReadOnlyStruct s DEBUG_UVM_REMOVAL_ARGUMENT) const {    
    bool returnDevice = true;
    {
      auto tmp = dualView.template view<TargetDeviceType>();
      if (tmp == this->dualView.view_host()) returnDevice = false;
    }

    if(returnDevice) {
      DEBUG_UVM_REMOVAL_PRINT_CALLER("getView<Device>ReadOnly");
      if(needsSyncPath()) {
	throwIfHostViewAlive();
	impl::sync_device(originalDualView);
      }
    }
    else {
      DEBUG_UVM_REMOVAL_PRINT_CALLER("getView<Host>ReadOnly");
      if(needsSyncPath()) {
	throwIfDeviceViewAlive();
	impl::sync_host(originalDualView);
      }
    }
    
    return dualView.template view<TargetDeviceType>();
  } 


  template<class TargetDeviceType>
  typename std::remove_reference<decltype(std::declval<DualViewType>().template view<TargetDeviceType>())>::type
  getView (Access::ReadWriteStruct s DEBUG_UVM_REMOVAL_ARGUMENT) const {    
    bool returnDevice = true;
    {
      auto tmp = dualView.template view<TargetDeviceType>();
      if (tmp == this->dualView.view_host()) returnDevice = false;
    }

    if(returnDevice) {
      DEBUG_UVM_REMOVAL_PRINT_CALLER("getView<Device>ReadWrite");
      static_assert(dualViewHasNonConstData,
                    "ReadWrite views are not available for DualView with const data");
      if(needsSyncPath()) {
	throwIfHostViewAlive();
	impl::sync_device(originalDualView);
	originalDualView.modify_device();
      }
    }
    else {
      DEBUG_UVM_REMOVAL_PRINT_CALLER("getView<Host>ReadWrite");
      static_assert(dualViewHasNonConstData,
                    "ReadWrite views are not available for DualView with const data");
      if(needsSyncPath()) {
	throwIfDeviceViewAlive();
	impl::sync_host(originalDualView);
	originalDualView.modify_host();      
      }
    }
    
    return dualView.template view<TargetDeviceType>();
  } 


  template<class TargetDeviceType>
  typename std::remove_reference<decltype(std::declval<DualViewType>().template view<TargetDeviceType>())>::type
  getView (Access::OverwriteAllStruct s DEBUG_UVM_REMOVAL_ARGUMENT) const {    
    if (iAmASubview())
      return getView<TargetDeviceType>(Access::ReadWrite);

    bool returnDevice = true;
    {
      auto tmp = dualView.template view<TargetDeviceType>();
      if (tmp == this->dualView.view_host()) returnDevice = false;
    }
    
    if(returnDevice) {
      DEBUG_UVM_REMOVAL_PRINT_CALLER("getView<Device>OverwriteAll");
      static_assert(dualViewHasNonConstData,
                    "OverwriteAll views are not available for DualView with const data");
      if(needsSyncPath()) {
	throwIfHostViewAlive();
	dualView.clear_sync_state();
	dualView.modify_host();
      }
    }
    else { 
      DEBUG_UVM_REMOVAL_PRINT_CALLER("getView<Host>OverwriteAll");
      static_assert(dualViewHasNonConstData,
                    "OverwriteAll views are not available for DualView with const data");
      if(needsSyncPath()) {
	throwIfDeviceViewAlive();
	dualView.clear_sync_state();
	dualView.modify_device();     
      }
    }
    
    return dualView.template view<TargetDeviceType>();
  } 


  typename HostViewType::const_type
  getHostSubview(int offset, int numEntries, Access::ReadOnlyStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  ) const 
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getHostSubviewReadOnly");
    if(needsSyncPath()) {
      throwIfDeviceViewAlive();
      impl::sync_host(originalDualView);
    }
    return getSubview(dualView.view_host(), offset, numEntries);
  }

  HostViewType
  getHostSubview(int offset, int numEntries, Access::ReadWriteStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  ) 
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getHostSubviewReadWrite");
    static_assert(dualViewHasNonConstData,
        "ReadWrite views are not available for DualView with const data");
    if(needsSyncPath()) {
      throwIfDeviceViewAlive();
      impl::sync_host(originalDualView);
      originalDualView.modify_host();
    }
    return getSubview(dualView.view_host(), offset, numEntries);
  }

  HostViewType
  getHostSubview(int offset, int numEntries, Access::OverwriteAllStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  ) 
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getHostSubviewOverwriteAll");
    static_assert(dualViewHasNonConstData,
        "OverwriteAll views are not available for DualView with const data");
    return getHostSubview(offset, numEntries, Access::ReadWrite);
  }

  typename DeviceViewType::const_type
  getDeviceSubview(int offset, int numEntries, Access::ReadOnlyStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  ) const
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getDeviceSubviewReadOnly");
    if(needsSyncPath()) {
      throwIfHostViewAlive();
      impl::sync_device(originalDualView);
    }
    return getSubview(dualView.view_device(), offset, numEntries);
  }

  DeviceViewType
  getDeviceSubview(int offset, int numEntries, Access::ReadWriteStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  ) 
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getDeviceSubviewReadWrite");
    static_assert(dualViewHasNonConstData,
        "ReadWrite views are not available for DualView with const data");
    if(needsSyncPath()) {
      throwIfHostViewAlive();
      impl::sync_device(originalDualView);
      originalDualView.modify_device();
    }
    return getSubview(dualView.view_device(), offset, numEntries);
  }

  DeviceViewType
  getDeviceSubview(int offset, int numEntries, Access::OverwriteAllStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  ) 
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getDeviceSubviewOverwriteAll");
    static_assert(dualViewHasNonConstData,
        "OverwriteAll views are not available for DualView with const data");
    return getDeviceSubview(offset, numEntries, Access::ReadWrite);
  }


  // A Kokkos implementation of WrappedDualView will have to make these
  // functions publically accessable, but in the Tpetra version, I'm
  // not sure we want this.  There are two options on proceeding here:
  // 1) Mark these as "Expert" only in the comments.  This will be consistent
  //    with a future Kokkos version.
  // 2) Make these protected and then friend Vector/MultiVector.  This will
  //    keep out users from shooting themselves in the feet.  
  const  DualViewType getOriginalDualView() const {
    return originalDualView;
  }

  const DualViewType getDualView() const {
    return dualView;
  }



private:
  template <typename ViewType>
  ViewType getSubview(ViewType view, int offset, int numEntries) const {
    return Kokkos::subview(view, Kokkos::pair<int, int>(offset, offset+numEntries));
  }

  template <typename ViewType,typename int_type>
  ViewType getSubview2D(ViewType view, Kokkos::pair<int_type,int_type> offset0, const Kokkos::Impl::ALL_t&) const {
    return Kokkos::subview(view,offset0,Kokkos::ALL());
  }

  template <typename ViewType,typename int_type>
  ViewType getSubview2D(ViewType view, const Kokkos::Impl::ALL_t&, Kokkos::pair<int_type,int_type> offset1) const {
    return Kokkos::subview(view,Kokkos::ALL(),offset1);
  }

  template <typename ViewType,typename int_type>
  ViewType getSubview2D(ViewType view, Kokkos::pair<int_type,int_type> offset0, Kokkos::pair<int_type,int_type> offset1) const {
    return Kokkos::subview(view,offset0,offset1);
  }


  bool memoryIsAliased() const {
    return deviceMemoryIsHostAccessible && dualView.h_view.data() == dualView.d_view.data();
  }

  bool needsSyncPath() const {
#ifdef KOKKOS_ENABLE_CUDA
    return std::is_same<typename t_dev::memory_space,Kokkos::CudaUVMSpace>::value || !memoryIsAliased();
#else
    return !memoryIsAliased();
#endif
  }

  void throwIfHostViewAlive() const {
    if (dualView.h_view.use_count() > dualView.d_view.use_count()) {
      std::ostringstream msg;
      msg << "Tpetra::Details::WrappedDualView (name = " << dualView.d_view.label() 
          << "; host use_count = " << dualView.h_view.use_count()
          << "; device use_count = " << dualView.d_view.use_count() << "): "
          << "Cannot access data on device while a host view is alive";
      throw std::runtime_error(msg.str());
    }
  }

  void throwIfDeviceViewAlive() const {  
    if (dualView.d_view.use_count() > dualView.h_view.use_count()) {
      std::ostringstream msg;
      msg << "Tpetra::Details::WrappedDualView (name = " << dualView.d_view.label()
          << "; host use_count = " << dualView.h_view.use_count()
          << "; device use_count = " << dualView.d_view.use_count() << "): "
          << "Cannot access data on host while a device view is alive";
      throw std::runtime_error(msg.str());
    }
  }

  bool iAmASubview() {
    return originalDualView.h_view != dualView.h_view;
  }

  mutable DualViewType originalDualView;
  mutable DualViewType dualView;
};

} // namespace Details

} // namespace Tpetra

#endif
