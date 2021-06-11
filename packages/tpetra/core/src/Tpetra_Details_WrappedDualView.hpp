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

// #define DEBUG_UVM_REMOVAL  // Works only with gcc > 4.8

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

private:
  static constexpr bool dualViewHasNonConstData = !impl::hasConstData<DualViewType>::value;
  static constexpr bool deviceMemoryIsHostAccessible =
    Kokkos::SpaceAccessibility<Kokkos::Serial, typename DeviceViewType::memory_space>::accessible;

public:
  WrappedDualView() {}

  WrappedDualView(DualViewType dualV)
    : originalDualView(dualV),
      dualView(originalDualView)
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

  size_t extent(const int i) const {
    return dualView.extent(i);
  }

  typename HostViewType::const_type
  getHostView(Access::ReadOnlyStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  ) const 
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getHostViewReadOnly");
    throwIfDeviceViewAlive();
    impl::sync_host(originalDualView);
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
    throwIfDeviceViewAlive();
    impl::sync_host(originalDualView);
    originalDualView.modify_host();
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
    throwIfDeviceViewAlive();
    if (deviceMemoryIsHostAccessible) Kokkos::fence();
    dualView.clear_sync_state();
    dualView.modify_host();
    return dualView.view_host();
  }

  typename DeviceViewType::const_type
  getDeviceView(Access::ReadOnlyStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  ) const 
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getDeviceViewReadOnly");
    throwIfHostViewAlive();
    impl::sync_device(originalDualView);
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
    throwIfHostViewAlive();
    impl::sync_device(originalDualView);
    originalDualView.modify_device();
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
    throwIfHostViewAlive();
    dualView.clear_sync_state();
    dualView.modify_device();
    return dualView.view_device();
  }

  typename HostViewType::const_type
  getHostSubview(int offset, int numEntries, Access::ReadOnlyStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  ) const 
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getHostSubviewReadOnly");
    throwIfDeviceViewAlive();
    impl::sync_host(originalDualView);
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
    throwIfDeviceViewAlive();
    impl::sync_host(originalDualView);
    originalDualView.modify_host();
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
    throwIfHostViewAlive();
    impl::sync_device(originalDualView);
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
    throwIfHostViewAlive();
    impl::sync_device(originalDualView);
    originalDualView.modify_device();
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

private:
  template <typename ViewType>
  ViewType getSubview(ViewType view, int offset, int numEntries) const {
    return Kokkos::subview(view, Kokkos::pair<int, int>(offset, offset+numEntries));
  }

  void throwIfHostViewAlive() const {
    if( deviceMemoryIsHostAccessible && dualView.h_view.data() == dualView.d_view.data()) return;

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
    if(deviceMemoryIsHostAccessible && dualView.h_view.data() == dualView.d_view.data()) return;

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
