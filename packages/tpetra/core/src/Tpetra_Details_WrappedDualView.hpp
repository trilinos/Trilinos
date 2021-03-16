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
#include <sstream>

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
private:
  static constexpr bool dualViewHasNonConstData = !impl::hasConstData<DualViewType>::value;
  using HostViewType = typename DualViewType::t_host;
  using DeviceViewType = typename DualViewType::t_dev;

public:
  WrappedDualView() {}

  WrappedDualView(DualViewType dualV)
    : dualView(dualV)
  { }

  WrappedDualView(const DeviceViewType devView)
  {
     HostViewType hostView =
       Kokkos::create_mirror_view_and_copy(typename HostViewType::memory_space(), devView);
     dualView = DualViewType(devView, hostView);
  }

  size_t extent(const int i) const {
    return dualView.extent(i);
  }

  typename DualViewType::t_host::const_type
  getHostView(Access::ReadOnlyStruct) const {
    throwIfDeviceViewAlive();
    impl::sync_host(dualView);
    return dualView.view_host();
  }

  typename DualViewType::t_host
  getHostView(Access::ReadWriteStruct) {
    static_assert(dualViewHasNonConstData,
        "ReadWrite views are not available for DualView with const data");
    throwIfDeviceViewAlive();
    dualView.sync_host();
    dualView.modify_host();
    return dualView.view_host();
  }

  typename DualViewType::t_host
  getHostView(Access::WriteOnlyStruct) {
    static_assert(dualViewHasNonConstData,
        "WriteOnly views are not available for DualView with const data");
    throwIfDeviceViewAlive();
    dualView.clear_sync_state();
    dualView.modify_host();
    return dualView.view_host();
  }

  typename DualViewType::t_dev::const_type
  getDeviceView(Access::ReadOnlyStruct) const {
    throwIfHostViewAlive();
    impl::sync_device(dualView);
    return dualView.view_device();
  }

  typename DualViewType::t_dev
  getDeviceView(Access::ReadWriteStruct) {
    static_assert(dualViewHasNonConstData,
        "ReadWrite views are not available for DualView with const data");
    throwIfHostViewAlive();
    dualView.sync_device();
    dualView.modify_device();
    return dualView.view_device();
  }

  typename DualViewType::t_dev
  getDeviceView(Access::WriteOnlyStruct) {
    static_assert(dualViewHasNonConstData,
        "WriteOnly views are not available for DualView with const data");
    throwIfHostViewAlive();
    dualView.clear_sync_state();
    dualView.modify_device();
    return dualView.view_device();
  }

  typename DualViewType::t_host::const_type
  getHostSubview(int offset, int numEntries, Access::ReadOnlyStruct) const {
    throwIfDeviceViewAlive();
    impl::sync_host(dualView);
    return getSubview(dualView.view_host(), offset, numEntries);
  }

  typename DualViewType::t_host
  getHostSubview(int offset, int numEntries, Access::ReadWriteStruct) {
    static_assert(dualViewHasNonConstData,
        "ReadWrite views are not available for DualView with const data");
    throwIfDeviceViewAlive();
    dualView.sync_host();
    dualView.modify_host();
    return getSubview(dualView.view_host(), offset, numEntries);
  }

  typename DualViewType::t_host
  getHostSubview(int offset, int numEntries, Access::WriteOnlyStruct) {
    static_assert(dualViewHasNonConstData,
        "WriteOnly views are not available for DualView with const data");
    return getHostSubview(offset, numEntries, Access::ReadWrite);
  }

  typename DualViewType::t_dev::const_type
  getDeviceSubview(int offset, int numEntries, Access::ReadOnlyStruct) const {
    throwIfHostViewAlive();
    impl::sync_device(dualView);
    return getSubview(dualView.view_device(), offset, numEntries);
  }

  typename DualViewType::t_dev
  getDeviceSubview(int offset, int numEntries, Access::ReadWriteStruct) {
    static_assert(dualViewHasNonConstData,
        "ReadWrite views are not available for DualView with const data");
    throwIfHostViewAlive();
    dualView.sync_device();
    dualView.modify_device();
    return getSubview(dualView.view_device(), offset, numEntries);
  }

  typename DualViewType::t_dev
  getDeviceSubview(int offset, int numEntries, Access::WriteOnlyStruct) {
    static_assert(dualViewHasNonConstData,
        "WriteOnly views are not available for DualView with const data");
    return getDeviceSubview(offset, numEntries, Access::ReadWrite);
  }

private:
  template <typename ViewType>
  ViewType getSubview(ViewType view, int offset, int numEntries) const {
    return Kokkos::subview(view, Kokkos::pair<int, int>(offset, offset+numEntries));
  }

  void throwIfHostViewAlive() const {
    if (dualView.h_view.use_count() > dualView.d_view.use_count()) {
      std::ostringstream msg;
      msg << "Tpetra::Details::WrappedDualView (name = " << dualView.d_view.label() << "): "
          << "Cannot access data on device while a host view is alive";
      throw std::runtime_error(msg.str());
    }
  }

  void throwIfDeviceViewAlive() const {
    if (dualView.d_view.use_count() > dualView.h_view.use_count()) {
      std::ostringstream msg;
      msg << "Tpetra::Details::WrappedDualView (name = " << dualView.d_view.label() << "): "
          << "Cannot access data on host while a device view is alive";
      throw std::runtime_error(msg.str());
    }
  }

  mutable DualViewType dualView;
};

} // namespace Details

} // namespace Tpetra

#endif
