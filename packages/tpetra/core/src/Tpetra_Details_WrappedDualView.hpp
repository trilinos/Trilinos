// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_WRAPPEDDUALVIEW_HPP
#define TPETRA_DETAILS_WRAPPEDDUALVIEW_HPP

#include <Tpetra_Access.hpp>
#include <Tpetra_Details_temporaryViewUtils.hpp>
#include <Kokkos_DualView.hpp>
#include "Teuchos_TestForException.hpp"
#include "Tpetra_Details_ExecutionSpaces.hpp"
#include <sstream>

//#define DEBUG_UVM_REMOVAL  // Works only with gcc > 4.8

#ifdef DEBUG_UVM_REMOVAL

#define DEBUG_UVM_REMOVAL_ARGUMENT ,const char* callerstr = __builtin_FUNCTION(),const char * filestr=__builtin_FILE(),const int linnum = __builtin_LINE()

#define DEBUG_UVM_REMOVAL_PRINT_CALLER(fn) \
  { \
  auto envVarSet = std::getenv("TPETRA_UVM_REMOVAL"); \
  if (envVarSet && (std::strcmp(envVarSet,"1") == 0)) \
    std::cout << (fn) << " called from " << callerstr \
              << " at " << filestr << ":"<<linnum \
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

  // We really need this forward declaration here for friend to work
  template<typename SC, typename LO, typename GO, typename NO>
  class MultiVector;


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
  // This will sync, but only if needed
  dualView.sync_host();
}

template <typename DualViewType>
enableIfConstData<DualViewType>
sync_host(DualViewType dualView) { }

template <typename DualViewType>
enableIfNonConstData<DualViewType>
sync_device(DualViewType dualView) {
  // This will sync, but only if needed
    dualView.sync_device();
}

template <typename DualViewType>
enableIfConstData<DualViewType>
sync_device(DualViewType dualView) { }

}// end namespace Impl

/// \brief Whether WrappedDualView reference count checking is enabled. Initially true.
/// Since the DualView sync functions are not thread-safe, tracking should be disabled
/// during host-parallel regions where WrappedDualView is used.

extern bool wdvTrackingEnabled;

/// \brief Disable WrappedDualView reference-count tracking and syncing.
/// Call this before entering a host-parallel region that uses WrappedDualView.
/// For each WrappedDualView used in the parallel region, its view must be accessed
/// (e.g. getHostView...) before disabling the tracking, so that it may be synced and marked modified correctly.
void disableWDVTracking();

//! Enable WrappedDualView reference-count tracking and syncing. Call this after exiting a host-parallel region that uses WrappedDualView.
void enableWDVTracking();

/// \brief A wrapper around Kokkos::DualView to safely manage data
///        that might be replicated between host and device.
template <typename DualViewType>
class WrappedDualView {
public:

  using DVT = DualViewType;
  using t_host = typename DualViewType::t_host;
  using t_dev  = typename DualViewType::t_dev;

  using HostType   = typename t_host::device_type;
  using DeviceType = typename t_dev::device_type;

private:
  static constexpr bool dualViewHasNonConstData = !impl::hasConstData<DualViewType>::value;
  static constexpr bool deviceMemoryIsHostAccessible =
    Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace, typename t_dev::memory_space>::accessible;

private:
  template <typename>
  friend class WrappedDualView;

public:
  WrappedDualView() {}

  WrappedDualView(DualViewType dualV)
    : originalDualView(dualV),
      dualView(originalDualView)
  { }

  //! Conversion copy constructor.
  template <class SrcDualViewType>
  WrappedDualView(const WrappedDualView<SrcDualViewType>& src)
    : originalDualView(src.originalDualView),
      dualView(src.dualView)
  { }
  
  //! Conversion assignment operator.
  template <class SrcDualViewType>
  WrappedDualView& operator=(const WrappedDualView<SrcDualViewType>& src) {
    originalDualView = src.originalDualView;
    dualView = src.dualView;
    return *this;
  }

  // This is an expert-only constructor
  // For WrappedDualView to manage synchronizations correctly,
  // it must have an DualView which is not a subview to due the
  // sync's on.  This is what origDualV is for.  In this case,
  // dualV is a subview of origDualV.
  WrappedDualView(DualViewType dualV,DualViewType origDualV)
    : originalDualView(origDualV),
      dualView(dualV)
  { }


  WrappedDualView(const t_dev deviceView) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        deviceView.data() != nullptr && deviceView.use_count() == 0,
        std::invalid_argument,
        "Tpetra::Details::WrappedDualView: cannot construct with a device view that\n"
        "does not own its memory (i.e. constructed with a raw pointer and dimensions)\n"
        "because the WrappedDualView needs to assume ownership of the memory.");
    //If the provided view is default-constructed (null, 0 extent, 0 use count),
    //leave the host mirror default-constructed as well in order to have a matching use count of 0.
    t_host hostView;
    if(deviceView.use_count() != 0)
    {
      hostView = Kokkos::create_mirror_view(
          Kokkos::WithoutInitializing,
          typename t_host::memory_space(),
          deviceView);
    }
    originalDualView = DualViewType(deviceView, hostView);
    originalDualView.clear_sync_state();
    originalDualView.modify_device();
    dualView = originalDualView;
  }

  // 1D View constructors
  WrappedDualView(const WrappedDualView parent, int offset, int numEntries) {
    originalDualView = parent.originalDualView;
    dualView = getSubview(parent.dualView, offset, numEntries);
  }


  // 2D View Constructors
  WrappedDualView(const WrappedDualView parent,const Kokkos::pair<size_t,size_t>& rowRng, const Kokkos::ALL_t& colRng) {
    originalDualView = parent.originalDualView;
    dualView = getSubview2D(parent.dualView,rowRng,colRng);
  }

  WrappedDualView(const WrappedDualView parent,const Kokkos::ALL_t &rowRng, const Kokkos::pair<size_t,size_t>& colRng) {
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


  size_t origExtent(const int i) const {
    return originalDualView.h_view.extent(i);
  }

  const char * label() const {
    return dualView.d_view.label();
  }


  typename t_host::const_type
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

  t_host
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

  t_host
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
      throwIfDeviceViewAlive();
      if (deviceMemoryIsHostAccessible) Kokkos::fence("WrappedDualView::getHostView");
      dualView.clear_sync_state();
      dualView.modify_host();
    }
    return dualView.view_host();
  }

  typename t_dev::const_type
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

  t_dev
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

  t_dev
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
      throwIfHostViewAlive();
      if (deviceMemoryIsHostAccessible) Kokkos::fence("WrappedDualView::getDeviceView");
      dualView.clear_sync_state();
      dualView.modify_device();
    }
    return dualView.view_device();
  }

  template<class TargetDeviceType>
  typename std::remove_reference<decltype(std::declval<DualViewType>().template view<TargetDeviceType>())>::type::const_type
  getView (Access::ReadOnlyStruct s DEBUG_UVM_REMOVAL_ARGUMENT) const {
    using ReturnViewType = typename std::remove_reference<decltype(std::declval<DualViewType>().template view<TargetDeviceType>())>::type::const_type;
    using ReturnDeviceType = typename ReturnViewType::device_type;
    constexpr bool returnDevice = std::is_same<ReturnDeviceType, DeviceType>::value;
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
    using ReturnViewType = typename std::remove_reference<decltype(std::declval<DualViewType>().template view<TargetDeviceType>())>::type;
    using ReturnDeviceType = typename ReturnViewType::device_type;
    constexpr bool returnDevice = std::is_same<ReturnDeviceType, DeviceType>::value;

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
    using ReturnViewType = typename std::remove_reference<decltype(std::declval<DualViewType>().template view<TargetDeviceType>())>::type;
    using ReturnDeviceType = typename ReturnViewType::device_type;

    if (iAmASubview())
      return getView<TargetDeviceType>(Access::ReadWrite);

    constexpr bool returnDevice = std::is_same<ReturnDeviceType, DeviceType>::value;

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


  typename t_host::const_type
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

  t_host
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

  t_host
  getHostSubview(int offset, int numEntries, Access::OverwriteAllStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  )
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getHostSubviewOverwriteAll");
    static_assert(dualViewHasNonConstData,
        "OverwriteAll views are not available for DualView with const data");
    return getHostSubview(offset, numEntries, Access::ReadWrite);
  }

  typename t_dev::const_type
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

  t_dev
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

  t_dev
  getDeviceSubview(int offset, int numEntries, Access::OverwriteAllStruct
    DEBUG_UVM_REMOVAL_ARGUMENT
  )
  {
    DEBUG_UVM_REMOVAL_PRINT_CALLER("getDeviceSubviewOverwriteAll");
    static_assert(dualViewHasNonConstData,
        "OverwriteAll views are not available for DualView with const data");
    return getDeviceSubview(offset, numEntries, Access::ReadWrite);
  }


  // Debugging functions to get copies of the view state
  typename t_host::HostMirror getHostCopy() const {
    auto X_dev = dualView.view_host();
    if(X_dev.span_is_contiguous()) {
      auto mirror = Kokkos::create_mirror_view(X_dev);
      Kokkos::deep_copy(mirror,X_dev);
      return mirror;
    }
    else {
      auto X_contig = Tpetra::Details::TempView::toLayout<decltype(X_dev), Kokkos::LayoutLeft>(X_dev);
      auto mirror = Kokkos::create_mirror_view(X_contig);
      Kokkos::deep_copy(mirror,X_contig);
      return mirror;
    }
  }

  typename t_dev::HostMirror getDeviceCopy() const {
    auto X_dev = dualView.view_device();
    if(X_dev.span_is_contiguous()) {
      auto mirror = Kokkos::create_mirror_view(X_dev);
      Kokkos::deep_copy(mirror,X_dev);
      return mirror;
    }
    else {
      auto X_contig = Tpetra::Details::TempView::toLayout<decltype(X_dev), Kokkos::LayoutLeft>(X_dev);
      auto mirror = Kokkos::create_mirror_view(X_contig);
      Kokkos::deep_copy(mirror,X_contig);
      return mirror;
    }
  }

  // Debugging functions for validity checks
  bool is_valid_host() const {
    return dualView.view_host().size() == 0   || dualView.view_host().data();
  }

  bool is_valid_device() const {
    return dualView.view_device().size() == 0 || dualView.view_device().data();
  }


  bool need_sync_host() const {
    return originalDualView.need_sync_host();
  }

  bool need_sync_device() const {
    return originalDualView.need_sync_device();
  }

  int host_view_use_count() const {
    return originalDualView.h_view.use_count();
  }

  int device_view_use_count() const {
    return originalDualView.d_view.use_count();
  }


  // MultiVector really needs to get at the raw DualViews,
  // but we'd very much prefer that users not.
  template<typename SC, typename LO, typename GO, typename NO>
  friend class ::Tpetra::MultiVector;

private:
  // A Kokkos implementation of WrappedDualView will have to make these
  // functions publically accessable, but in the Tpetra version, we'd
  // really rather not.
  DualViewType getOriginalDualView() const {
    return originalDualView;
  }

  DualViewType getDualView() const {
    return dualView;
  }

  template <typename ViewType>
  ViewType getSubview(ViewType view, int offset, int numEntries) const {
    return Kokkos::subview(view, Kokkos::pair<int, int>(offset, offset+numEntries));
  }

  template <typename ViewType,typename int_type>
  ViewType getSubview2D(ViewType view, Kokkos::pair<int_type,int_type> offset0, const Kokkos::ALL_t&) const {
    return Kokkos::subview(view,offset0,Kokkos::ALL());
  }

  template <typename ViewType,typename int_type>
  ViewType getSubview2D(ViewType view, const Kokkos::ALL_t&, Kokkos::pair<int_type,int_type> offset1) const {
    return Kokkos::subview(view,Kokkos::ALL(),offset1);
  }

  template <typename ViewType,typename int_type>
  ViewType getSubview2D(ViewType view, Kokkos::pair<int_type,int_type> offset0, Kokkos::pair<int_type,int_type> offset1) const {
    return Kokkos::subview(view,offset0,offset1);
  }

  bool memoryIsAliased() const {
    return deviceMemoryIsHostAccessible && dualView.h_view.data() == dualView.d_view.data();
  }


  /// \brief needsSyncPath tells us whether we need the "sync path" where we (potentially) fence,
  ///        check use counts and take care of sync/modify for the underlying DualView.
  ///
  /// The logic is this:
  /// 1. If WrappedDualView tracking is disabled, then never take the sync path.
  /// 2. For non-GPU architectures where the host/device pointers are aliased
  ///    we don't need the "sync path."
  /// 3. For GPUs, we always need the "sync path."  For shared host/device memory (e.g. CudaUVM)
  ///    the Kokkos::deep_copy in the sync is a no-op, but the fence associated with it matters.
  ///
  ///
  /// Avoiding the "sync path" speeds up calculations on architectures where we can
  /// avoid it (e.g. SerialNode) by not not touching the modify flags.
  ///
  /// Note for the future: Memory spaces that can be addressed on both host and device
  /// that don't otherwise have an intrinsic fencing mechanism will need to trigger the
  /// "sync path"
  bool needsSyncPath() const {
    if(!wdvTrackingEnabled)
      return false;

    // We check to see if the memory is not aliased *or* if it is a supported
    // (heterogeneous memory) accelerator (for shared host/device memory).
    return !memoryIsAliased() || Spaces::is_gpu_exec_space<typename DualViewType::execution_space>();
  }


  void throwIfViewsAreDifferentSizes() const {    
    // Here we check *size* (the product of extents) rather than each extent individually.
    // This is mostly designed to catch people resizing one view, but not the other.
    if(dualView.d_view.size() != dualView.h_view.size()) {    
        std::ostringstream msg;
        msg << "Tpetra::Details::WrappedDualView (name = " << dualView.d_view.label()
            << "; host and device views are different sizes: "
            << dualView.h_view.size() << " vs " <<dualView.h_view.size();
        throw std::runtime_error(msg.str());
    }
  }

  void throwIfHostViewAlive() const {
    throwIfViewsAreDifferentSizes();
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
    throwIfViewsAreDifferentSizes();
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
