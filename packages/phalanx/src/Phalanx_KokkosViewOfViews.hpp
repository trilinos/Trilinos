// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_KOKKOS_VIEW_OF_VIEWS_HPP
#define PHALANX_KOKKOS_VIEW_OF_VIEWS_HPP

#include "Sacado.hpp" // for IsADType
#include <utility> // for declval

namespace PHX {

  // ****************************
  // Add pointer (used to construct the static data type)
  // ****************************

  namespace v_of_v_utils {
    template<typename Data,int N> struct add_pointer;

    template<typename Data,int N> struct add_pointer
    { using type = typename add_pointer<Data*,N-1>::type; };

    template<typename Data> struct add_pointer<Data,0>
    { using type = Data; };
  }

  // ****************************
  // Functions to free inner view memory (assumes outer views are on host)
  // ****************************
  template<typename InnerViewDataType,typename... InnerProps,typename... OuterProps>
  auto freeInnerViewsOfHostHostViewOfViews(Kokkos::View<Kokkos::View<InnerViewDataType,InnerProps...>*,OuterProps...>& v_of_v) {
    for (std::size_t i=0; i < v_of_v.extent(0); ++i) {
      v_of_v(i) = Kokkos::View<InnerViewDataType,InnerProps...>();
    }
  }

  template<typename InnerViewDataType,typename... InnerProps,typename... OuterProps>
  auto freeInnerViewsOfHostHostViewOfViews(Kokkos::View<Kokkos::View<InnerViewDataType,InnerProps...>**,OuterProps...>& v_of_v) {
    for (std::size_t i=0; i < v_of_v.extent(0); ++i) {
      for (std::size_t j=0; j < v_of_v.extent(1); ++j) {
        v_of_v(i,j) = Kokkos::View<InnerViewDataType,InnerProps...>();
      }
    }
  }

  template<typename InnerViewDataType,typename... InnerProps,typename... OuterProps>
  auto freeInnerViewsOfHostHostViewOfViews(Kokkos::View<Kokkos::View<InnerViewDataType,InnerProps...>***,OuterProps...>& v_of_v) {
    for (std::size_t i=0; i < v_of_v.extent(0); ++i) {
      for (std::size_t j=0; j < v_of_v.extent(1); ++j) {
        for (std::size_t k=0; k < v_of_v.extent(2); ++k) {
          v_of_v(i,j,k) = Kokkos::View<InnerViewDataType,InnerProps...>();
        }
      }
    }
  }

  // ****************************
  // ViewOfViews: third version (inner views are runtime unmanaged - no Unmanaged template parameter)
  // ****************************

  /** Wrapper class that correctly handles ViewOfViews construction
      and object lifetime. This class makes sure the host view stays
      in scope for the life of the device view and makes sure that the
      device is synced to host before use.

      Main restrictions:

      1. When UVM is not used in the outer view, we need to allocate
      the outer VofV on host and copy to device to initialize the
      inner views correctly (tracking object).

      2. Step 1 means that the host view must exist as long as the
      device view is being used, otherwise the views may go out of
      scope. This object exists to pair up the host and device view to
      make sure the inner views are not deleted early.

      3. Normally we use an unmanaged view (constructed with the
      Unmanaged template parameter) for the inner views to prevent
      double deletion. However, there are use cases where it's painful
      to pass around views built with the unmanaged template parameter
      (libraries with functions that block the unmanaged argument). We
      can generate an unmanged view without the template parameter by
      constructing the view with a raw pointer. This thrid
      implementation does that here.
  */
  template<int OuterViewRank,typename InnerViewType,typename... OuterViewProps>
  class ViewOfViews {

  public:
    // Layout of the outer view doesn't matter for performance so we
    // use a raw Kokkos::View instead of PHX::View. The inner views are
    // what is important for performance.
    using OuterDataType = typename PHX::v_of_v_utils::add_pointer<InnerViewType,OuterViewRank>::type;
    using OuterViewType = Kokkos::View<OuterDataType,OuterViewProps...>;

  private:
    // Inner views are mananged - used to prevent early deletion
    typename OuterViewType::HostMirror view_host_;
    // Inner views are unmanaged by runtime construction with pointer
    // (avoids template parameter). Used to correctly initialize outer
    // device view on device.
    typename OuterViewType::HostMirror view_host_unmanaged_;
    // Device view
    OuterViewType view_device_;
    // True if the host view has not been synced to device
    bool device_view_is_synced_;
    // True if the outer view has been initialized
    bool is_initialized_;
    // Use count of device view after initialization. This changes based on whether the view_device_ is accessible to host space. If not accessible, the use_count_ is 1. If it is accessible, the value is 2 due to using create_mirror_view in initialization if view_host_unmanaged_.
    int use_count_;
    // A safety check. If true, this makes sure there are no external references to the device view of views.
    bool check_use_count_;

  public:
    /// Ctor that uses the default execution space instance.
    template<typename... Extents>
    ViewOfViews(const std::string name,Extents... extents)
      : view_host_(name,extents...),
        view_device_(name,extents...),
        device_view_is_synced_(false),
        is_initialized_(true),
        use_count_(0),
        check_use_count_(true)
    {
      view_host_unmanaged_ = Kokkos::create_mirror_view(view_device_);
      use_count_ = view_device_.impl_track().use_count();
    }

    /// Ctor that uses a user specified execution space instance.
    /// NOTE: Consistent with Kokkos, when a user supplies the
    /// execution space instance, the function does not internally
    /// fence. Be sure to manually fence as needed.
    template<typename ExecSpace,typename... Extents>
    ViewOfViews(const ExecSpace& e_space,const std::string name,Extents... extents)
      : view_host_(Kokkos::view_alloc(typename OuterViewType::HostMirror::execution_space(),name),extents...),
        view_device_(Kokkos::view_alloc(e_space,name),extents...),
        device_view_is_synced_(false),
        is_initialized_(true),
        use_count_(0),
        check_use_count_(true)
    {
      view_host_unmanaged_ = Kokkos::create_mirror_view(Kokkos::view_alloc(typename OuterViewType::HostMirror::execution_space()),view_device_);
      use_count_ = view_device_.impl_track().use_count();
    }

    ViewOfViews()
      : device_view_is_synced_(false),
        is_initialized_(false),
        use_count_(0),
        check_use_count_(true)
    {}

    ViewOfViews(const ViewOfViews<OuterViewRank,InnerViewType,OuterViewProps...>& ) = default;
    ViewOfViews& operator=(const ViewOfViews<OuterViewRank,InnerViewType,OuterViewProps...>& ) = default;
    ViewOfViews(ViewOfViews<OuterViewRank,InnerViewType,OuterViewProps...>&& src) = default;
    ViewOfViews& operator=(ViewOfViews<OuterViewRank,InnerViewType,OuterViewProps...>&& ) = default;

    // Making this a kokkos function eliminates cuda compiler warnings
    // in objects that contain ViewOfViews that are copied to device.
    KOKKOS_INLINE_FUNCTION
    ~ViewOfViews()
    {
      // Make sure there is not another object pointing to the device
      // view if the host view is about to be deleted. The host view
      // may delete the inner views if it is the last owner.
      KOKKOS_IF_ON_HOST((
        if ( check_use_count_ ) {
          if ( (view_host_.impl_track().use_count() == 1) && (view_device_.impl_track().use_count() > use_count_) ) {
            Kokkos::abort("\n ERROR - PHX::ViewOfViews - please free all instances of device ViewOfView \n before deleting the host ViewOfView!\n\n");
          }
        }
        // We must manually delete the inner views now. The Serial
        // backend has been updated for thread safety, causing
        // deadlock from a view dtor mutex - nested views cause nested
        // parallel_for calls. Details:
        // Commit to kokkos that caused issue: https://github.com/kokkos/kokkos/pull/7080
        // Kokkos issue that reports failures: https://github.com/kokkos/kokkos/issues/7092
        // Phalanx: fixes: https://github.com/trilinos/Trilinos/pull/13188
        if (view_host_.impl_track().use_count() == 1) {
          PHX::freeInnerViewsOfHostHostViewOfViews(view_host_);
        }
      ))
    }

    /// Enable safety check in dtor for external references.
    void enableSafetyCheck() { check_use_count_ = true; }

    /// Disable safety check in dtor for external references.
    void disableSafetyCheck() { check_use_count_ = false; }

    /// Allocate the out view objects. Extents are for the outer view. Uses the default execution space.
    template<typename... Extents>
    void initialize(const std::string name,Extents... extents)
    {
      view_host_ = typename OuterViewType::HostMirror(name,extents...);
      view_device_ = OuterViewType(name,extents...);
      view_host_unmanaged_ = Kokkos::create_mirror_view(view_device_);
      device_view_is_synced_ = false;
      is_initialized_ = true;
      use_count_ = view_device_.impl_track().use_count();
    }

    /// Allocate the out view objects. Extents are for the outer
    /// view. Uses a user supplied execution space.  NOTE: Consistent
    /// with Kokkos, when a user supplies the execution space
    /// instance, the function does not internally fence. Be sure to
    /// manually fence as needed.
    template<typename ExecSpace,typename... Extents>
    void initialize(const ExecSpace& e_space,const std::string name,Extents... extents)
    {
      view_host_ = typename OuterViewType::HostMirror(Kokkos::view_alloc(typename OuterViewType::HostMirror::execution_space(),name),extents...);
      view_device_ = OuterViewType(Kokkos::view_alloc(e_space,name),extents...);
      view_host_unmanaged_ = Kokkos::create_mirror_view(Kokkos::view_alloc(typename OuterViewType::HostMirror::execution_space()),view_device_);
      device_view_is_synced_ = false;
      is_initialized_ = true;
      use_count_ = view_device_.impl_track().use_count();
    }

    // Returns true if the outer view has been initialized.
    bool isInitialized() const {return is_initialized_;}

    bool deviceViewIsSynced() const {return device_view_is_synced_;}

    bool safetyCheck() const {return check_use_count_;}

    template<typename... Indices>
    void addView(InnerViewType v,Indices... i)
    {
      static_assert(sizeof...(Indices)==OuterViewRank,
        "Error: PHX::ViewOfViews::addView() - the number of indices must match the outer view rank!");

      TEUCHOS_ASSERT(is_initialized_);

      // Store the managed version so inner views don't get deleted.
      view_host_(i...) = v;

      // Store a runtime unmanaged view for deep_copy to
      // device. Unmanaged is required to prevent double deletion on
      // device. For FAD types, we need to adjust the layout to
      // account for the derivative array. The layout() method reutrns
      // the non-fad adjusted layout.
      if (Sacado::IsADType<typename InnerViewType::value_type>::value) {
        auto layout = v.layout();
        layout.dimension[InnerViewType::rank] = Kokkos::dimension_scalar(v);
        view_host_unmanaged_(i...) = InnerViewType(v.data(),layout);
      }
      else
        view_host_unmanaged_(i...) = InnerViewType(v.data(),v.layout());

      device_view_is_synced_ = false;
    }

    /// deep_copy the outer view to device. Uses default execution
    /// space. Note this only syncs the outer view. The inner views
    /// are assumed to be on device for both host and device outer
    /// views.
    void syncHostToDevice()
    {
      TEUCHOS_ASSERT(is_initialized_);
      Kokkos::deep_copy(view_device_,view_host_unmanaged_);
      device_view_is_synced_ = true;
    }

    /// deep_copy the outer view to device. Uses a user supplied
    /// execution space.  Note this only syncs the outer view. The
    /// inner views are assumed to be on device for both host and
    /// device outer views.  NOTE: Consistent with Kokkos, when a user
    /// supplies the execution space instance, the function does not
    /// internally fence. Be sure to manually fence as needed.
    template<typename ExecSpace>
    void syncHostToDevice(const ExecSpace& e_space)
    {
      TEUCHOS_ASSERT(is_initialized_);
      Kokkos::deep_copy(e_space,view_device_,view_host_unmanaged_);
      device_view_is_synced_ = true;
    }

    /// Returns a host mirror view for the outer view, where the inner
    /// views are still on device.
    auto getViewHost()
    {
      TEUCHOS_ASSERT(is_initialized_);
      return view_host_;
    }

    /// Returns a host mirror view for the outer view, where the inner
    /// views are still on device.
    auto getViewHost() const
    {
      TEUCHOS_ASSERT(is_initialized_);
      return view_host_;
    }

    /// Returns device view of views
    auto getViewDevice()
    {
      KOKKOS_ASSERT(device_view_is_synced_);
      return view_device_;
    }

    /// Returns device view of views
    auto getViewDevice() const
    {
      KOKKOS_ASSERT(device_view_is_synced_);
      return view_device_;
    }
  };

  /// Returns a rank-1 view of views where both the outer and inner views are on host. Values are deep_copied from input v_of_v.
  template<typename InnerViewType,typename... OuterViewProps>
  auto createHostHostViewOfViews(const PHX::ViewOfViews<1,InnerViewType,OuterViewProps...>& vov)
  {
    auto outer_host_inner_device_vov = vov.getViewHost();
    using HostInnerViewType = decltype(Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),outer_host_inner_device_vov(0)));
    PHX::ViewOfViews<1,HostInnerViewType,OuterViewProps...> host_host_vov("OuterHostInnerHost VOV: "+vov.getViewHost().label(),
                                                                          vov.getViewHost().extent(0));
    for (size_t i=0; i < outer_host_inner_device_vov.extent(0); ++i) {
      host_host_vov.addView(Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),outer_host_inner_device_vov(i)),i);
    }
    return host_host_vov;
  }

  /// Returns a rank-2 view of views where both the outer and inner views are on host. Values are deep_copied from input v_of_v.
  template<typename InnerViewType,typename... OuterViewProps>
  auto createHostHostViewOfViews(const PHX::ViewOfViews<2,InnerViewType,OuterViewProps...>& vov)
  {
    auto outer_host_inner_device_vov = vov.getViewHost();
    using HostInnerViewType = decltype(Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),outer_host_inner_device_vov(0,0)));
    PHX::ViewOfViews<2,HostInnerViewType,OuterViewProps...> host_host_vov("OuterHostInnerHost VOV: "+vov.getViewHost().label(),
                                                                          vov.getViewHost().extent(0),
                                                                          vov.getViewHost().extent(1));
    for (size_t i=0; i < outer_host_inner_device_vov.extent(0); ++i) {
      for (size_t j=0; j < outer_host_inner_device_vov.extent(1); ++j) {
        host_host_vov.addView(Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),outer_host_inner_device_vov(i,j)),i,j);
      }
    }
    return host_host_vov;
  }

  /// Returns a rank-3 view of views where both the outer and inner views are on host. Values are deep_copied from input v_of_v.
  template<typename InnerViewType,typename... OuterViewProps>
  auto createHostHostViewOfViews(const PHX::ViewOfViews<3,InnerViewType,OuterViewProps...>& vov)
  {
    auto outer_host_inner_device_vov = vov.getViewHost();
    using HostInnerViewType = decltype(Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),outer_host_inner_device_vov(0,0,0)));
    PHX::ViewOfViews<3,HostInnerViewType,OuterViewProps...> host_host_vov("OuterHostInnerHost VOV: "+vov.getViewHost().label(),
                                                                          vov.getViewHost().extent(0),
                                                                          vov.getViewHost().extent(1),
                                                                          vov.getViewHost().extent(2));
    for (size_t i=0; i < outer_host_inner_device_vov.extent(0); ++i) {
      for (size_t j=0; j < outer_host_inner_device_vov.extent(1); ++j) {
        for (size_t k=0; k < outer_host_inner_device_vov.extent(2); ++k) {
          host_host_vov.addView(Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),outer_host_inner_device_vov(i,j,k)),i,j,k);
        }
      }
    }
    return host_host_vov;
  }

  // Old declaration for backwards compatibility
  template<int OuterViewRank,typename InnerViewType,typename... OuterViewProps>
  using ViewOfViews3 = ViewOfViews<OuterViewRank,InnerViewType,OuterViewProps...>;

} // namespace PHX

#endif
