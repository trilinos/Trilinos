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
#include <memory> // for shared_ptr

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

  namespace details {
    struct ViewOfViewsDeleter {
      bool do_safety_check_ = false;

      template <class D, class... P>
      std::enable_if_t<
          Kokkos::is_view_v<typename Kokkos::View<D, P...>::value_type>>
      operator()(Kokkos::View<D, P...>* vov) const {
        Kokkos::fence("PHX:ViewOfViewsDeleter: fence before host cleanup of View of Views");
        if (do_safety_check_) {
          if (vov->use_count() > 1) {
            Kokkos::abort("\n\n********\n ERROR - PHX::ViewOfViews - please free all instances of device Kokkos::View<View<...>> \n before deleting the host ViewOfView!\n********\n\n");
          }
        }

        // destroy inner views on host, outside of a parallel region
        constexpr size_t rank = Kokkos::View<D, P...>::rank();
        static constexpr bool device_view_is_accessible_from_host = Kokkos::SpaceAccessibility<Kokkos::HostSpace, typename Kokkos::View<D, P...>::memory_space>::accessible;
        if (device_view_is_accessible_from_host) {
          if constexpr (rank == 0) {
            (*vov)() = {};
          } else if constexpr (rank == 1) {
            for (size_t i = 0; i < vov->extent(0); ++i) {
              (*vov)(i) = {};
            }
          } else if constexpr (rank == 2) {
            for (size_t i = 0; i < vov->extent(0); ++i) {
              for (size_t j = 0; j < vov->extent(1); ++j) {
                (*vov)(i, j) = {};
              }
            }
          } else if constexpr (rank == 3) {
            for (size_t i = 0; i < vov->extent(0); ++i) {
              for (size_t j = 0; j < vov->extent(1); ++j) {
                for (size_t k = 0; k < vov->extent(2); ++k) {
                  (*vov)(i, j, k) = {};
                }
              }
            }
          } else {
            static_assert(std::is_void_v<decltype(vov)>, "\n\n********\n Error - PHX::ViewOfViews - ViewOfView with outer view greater than rank 3 is not supported!\n********\n\n");
          }
        }
        // dispose of the outer view
        delete vov;
      }
    };

    template <class VoV>
    struct ViewOfViewsMaker {
      static_assert(Kokkos::is_view_v<VoV>);
    };
    template <class D, class... P>
    struct ViewOfViewsMaker<Kokkos::View<D, P...>> {
      static_assert(Kokkos::is_view_v<typename Kokkos::View<D, P...>::value_type>);
      template <class... Args>
      static auto make_shared(Args&&... args) {
        return std::shared_ptr<Kokkos::View<D, P...>>(
            new Kokkos::View<D, P...>((Args &&) args...), ViewOfViewsDeleter());
      }
    };
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
    // Inner views are managed - used to prevent early deletion
    std::shared_ptr<typename OuterViewType::HostMirror> view_host_;
    // Device view
    std::shared_ptr<OuterViewType> view_device_;
    // Inner views are unmanaged by runtime construction with pointer
    // (avoids template parameter). Used to correctly initialize outer
    // device view on device.
    static constexpr bool device_view_is_accessible_from_host = Kokkos::SpaceAccessibility<Kokkos::HostSpace, typename OuterViewType::memory_space>::accessible;
    std::conditional_t<device_view_is_accessible_from_host, std::shared_ptr<OuterViewType>, std::shared_ptr<typename OuterViewType::HostMirror>> view_host_unmanaged_;
    // True if the host view has not been synced to device
    bool device_view_is_synced_;
    // True if the outer view has been initialized
    bool is_initialized_;
    // A safety check. If true, this makes sure there are no external references to the device view of views.
    bool check_use_count_;

  public:
    /// Ctor that uses the default execution space instance.
    template<typename... Extents>
    ViewOfViews(const std::string name,Extents... extents)
      : view_host_(details::ViewOfViewsMaker<typename OuterViewType::HostMirror>::make_shared(name,extents...)),
        view_device_(details::ViewOfViewsMaker<OuterViewType>::make_shared(name,extents...)),
        device_view_is_synced_(false),
        is_initialized_(true),
        check_use_count_(true)
    {
      if constexpr (device_view_is_accessible_from_host) {
        view_host_unmanaged_ = view_device_;
      } else {
        view_host_unmanaged_ = details::ViewOfViewsMaker<typename OuterViewType::HostMirror>::make_shared(Kokkos::create_mirror_view(*view_device_));
      }

      std::get_deleter<details::ViewOfViewsDeleter>(view_device_)->do_safety_check_ = check_use_count_;
    }

    /// Ctor that uses a user specified execution space instance.
    /// NOTE: Consistent with Kokkos, when a user supplies the
    /// execution space instance, the function does not internally
    /// fence. Be sure to manually fence as needed.
    template<typename ExecSpace,typename... Extents>
    ViewOfViews(const ExecSpace& e_space,const std::string name,Extents... extents)
      : view_host_(details::ViewOfViewsMaker<typename OuterViewType::HostMirror>::make_shared(Kokkos::view_alloc(typename OuterViewType::HostMirror::execution_space(),name),extents...)),
        view_device_(details::ViewOfViewsMaker<OuterViewType>::make_shared(Kokkos::view_alloc(e_space,name),extents...)),
        device_view_is_synced_(false),
        is_initialized_(true),
        check_use_count_(true)
    {
      if constexpr (device_view_is_accessible_from_host) {
        view_host_unmanaged_ = view_device_;
      } else {
        view_host_unmanaged_ = details::ViewOfViewsMaker<typename OuterViewType::HostMirror>::make_shared(Kokkos::create_mirror_view(Kokkos::view_alloc(typename OuterViewType::HostMirror::execution_space()),*view_device_));
      }

      std::get_deleter<details::ViewOfViewsDeleter>(view_device_)->do_safety_check_ = check_use_count_;
    }

    ViewOfViews()
      : device_view_is_synced_(false),
        is_initialized_(false),
        check_use_count_(true)
    {}

    /// Enable safety check in dtor for external references.
    void enableSafetyCheck()
    {
      check_use_count_ = true;
      if (this->isInitialized())
        std::get_deleter<details::ViewOfViewsDeleter>(view_device_)->do_safety_check_ = check_use_count_;
    }

    /// Disable safety check in dtor for external references.
    void disableSafetyCheck()
    {
      check_use_count_ = false;
      if (this->isInitialized())
        std::get_deleter<details::ViewOfViewsDeleter>(view_device_)->do_safety_check_ = check_use_count_;
    }

    /// Allocate the out view objects. Extents are for the outer view. Uses the default execution space.
    template<typename... Extents>
    void initialize(const std::string name,Extents... extents)
    {
      view_host_ = details::ViewOfViewsMaker<typename OuterViewType::HostMirror>::make_shared(name,extents...);
      view_device_ = details::ViewOfViewsMaker<OuterViewType>::make_shared(name,extents...);
      if constexpr (device_view_is_accessible_from_host) {
        view_host_unmanaged_ = view_device_;
      } else {
        view_host_unmanaged_ = details::ViewOfViewsMaker<typename OuterViewType::HostMirror>::make_shared(Kokkos::create_mirror_view(*view_device_));
      }

      std::get_deleter<details::ViewOfViewsDeleter>(view_device_)->do_safety_check_ = check_use_count_;

      device_view_is_synced_ = false;
      is_initialized_ = true;
    }

    /// Allocate the out view objects. Extents are for the outer
    /// view. Uses a user supplied execution space.  NOTE: Consistent
    /// with Kokkos, when a user supplies the execution space
    /// instance, the function does not internally fence. Be sure to
    /// manually fence as needed.
    template<typename ExecSpace,typename... Extents>
    void initialize(const ExecSpace& e_space,const std::string name,Extents... extents)
    {
      view_host_ = details::ViewOfViewsMaker<typename OuterViewType::HostMirror>::make_shared(Kokkos::view_alloc(typename OuterViewType::HostMirror::execution_space(),name),extents...);
      view_device_ = details::ViewOfViewsMaker<OuterViewType>::make_shared(Kokkos::view_alloc(e_space,name),extents...);
      if constexpr (device_view_is_accessible_from_host) {
        view_host_unmanaged_ = view_device_;
      } else {
        view_host_unmanaged_ = details::ViewOfViewsMaker<typename OuterViewType::HostMirror>::make_shared(Kokkos::create_mirror_view(Kokkos::view_alloc(typename OuterViewType::HostMirror::execution_space()),*view_device_));
      }

      std::get_deleter<details::ViewOfViewsDeleter>(view_device_)->do_safety_check_ = check_use_count_;

      device_view_is_synced_ = false;
      is_initialized_ = true;
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
      (*view_host_)(i...) = v;

      // Store a runtime unmanaged view for deep_copy to
      // device. Unmanaged is required to prevent double deletion on
      // device. For FAD types, we need to adjust the layout to
      // account for the derivative array. The layout() method reutrns
      // the non-fad adjusted layout.
      if (Sacado::IsADType<typename InnerViewType::value_type>::value) {
        auto layout = v.layout();
        layout.dimension[InnerViewType::rank] = Kokkos::dimension_scalar(v);
        (*view_host_unmanaged_)(i...) = InnerViewType(v.data(),layout);
      }
      else
        (*view_host_unmanaged_)(i...) = InnerViewType(v.data(),v.layout());

      device_view_is_synced_ = false;
    }

    /// deep_copy the outer view to device. Uses default execution
    /// space. Note this only syncs the outer view. The inner views
    /// are assumed to be on device for both host and device outer
    /// views.
    void syncHostToDevice()
    {
      TEUCHOS_ASSERT(is_initialized_);
      Kokkos::deep_copy(*view_device_,*view_host_unmanaged_);
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
      Kokkos::deep_copy(e_space,*view_device_,*view_host_unmanaged_);
      device_view_is_synced_ = true;
    }

    /// Returns a host mirror view for the outer view, where the inner
    /// views are still on device.
    auto getViewHost() const
    {
      TEUCHOS_ASSERT(is_initialized_);
      return *view_host_;
    }

    /// Returns device view of views
    auto getViewDevice() const
    {
      KOKKOS_ASSERT(device_view_is_synced_);
      return *view_device_;
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
