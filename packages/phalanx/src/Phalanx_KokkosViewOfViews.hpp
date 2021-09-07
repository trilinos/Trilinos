#ifndef PHALANX_KOKKOS_VIEW_OF_VIEWS_HPP
#define PHALANX_KOKKOS_VIEW_OF_VIEWS_HPP

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

  /** Wrapper class that correctly handles ViewOfViews construction
      and object lifetime. Can be used for more than just Views as the
      inner object. This class makes sure the host view stays in scope
      for the life of the device view and makes sure that the device
      is synced to host before use.

      Main restrictions:

      1. When UVM is not used in the outer view, we need to allocate
      the VofV on host and copy to device to set up the inner views
      correctly.

      2. Step 1 means that the host view must exist as long as the
      device view is being used, otherwise the views may go out of
      scope. This object exists to pair up the host and device view to
      make sure the inner views are not deleted early.
  */
  template<int OuterViewRank,typename InnerViewType,typename MemorySpace>
  class ViewOfViews {

  public:
    // Layout of the outer view doesn't matter for performance so we
    // use a raw Kokkos::View instead of PHX::View. The inner views are
    // what is important for performance.
    using OuterDataType = typename PHX::v_of_v_utils::add_pointer<InnerViewType,OuterViewRank>::type;
    using OuterViewType = Kokkos::View<OuterDataType,MemorySpace>;

  private:
    typename OuterViewType::HostMirror view_host_;
    OuterViewType view_device_;
    bool device_view_is_synced_;

  public:
    template<typename... Extents>
    ViewOfViews(const std::string name,
                Extents... extents)
      : view_host_(name,extents...),
        view_device_(name,extents...),
        device_view_is_synced_(false)
    {}

    ~ViewOfViews()
    {
      // Make sure there is not another object pointing to device view
      // since the host view will delete the inner views on exit.
      if (view_device_.impl_track().use_count() != 1)
        Kokkos::abort("\n ERROR - PHX::ViewOfViews - please free all instances of device ViewOfView \n before deleting the host ViewOfView!\n\n");
    }

    template<typename... Indices>
    void addView(InnerViewType v,Indices... i)
    {
      view_host_(i...) = v;
      device_view_is_synced_ = false;
    }

    void syncHostToDevice()
    {
      Kokkos::deep_copy(view_device_,view_host_);
      device_view_is_synced_ = true;
    }

    auto getViewHost()
    {
      return view_host_;
    }

    auto getViewDevice()
    {
      TEUCHOS_ASSERT(device_view_is_synced_);
      return view_device_;
    }
  };

}

#endif
