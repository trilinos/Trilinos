template<typename Scalar , class DeviceType>
struct minimum_stable_time_step;

template<typename Scalar>
struct minimum_stable_time_step<Scalar ,KOKKOS_MACRO_DEVICE>{

  typedef KOKKOS_MACRO_DEVICE       device_type;
  typedef device_type::size_type    size_type;

  typedef Scalar value_type;

  typedef Region<Scalar,device_type> MyRegion;

    minimum_stable_time_step( const MyRegion  & arg_region)
       : region(arg_region)
      {}


    KOKKOS_MACRO_DEVICE_FUNCTION
    static void init(value_type &update) {
      update = 1.0e32;
    }

    KOKKOS_MACRO_DEVICE_FUNCTION
    static void join(volatile value_type &update, const volatile value_type & source) {
      update = update < source ? update : source;
    }


    KOKKOS_MACRO_DEVICE_FUNCTION
    void operator()(int ielem, value_type & update) const {
      value_type source = region.elem_t_step(ielem);
      update = update < source ? update : source;
    }

    MyRegion   region;

}; //minimum_stable_time_step




template<typename Scalar , class DeviceType>
struct set_next_time_step;

template<typename Scalar>
struct set_next_time_step<Scalar ,KOKKOS_MACRO_DEVICE>{

  typedef KOKKOS_MACRO_DEVICE       device_type;
  typedef device_type::size_type    size_type;

  typedef Scalar value_type;

  typedef Region<Scalar,device_type> MyRegion;

    set_next_time_step(
                const MyRegion  & arg_region,
                const int       arg_next_state)
       : region(arg_region)
       , next_state(arg_next_state)
      {}


    KOKKOS_MACRO_DEVICE_FUNCTION
    void operator()(Scalar & result) const {
      region.delta_t(next_state) = result;
    }

    MyRegion   region;
    const int  next_state;

}; //minimum_stable_time_step

