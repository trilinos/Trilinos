#include <impl/Kokkos_Timer.hpp>
#include <iostream>
#include <algorithm>

template< class Scalar, class DeviceType >
struct TrivialForFunctor;

template< class Scalar >
struct TrivialForFunctor< Scalar, KOKKOS_MACRO_DEVICE > {
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type ;
  typedef Kokkos::ValueView<Scalar, device_type > value_type ;

  value_type X , Y , Z;

  TrivialForFunctor( 
    const value_type & argX , 
    const value_type & argY , 
    const value_type & argZ )
  : X( argX ) , Y( argY ) , Z( argZ ) {}


  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type iwork ) const
  { *Z = *X + *Y; }
};


template<class Scalar, class Device>
void test_run(const char *test_name, unsigned int runs){
  if(runs == 0){
    return;
  }
  std::cout << test_name << std::endl;

  double *launch_times = new double[runs];
  double *copy_times = new double[runs];

  typedef Kokkos::ValueView<Scalar, Device> value_type;
  value_type x = Kokkos::create_value<Scalar, Device>();
  value_type y = Kokkos::create_value<Scalar, Device>();
  value_type z = Kokkos::create_value<Scalar, Device>();
  Kokkos::deep_copy(x, Scalar(4));
  Kokkos::deep_copy(y, Scalar(5));
  for(unsigned int i= 0; i< runs; ++i){
    Kokkos::deep_copy(z, Scalar(0));  
    Device::wait_functor_completion();
    Kokkos::Impl::Timer wall_clock ;
    Kokkos::parallel_for(1, TrivialForFunctor<Scalar, Device>(x, y, z));
    Device::wait_functor_completion();
    launch_times[i] = wall_clock.seconds();
  }

  std::sort(launch_times, launch_times+runs);
  double median = launch_times[runs/2];
  double max_time=launch_times[runs-1];
  double min_time=launch_times[0];
  std::cout << "Kernel Launch Latency" << std::endl;
  std::cout << "Min: " << min_time << std::endl; 
  std::cout << "Max: " << max_time << std::endl; 
  std::cout << "Median: " << median << std::endl; 
  std::cout << "Times: [";
  for(unsigned int i =0; i<runs; ++i){
    if(i == runs-1){
      std::cout << launch_times[i] << "]" << std::endl;
    }
    else{
      std::cout << launch_times[i] << ", ";
    }
  }


  double target = 0;
  for(unsigned int i =0; i<runs; ++i){
    Kokkos::Impl::Timer wall_clock ;
    Kokkos::deep_copy(target, x);
    copy_times[i] = wall_clock.seconds();
  }
  
  std::sort(copy_times, copy_times+runs);
  median = copy_times[runs/2];
  max_time=copy_times[runs-1];
  min_time=copy_times[0];
  std::cout << "Copy Latency" << std::endl;
  std::cout << "Min: " << min_time << std::endl; 
  std::cout << "Max: " << max_time << std::endl; 
  std::cout << "Median: " << median << std::endl; 
  std::cout << "Times: [";
  for(unsigned int i =0; i<runs; ++i){
    if(i == runs-1){
      std::cout << copy_times[i] << "]" << std::endl;
    }
    else{
      std::cout << copy_times[i] << ", ";
    }
  }

}



