// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Stokhos_Sacado_Kokkos.hpp"
#include "Kokkos_View.hpp"

// The function to compute the polynomial chaos expansion of,
// written as a template function
template <typename ScalarType, typename InputScalarType, typename OutputScalarType>
KOKKOS_INLINE_FUNCTION
void simple_function(const InputScalarType& x, OutputScalarType& y) {
  ScalarType u = x*x;
  ScalarType v = std::pow(std::log(u),2.0);
  y = v/(x + 1.0);
}

template <typename ScalarType, typename ScalarView >
KOKKOS_INLINE_FUNCTION
void simple_function(const ScalarView & x, const ScalarView & y)
{
  ScalarType u = x*x;
  ScalarType v = std::pow(std::log(u),2.0);
  y = v/(x + 1.0);
}

template < typename Scalar,
           typename ArrayVector,
           typename ScalarVector,
           typename Device >
struct vector_kernel;

template < typename Scalar,
           typename ArrayVector,
           typename ScalarVector,
           typename Device >
struct vector_kernel {
  typedef Scalar scalar_type;
  typedef ArrayVector array_vector_type;
  typedef ScalarVector scalar_vector_type;
  typedef Device device_type;
  typedef typename array_vector_type::storage_type storage_type;
  typedef Kokkos::LayoutRight layout_type;
  typedef Kokkos::View<scalar_type**, layout_type, device_type> view_type;

  view_type dev_x, dev_y;
  bool reset, print;

  KOKKOS_INLINE_FUNCTION
  void operator() (device_type device) const {
    int element = device.league_rank();
    int num_threads = device.team_size();
    int thread = device.team_rank();

    int num_samples = dev_x.dimension_1();
    int num_samples_per_thread = num_samples / num_threads;

    // multi-point expansions
    array_vector_type x(num_samples_per_thread, 0.0), y(num_samples_per_thread, 0.0);

    // Initialize x
    if (reset && storage_type::supports_reset) {
      storage_type& x_s = x.storage();
      storage_type& y_s = y.storage();
      x_s.shallowReset(&dev_x(element,thread),
                       num_samples_per_thread,
                       num_threads,
                       false);
      y_s.shallowReset(&dev_y(element,thread),
                       num_samples_per_thread,
                       num_threads,
                       false);
    }
    else {
      for (int sample=0; sample<num_samples_per_thread; ++sample)
        x.fastAccessCoeff(sample) = dev_x(element,thread+sample*num_threads);
    }

    simple_function<scalar_vector_type>(x,y);

    // Print x and y
    if (print) {
      for (int tidx = 0; tidx<num_threads; tidx++) {
        if (thread == tidx) {
          printf("x(%i) = [ ",tidx);
          for (int sample=0; sample<num_samples_per_thread; sample++)
            printf("%g ", x.coeff(sample));
          printf("]\n\n");
        }
        device.team_barrier();
      }

      for (int tidx = 0; tidx<num_threads; tidx++) {
        if (thread == tidx) {
          printf("y(%i) = [ ",tidx);
          for (int sample=0; sample<num_samples_per_thread; sample++)
            printf("%g ", y.coeff(sample));
          printf("]\n\n");
        }
        device.team_barrier();
      }
    }

    // Return result
    if (!(reset && storage_type::supports_reset)) {
      for (int sample=0; sample<num_samples_per_thread; ++sample)
        dev_y(element,thread+sample*num_threads) = y.fastAccessCoeff(sample);
    }
  }

  static void run(Kokkos::ParallelWorkRequest config,
                  const view_type& x, const view_type& y,
                  bool reset, bool print) {
    vector_kernel kernel;
    kernel.dev_x = x;
    kernel.dev_y = y;
    kernel.reset = reset;
    kernel.print = print;

    Kokkos::parallel_for(config, kernel);
    device_type::fence();
  }

};

template < typename Scalar,
           typename Ordinal,
           typename ScalarVector,
           typename Device >
struct vector_kernel< Scalar,
                      Sacado::MP::Vector< Stokhos::ViewStridedStorage<Ordinal,Scalar,Device> >,
                      ScalarVector,
                      Device > {
  typedef Scalar scalar_type;
  typedef Sacado::MP::Vector< Stokhos::ViewStridedStorage<Ordinal,Scalar,Device> > array_vector_type;
  typedef ScalarVector scalar_vector_type;
  typedef Device device_type;
  typedef typename array_vector_type::storage_type storage_type;
  typedef Kokkos::LayoutRight layout_type;
  typedef Kokkos::View<scalar_type**, layout_type, device_type> view_type;

  view_type dev_x, dev_y;
  bool reset, print;

  KOKKOS_INLINE_FUNCTION
  void operator() (device_type device) const {
    int element = device.league_rank();
    int num_threads = device.team_size();
    int thread = device.team_rank();

    int num_samples = dev_x.dimension_1();
    int num_samples_per_thread = num_samples / num_threads;

    // Initialize x
    storage_type x_s(&dev_x(element,thread),
                     num_samples_per_thread,
                     num_threads);
    storage_type y_s(&dev_y(element,thread),
                     num_samples_per_thread,
                     num_threads);
    array_vector_type x(x_s), y(y_s);

    simple_function<scalar_vector_type>(x,y);

    // Print x and y
    if (print) {
      for (int tidx = 0; tidx<num_threads; tidx++) {
        if (thread == tidx) {
          printf("x(%i) = [ ",tidx);
          for (int sample=0; sample<num_samples_per_thread; sample++)
            printf("%g ", x.coeff(sample));
          printf("]\n\n");
        }
        device.team_barrier();
      }

      for (int tidx = 0; tidx<num_threads; tidx++) {
        if (thread == tidx) {
          printf("y(%i) = [ ",tidx);
          for (int sample=0; sample<num_samples_per_thread; sample++)
            printf("%g ", y.coeff(sample));
          printf("]\n\n");
        }
        device.team_barrier();
      }
    }
  }

  static void run(Kokkos::ParallelWorkRequest config,
                  const view_type& x, const view_type& y,
                  bool reset, bool print) {
    vector_kernel kernel;
    kernel.dev_x = x;
    kernel.dev_y = y;
    kernel.reset = reset;
    kernel.print = print;

    Kokkos::parallel_for(config, kernel);
    device_type::fence();
  }

};


template < typename ViewType >
struct view_kernel
{
  typedef typename Kokkos::Impl::if_c< Kokkos::is_view< ViewType >::value , ViewType , void >::type view_type ;
  typedef typename view_type::value_type   value_type ;  ///<  == Sacado::MP::Vector< Storage >
  typedef typename view_type::scalar_type  scalar_type ; ///<  == Storage::value_type
  typedef typename view_type::device_type  device_type ; ///<  == Storage::device_type

  typedef Kokkos::View< typename view_type::data_type ,
                        typename view_type::array_layout ,
                        device_type ,
                        Kokkos::MemoryUnmanaged > local_view_type ;

  view_type dev_x, dev_y;
  bool reset, print;

  KOKKOS_INLINE_FUNCTION
  void operator() (device_type device) const
  {
    typename local_view_type::Partition part( device.team_rank() , device.team_size() );
    const local_view_type local_x( dev_x , part );
    const local_view_type local_y( dev_y , part );
    const int element = device.league_rank();

    // Apply evaluation function to this thread's fix-sized UQ sample set.
    simple_function< value_type >( local_x(element) , local_y(element) );

    // Print x and y
    if (print) {

      device.team_barrier();

      if ( ! device.league_rank() && ! device.team_rank() ) {
        printf("view_kernel league(%d:%d) team_size(%d) dim(%d) size(%d)\n",
               device.league_rank(), device.league_size(),device.team_size(),
               int(dev_x.dimension_1()), int(local_x(element).size()) );
      }

      if ( ! device.team_rank() ) {
        printf("x(%i) = { ",element);
        for (int sample=0; sample< int(dev_x.dimension_1()); ++sample) {
          // printf("%g ", dev_x(element,sample));
          printf("%g ", dev_x(element).coeff(sample));
        }
        printf("}\n\n");
        printf("y(%i) = { ",element);
        for (int sample=0; sample< int(dev_y.dimension_1()); ++sample) {
          // printf("%g ", dev_y(element,sample));
          printf("%g ", dev_y(element).coeff(sample));
        }
        printf("}\n\n");
      }

      device.team_barrier();
    }
  }

  static void run(Kokkos::ParallelWorkRequest config,
                  const view_type& x, const view_type& y,
                  bool reset, bool print) {
    view_kernel kernel;
    kernel.dev_x = x;
    kernel.dev_y = y;
    kernel.reset = reset;
    kernel.print = print;

    Kokkos::parallel_for(config, kernel);
    device_type::fence();
  }

};

template <typename Scalar, typename Device>
struct scalar_kernel {
  typedef Scalar scalar_type;
  typedef Device device_type;
  typedef Kokkos::LayoutRight layout_type;
  typedef Kokkos::View<scalar_type**, layout_type, device_type> view_type;

  view_type dev_x, dev_y;
  bool reset, print;

  KOKKOS_INLINE_FUNCTION
  void operator() (device_type device) const {
    int element = device.league_rank();
    int num_threads = device.team_size();
    int thread = device.team_rank();

    int num_samples = dev_x.dimension_1();

    scalar_type x, y;
    for (int sample=thread; sample<num_samples; sample+=num_threads) {

      // Initialize x
      x = dev_x(element, sample);

      // Compute function
      simple_function<scalar_type>(x,y);

      // Return result
      dev_y(element, sample) = y;

    }
  }

  static void run(Kokkos::ParallelWorkRequest config,
                  const view_type& x, const view_type& y,
                  bool reset, bool print) {
    scalar_kernel kernel;
    kernel.dev_x = x;
    kernel.dev_y = y;
    kernel.reset = reset;
    kernel.print = print;

    Kokkos::parallel_for(config, kernel);
    device_type::fence();
  }
};

template <typename Scalar, typename ArrayVector, typename ScalarVector, typename Device>
bool run_kernels(Kokkos::ParallelWorkRequest config,
                 int num_elements, int num_samples, bool reset, bool print,
                 const std::string& device_name) {
  typedef vector_kernel<Scalar,ArrayVector,ScalarVector,Device> vec_kernel;
  typedef scalar_kernel<Scalar,Device> sca_kernel;
  typedef typename vec_kernel::view_type view_type;
  typedef typename view_type::HostMirror host_view_type;

  view_type x     = view_type("x", num_elements, num_samples);
  view_type y_vec = view_type("y", num_elements, num_samples);
  view_type y_sca = view_type("y", num_elements, num_samples);

  host_view_type hx = Kokkos::create_mirror(x);

  // Initialize x
  for (int element=0; element<num_elements; element++) {
    for (int sample=0; sample<num_samples; sample++) {
      hx(element,sample) =
        static_cast<Scalar>(element+sample+1) /
        static_cast<Scalar>(num_elements*num_samples);
    }
  }

  // Copy x to device
  Kokkos::deep_copy(x, hx);

  // Run vector and scalar kernels
  {
    TEUCHOS_FUNC_TIME_MONITOR(device_name + " calculation");
    sca_kernel::run(config, x, y_sca, reset, print);
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR(device_name + " calculation (MP)");
    vec_kernel::run(config, x, y_vec, reset, print);
  }

  // Copy results back to host
  host_view_type hy_vec = Kokkos::create_mirror(y_vec);
  host_view_type hy_sca = Kokkos::create_mirror(y_sca);
  Kokkos::deep_copy(hy_vec, y_vec);
  Kokkos::deep_copy(hy_sca, y_sca);

  // Check results agree
  double rtol = 1e-15;
  double atol = 1e-15;
  bool agree = true;
  for (int e=0; e<num_elements; e++) {
    for (int s=0; s<num_samples; s++) {
      if (std::abs(hy_vec(e,s)-hy_sca(e,s)) > std::abs(hy_sca(e,s))*rtol+atol) {
        agree = false;
        break;
      }
    }
  }

  // Print results if requested
  if (print) {
    std::cout << "x    = [ ";
    for (int e=0; e<num_elements; e++) {
      for (int s=0; s<num_samples; s++)
        std::cout << hx(e,s) << " ";
      std::cout << ";" << std::endl;
    }
    std::cout << "]" << std::endl;

    std::cout << "y      [ ";
    for (int e=0; e<num_elements; e++) {
      for (int s=0; s<num_samples; s++)
        std::cout << hy_sca(e,s) << " ";
      std::cout << ";" << std::endl;
    }
    std::cout << "]" << std::endl;

    std::cout << "y_mp = [ ";
    for (int e=0; e<num_elements; e++) {
      for (int s=0; s<num_samples; s++)
        std::cout << hy_vec(e,s) << " ";
      std::cout << ";" << std::endl;
    }
    std::cout << "]" << std::endl;
  }

  return agree;
}

template < typename Scalar , unsigned SamplesPerThread , typename Device >
bool run_view_kernel( Kokkos::ParallelWorkRequest config,
                      int num_elements ,
                      int num_samples ,
                      bool reset, bool print,
                      const std::string& device_name )
{
  typedef typename
    Kokkos::Impl::if_c< SamplesPerThread ,
                        Stokhos::StaticFixedStorage< int , Scalar , SamplesPerThread , Device > ,
                        Stokhos::DynamicStorage<     int , Scalar , Device > >::type
      storage_type ;

  typedef Kokkos::View< Sacado::MP::Vector< storage_type > * , Device > view_type ;

  typedef view_kernel< view_type > kernel;
  typedef typename view_type::HostMirror host_view_type;

  const bool input_error = SamplesPerThread && ( num_samples != int(config.team_size * SamplesPerThread) );

  if ( print || input_error ) {
    std::cout << "run_view_kernel" ;
    if ( input_error ) {
      std::cout << " : ERROR num_samples != static(SamplesPerThread) * team_size : " ;
    }
    std::cout << " num_elements(" << num_elements << ")"
              << " num_samples(" << num_samples << ")"
              << " league_size(" << config.league_size << ")"
              << " team_size(" << config.team_size << ")" ;
    if ( SamplesPerThread ) { std::cout << " static(" << SamplesPerThread << ")" ; }
    else { std::cout << " dynamic" ; }
    std::cout << std::endl ;

    if ( input_error ) { return false ; }
  }

  view_type x( "x", num_elements, num_samples);
  view_type y( "y", num_elements, num_samples);

  host_view_type hy_ans("y_ans", num_elements, num_samples);

  host_view_type hx = Kokkos::create_mirror(x);

  typename host_view_type::array_type hay_ans = hy_ans ;
  typename host_view_type::array_type hax = hx ;

  // Initialize x
  for (int element=0; element<num_elements; element++) {
    for (int sample=0; sample<num_samples; sample++) {
      hax(element,sample) =
        static_cast<Scalar>(element+sample+1) /
        static_cast<Scalar>(num_elements*num_samples);
      simple_function<Scalar>( hax(element,sample) , hay_ans(element,sample) );
    }
  }

  // Copy x to device
  Kokkos::deep_copy(x, hx);

  // Run vector and scalar kernels
  {
    TEUCHOS_FUNC_TIME_MONITOR(device_name + " calculation");
    kernel::run(config, x, y, reset, print);
  }

  // Copy results back to host
  host_view_type hy = Kokkos::create_mirror(y);

  Kokkos::deep_copy(hy, y);

  typename host_view_type::array_type hay = hy ;

  // Check results agree
  double rtol = 1e-15;
  double atol = 1e-15;
  bool agree = true;
  for (int e=0; e<num_elements; e++) {
    for (int s=0; s<num_samples; s++) {
      if (std::abs(hay(e,s)-hay_ans(e,s)) > std::abs(hay_ans(e,s))*rtol+atol) {
        agree = false;
        break;
      }
    }
  }

  // Print results if requested
  if (print) {
    std::cout << "x = [ ";
    for (int e=0; e<num_elements; e++) {
      for (int s=0; s<num_samples; s++)
        std::cout << hax(e,s) << " ";
      std::cout << ";" << std::endl;
    }
    std::cout << "]" << std::endl;

    std::cout << "y = [ ";
    for (int e=0; e<num_elements; e++) {
      for (int s=0; s<num_samples; s++)
        std::cout << hay(e,s) << " ";
      std::cout << ";" << std::endl;
    }
    std::cout << "]" << std::endl;

    std::cout << "hy_ans = [ ";
    for (int e=0; e<num_elements; e++) {
      for (int s=0; s<num_samples; s++)
        std::cout << hay_ans(e,s) << " ";
      std::cout << ";" << std::endl;
    }
    std::cout << "]" << std::endl;
  }

  return agree;
}

// storage options
enum Storage_Method { STATIC,
                      STATIC_FIXED,
                      LOCAL,
                      DYNAMIC,
                      DYNAMIC_STRIDED,
                      DYNAMIC_THREADED,
                      VIEW_STATIC,
                      VIEW_STATIC_FIXED,
                      VIEW_LOCAL,
                      VIEW_DYNAMIC,
                      VIEW_DYNAMIC_STRIDED,
                      VIEW_DYNAMIC_THREADED };

template <int MaxSize, typename scalar_type, typename device_type>
struct MPVectorTypes {
  // Storage types
  typedef Stokhos::StaticStorage<int,scalar_type,MaxSize,device_type> static_storage;
  typedef Stokhos::StaticFixedStorage<int,scalar_type,1,device_type> static_fixed_storage_1;
  typedef Stokhos::StaticFixedStorage<int,scalar_type,2,device_type> static_fixed_storage_2;
  typedef Stokhos::StaticFixedStorage<int,scalar_type,4,device_type> static_fixed_storage_4;
  typedef Stokhos::StaticFixedStorage<int,scalar_type,8,device_type> static_fixed_storage_8;
  typedef Stokhos::LocalStorage<int,scalar_type,1,device_type> local_storage_1;
  typedef Stokhos::LocalStorage<int,scalar_type,2,device_type> local_storage_2;
  typedef Stokhos::LocalStorage<int,scalar_type,4,device_type> local_storage_4;
  typedef Stokhos::LocalStorage<int,scalar_type,8,device_type> local_storage_8;
  typedef Stokhos::DynamicStorage<int,scalar_type,device_type> dynamic_storage;
  typedef Stokhos::DynamicStridedStorage<int,scalar_type,device_type> dynamic_strided_storage;
  typedef Stokhos::DynamicThreadedStorage<int,scalar_type,device_type> dynamic_threaded_storage;
  typedef Stokhos::ViewStridedStorage<int,scalar_type,device_type> view_storage;

  // Vector types
  typedef Sacado::MP::Vector<static_storage> static_vector;
  typedef Sacado::MP::Vector<static_fixed_storage_1> static_fixed_vector_1;
  typedef Sacado::MP::Vector<static_fixed_storage_2> static_fixed_vector_2;
  typedef Sacado::MP::Vector<static_fixed_storage_4> static_fixed_vector_4;
  typedef Sacado::MP::Vector<static_fixed_storage_8> static_fixed_vector_8;
  typedef Sacado::MP::Vector<local_storage_1> local_vector_1;
  typedef Sacado::MP::Vector<local_storage_2> local_vector_2;
  typedef Sacado::MP::Vector<local_storage_4> local_vector_4;
  typedef Sacado::MP::Vector<local_storage_8> local_vector_8;
  typedef Sacado::MP::Vector<dynamic_storage> dynamic_vector;
  typedef Sacado::MP::Vector<dynamic_strided_storage> dynamic_strided_vector;
  typedef Sacado::MP::Vector<dynamic_threaded_storage> dynamic_threaded_vector;
  typedef Sacado::MP::Vector<view_storage> view_vector;

  // View to array of vector types
  typedef Kokkos::View< static_fixed_vector_1 * , Kokkos::LayoutLeft , device_type > static_fixed_view_left_1 ;
  typedef Kokkos::View< static_fixed_vector_2 * , Kokkos::LayoutLeft , device_type > static_fixed_view_left_2 ;
  typedef Kokkos::View< static_fixed_vector_4 * , Kokkos::LayoutLeft , device_type > static_fixed_view_left_4 ;
  typedef Kokkos::View< static_fixed_vector_8 * , Kokkos::LayoutLeft , device_type > static_fixed_view_left_8 ;

  typedef Kokkos::View< static_fixed_vector_1 * , Kokkos::LayoutRight , device_type > static_fixed_view_right_1 ;
  typedef Kokkos::View< static_fixed_vector_2 * , Kokkos::LayoutRight , device_type > static_fixed_view_right_2 ;
  typedef Kokkos::View< static_fixed_vector_4 * , Kokkos::LayoutRight , device_type > static_fixed_view_right_4 ;
  typedef Kokkos::View< static_fixed_vector_8 * , Kokkos::LayoutRight , device_type > static_fixed_view_right_8 ;
};

template <int MaxSize, typename Scalar, typename Device>
struct MPVectorExample {
  static bool run(Storage_Method storage_method,
                  int num_elements, int num_samples,
                  int team_size, int league_size, bool reset, bool print);
};

// Maximum size of expansion -- currently 2, 4, or 8 for LocalStorage
const int MaxSize = 8;

// Scalar type used in kernels
typedef double Scalar;

