/*
 * random_numbers.cpp
 *
 *  Created on: Jul 17, 2014
 *      Author: crtrott
 */

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_DualView.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <cstdlib>

typedef Kokkos::Impl::DefaultDeviceType::host_mirror_device_type DefaultHostType;

template<class GeneratorPool>
struct generate_random {
  GeneratorPool rand_pool;
  Kokkos::View<uint64_t*> vals;
  int samples;

  generate_random(Kokkos::View<uint64_t*> vals_,
                       GeneratorPool rand_pool_,
                       int samples_):
                       vals(vals_),rand_pool(rand_pool_),samples(samples_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    // Get a random number state from the pool
    typename GeneratorPool::generator_type rand_gen = rand_pool.get_state();

    // Draw samples numbers from the pool
    for(int k = 0;k<samples;k++)
      vals(i*samples+k) = rand_gen.urand64();


    rand_pool.free_state(rand_gen);
  }
};




int main(int argc, char* args[]) {
  srand(5374857);
  Kokkos::initialize(argc,args);
  int size = atoi(args[1]);
  int samples = atoi(args[2]);

  Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
  Kokkos::Random_XorShift1024_Pool<> rand_pool1024(5374857);
  Kokkos::DualView<uint64_t*> vals("Vals",size*samples);

  // Run some performance comparisons of the virtual and non virtual variant
  Kokkos::Impl::Timer timer;
  Kokkos::parallel_for(size,generate_random<Kokkos::Random_XorShift64_Pool<> >(vals.d_view,rand_pool64,samples));
  Kokkos::fence();

  timer.reset();
  Kokkos::parallel_for(size,generate_random<Kokkos::Random_XorShift64_Pool<> >(vals.d_view,rand_pool64,samples));
  Kokkos::fence();
  double time_64 = timer.seconds();

  Kokkos::parallel_for(size,generate_random<Kokkos::Random_XorShift1024_Pool<> >(vals.d_view,rand_pool1024,samples));
  Kokkos::fence();

  timer.reset();
  Kokkos::parallel_for(size,generate_random<Kokkos::Random_XorShift1024_Pool<> >(vals.d_view,rand_pool1024,samples));
  Kokkos::fence();
  double time_1024 = timer.seconds();

  printf("#Time XorShift64*:   %lf %lf\n",time_64,1.0e-9*samples*size/time_64 );
  printf("#Time XorShift1024*: %lf %lf\n",time_1024,1.0e-9*samples*size/time_1024 );

  Kokkos::deep_copy(vals.h_view,vals.d_view);

  Kokkos::finalize();
  return 0;
}


