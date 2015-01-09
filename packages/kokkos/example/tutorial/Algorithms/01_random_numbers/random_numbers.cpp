/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_DualView.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <cstdlib>

typedef Kokkos::HostSpace::execution_space DefaultHostType;

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
  if (argc != 3){
	printf("Please pass two integers on the command line\n");
  }
  else {
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
  }
  return 0;
}


