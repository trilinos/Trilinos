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
#include <cstdio>

#ifdef KOKKOS_HAVE_CXX11

// Using default execution space:
typedef typename Kokkos::TeamPolicy<>::member_type team_member ;

struct SomeCorrelation {
  typedef int value_type; //Specify value type for reduction target, sum
  typedef Kokkos::DefaultExecutionSpace::scratch_memory_space shared_space;
  typedef Kokkos::View<int*,shared_space,Kokkos::MemoryUnmanaged> shared_1d_int;

  Kokkos::View<const int***,Kokkos::LayoutRight> data;
  Kokkos::View<int> gsum;

  SomeCorrelation(Kokkos::View<int***,Kokkos::LayoutRight> data_in,
                  Kokkos::View<int> sum):data(data_in),gsum(sum){}

  KOKKOS_INLINE_FUNCTION
  void operator() ( const team_member & thread) const {
    int i = thread.league_rank();

    shared_1d_int count(thread.team_shmem(),data.dimension_1());

    // With each team run a parallel_for with its threads
    Kokkos::parallel_for(Kokkos::TeamThreadLoop(thread,data.dimension_1()), [=] (const int& j) {
      int tsum;
      // Run a vector loop reduction over the inner dimension of data
      // Count how many values are multiples of 4
      // Every vector lane gets the same reduction value (tsum) back
      Kokkos::parallel_reduce(Kokkos::ThreadVectorLoop(thread,data.dimension_2()), [=] (const int& k, int & vsum) {
        vsum+= (data(i,j,k) % 4 == 0)?1:0;
      },tsum);

      // Make sure only one vector lane adds the reduction value to the shared array
      Kokkos::single(Kokkos::PerThread(thread),[=] () {
        count(j) = tsum;
      });
    });

    // Wait for all threads to finish the parallel_for so that all shared memory writes are done
    thread.team_barrier();

    // Check with one vector lane from each thread how many consecutive
    // data segments have the same number of values divisible by 4
    int team_sum = 0;
    Kokkos::parallel_reduce(Kokkos::TeamThreadLoop(thread, data.dimension_1()-1), [=] (const int& j, int& thread_sum) {
      // It is not valid to directly add to thread_sum
      // Use a single function with broadcast instead
      // team_sum will be used as input to the operator (i.e. it is used to initialize sum)
      // the end value of sum will be broadcast to all vector lanes in the thread.
      Kokkos::single(Kokkos::PerThread(thread),[=] (int& sum) {
        if(count(j)==count(j+1)) sum++;
      },thread_sum);
    },team_sum);

    // Add with one thread of the team the team_sum to the global value
    Kokkos::single(Kokkos::PerTeam(thread),[=] () {
      Kokkos::atomic_add(&gsum(),team_sum);
    });
  }

  size_t team_shmem_size( int team_size ) const {
    return shared_1d_int::shmem_size(data.dimension_1());
  }
};

int main(int narg, char* args[]) {
  Kokkos::initialize(narg,args);

  // Produce some 3D random data
  Kokkos::View<int***,Kokkos::LayoutRight> data("Data",512,512,32);
  Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
  Kokkos::fill_random(data,rand_pool64,100);

  // A global value to put the result in
  Kokkos::View<int> gsum("Sum");

  // Each team handles a slice of the data
  // Set up TeamPolicy with 512 teams with maximum number of threads per team and 16 vector lanes
  // The team_size_max function will determine the maximum number of threads taking into account
  // shared memory requirements of the Functor
  const Kokkos::TeamPolicy<> policy( 512 , Kokkos::TeamPolicy<>::team_size_max(SomeCorrelation(data,gsum)) , 16);

  Kokkos::parallel_for( policy , SomeCorrelation(data,gsum) );

  Kokkos::fence();

  // Copy result value back
  int sum = 0;
  Kokkos::deep_copy(sum,gsum);
  printf("Result %i\n",sum);

  Kokkos::finalize();
}

#endif //KOKKOS_HAVE_CXX11
