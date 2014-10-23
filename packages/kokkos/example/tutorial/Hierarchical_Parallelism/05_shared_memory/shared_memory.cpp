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
#include <Kokkos_DualView.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <cstdio>
#include <cstdlib>

typedef Kokkos::DefaultExecutionSpace         ExecutionSpace ;
typedef Kokkos::HostSpace::execution_space    HostExecutionSpace ;
typedef Kokkos::TeamPolicy< ExecutionSpace >  team_policy ;
typedef team_policy::member_type              team_member ;

#define TEAM_SIZE 16

struct find_2_tuples {
  int chunk_size;
  Kokkos::View<const int*> data;
  Kokkos::View<int**> histogram;

  typedef ExecutionSpace::scratch_memory_space SharedSpace;

  typedef Kokkos::View<int**,SharedSpace,Kokkos::MemoryUnmanaged> shared_2d_int;
  typedef Kokkos::View<int*,SharedSpace,Kokkos::MemoryUnmanaged> shared_1d_int;


  find_2_tuples(int chunk_size_, Kokkos::DualView<int*> data_,
                Kokkos::DualView<int**> histogram_):chunk_size(chunk_size_),
                data(data_.d_view),histogram(histogram_.d_view) {
      data_.sync<ExecutionSpace>();
      histogram_.sync<ExecutionSpace>();
      histogram_.modify<ExecutionSpace>();
  }

  KOKKOS_INLINE_FUNCTION
  void operator() ( const team_member & thread) const {
    // FIXME (mfh 23 Oct 2014) It seems like we should use
    // thread.team_size() here, instead of TEAM_SIZE.  However, the
    // example still fails in that case.  Not sure what to do.

    shared_2d_int l_histogram(thread.team_shmem(),TEAM_SIZE,TEAM_SIZE);
    shared_1d_int l_data(thread.team_shmem(),chunk_size+1);

    const int i = thread.league_rank() * chunk_size;
    thread.team_par_for( chunk_size , [&] (int& j) {
      l_data(j) = data(i+j);
    });

    thread.team_par_for( chunk_size , [&] (int& k) {
      for(int l = 0; l < TEAM_SIZE; l++)
        l_histogram(k,l) = 0;
    });

    thread.team_barrier();

    for(int j = 0; j<chunk_size; j++) {
      thread.team_par_for( chunk_size , [&] (int& k) {
        for(int l = 0; l < TEAM_SIZE; l++) {
          if((l_data(j) == k) && (l_data(j+1)==l))
            l_histogram(k,l)++;
        }
      });
    }

    thread.team_par_for( chunk_size , [&] (int& k) {
      for(int l = 0; l < TEAM_SIZE; l++) {
        Kokkos::atomic_fetch_add(&histogram(k,l),l_histogram(k,l));
      }
    });
    thread.team_barrier();
  }

  size_t team_shmem_size( int team_size ) const {
    return shared_2d_int::shmem_size(team_size,team_size) +
           shared_1d_int::shmem_size(chunk_size+1);
  }
};

int main(int narg, char* args[]) {
  Kokkos::initialize(narg,args);

  int chunk_size = 1024;
  int nchunks = 100000; //1024*1024;
  Kokkos::DualView<int*> data("data",nchunks*chunk_size+1);

  srand(1231093);

  for(int i = 0; i < data.dimension_0(); i++) {
    data.h_view(i) = rand()%TEAM_SIZE;
  }
  data.modify<Kokkos::HostSpace>();
  data.sync<ExecutionSpace>();

  Kokkos::DualView<int**> histogram("histogram",TEAM_SIZE,TEAM_SIZE);


  Kokkos::Impl::Timer timer;
  // Threads/team (TEAM_SIZE) is automically limited to the maximum supported by the device.
  Kokkos::parallel_for( team_policy( nchunks , TEAM_SIZE )
                      , find_2_tuples(chunk_size,data,histogram) );
  Kokkos::fence();
  double time = timer.seconds();

  histogram.sync<Kokkos::HostSpace>();

  printf("Time: %lf \n\n",time);
  int sum = 0;
  for(int k=0; k<TEAM_SIZE; k++) {
    for(int l=0; l<TEAM_SIZE; l++) {
      printf("%i ",histogram.h_view(k,l));
      sum += histogram.h_view(k,l);
    }
    printf("\n");
  }
  printf("Result: %i %i\n",sum,chunk_size*nchunks);
  Kokkos::finalize();
}

