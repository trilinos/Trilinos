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
#include <cstdio>
#include <cstdlib>
#include <cmath>

typedef Kokkos::View<int*> view_type;
typedef view_type::HostMirror host_view_type;
typedef Kokkos::View<int> count_type;
typedef count_type::HostMirror host_count_type;

struct findprimes {
  view_type data;
  view_type result;

  count_type count;
  findprimes(view_type data_, view_type result_, count_type count_):data(data_),result(result_),count(count_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    const int number = data(i);
    const int upper_bound = sqrt(1.0*number)+1;
    bool is_prime = !(number%2 == 0);
    int k = 3;
    while(k<upper_bound && is_prime) {
      is_prime = !(number%k == 0);
      k+=2;
    }

    if(is_prime) {
      int idx = Kokkos::atomic_fetch_add(&count(),1);
      result(idx) = number;
    }
  }

};
 
int main() {
  Kokkos::initialize();

  srand(61391);

  int nnumbers = 100000;
  view_type data("RND",nnumbers);
  view_type result("Prime",nnumbers);
  count_type count("Count");

  host_view_type h_data = Kokkos::create_mirror_view(data);
  host_view_type h_result = Kokkos::create_mirror_view(result);
  host_count_type h_count = Kokkos::create_mirror_view(count);
  
  for(int i = 0; i < data.dimension_0(); i++)
      h_data(i) = rand()%100000;

  Kokkos::deep_copy(data,h_data);

  int sum = 0;
  Kokkos::parallel_for(data.dimension_0(),findprimes(data,result,count));
  Kokkos::deep_copy(h_count,count);

  printf("Found %i prime numbers in %i random numbers\n",h_count(),nnumbers);
  Kokkos::finalize();
}

