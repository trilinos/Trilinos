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

typedef Kokkos::View<double*[3]> view_type;
typedef view_type::HostMirror host_view_type;

struct squaresum {
  view_type a;
  squaresum(view_type a_):a(a_) {}
  typedef int value_type; //Specify type for reduction value, lsum

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int &lsum) const {
    lsum += a(i,0)-a(i,1)+a(i,2);
  }
};
 
int main() {
  Kokkos::initialize();

  view_type a("A",10);
  host_view_type h_a = Kokkos::create_mirror_view(a); 
  
  for(int i = 0; i < 10; i++)
    for(int j = 0; j < 3; j++)
      h_a(i,j) = i*10 + j;

  Kokkos::deep_copy(a,h_a);

  int sum = 0;
  Kokkos::parallel_reduce(10,squaresum(a),sum);
  printf("Result is %i\n",sum);  

  Kokkos::finalize();
}

