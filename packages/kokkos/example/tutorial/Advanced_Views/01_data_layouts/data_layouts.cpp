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
#include <impl/Kokkos_Timer.hpp>
#include <cstdio>

typedef Kokkos::View<double**, Kokkos::LayoutLeft> left_type;
typedef Kokkos::View<double**, Kokkos::LayoutRight> right_type;
typedef Kokkos::View<double*> view_type;

template<class ViewType>
struct init_view {
  ViewType a;
  init_view(ViewType a_):a(a_) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < a.dimension_1(); j++)
      a(i,j) = 1.0*a.dimension_0()*i + 1.0*j;
  }
};

template<class ViewType1, class ViewType2>
struct contraction {
  view_type a;
  typename ViewType1::const_type v1;
  typename ViewType2::const_type v2;
  contraction(view_type a_, ViewType1 v1_, 
              ViewType2 v2_):a(a_),v1(v1_),v2(v2_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < v1.dimension_1(); j++)
      a(i) = v1(i,j)*v2(j,i);
  }
};

struct dot {
  view_type a;
  dot(view_type a_):a(a_) {};
  typedef double value_type; //Specify type for reduction target, lsum
  KOKKOS_INLINE_FUNCTION
  void operator() (int i, double &lsum) const {
    lsum+= a(i)*a(i);
  }
};
 
int main(int narg, char* arg[]) {
  Kokkos::initialize(narg,arg);
  
  int size = 10000;
  view_type a("A",size);
  left_type l("L",size,10000);
  right_type r("R",size,10000);


  Kokkos::parallel_for(size,init_view<left_type>(l));
  Kokkos::parallel_for(size,init_view<right_type>(r));
  Kokkos::fence();

  Kokkos::Impl::Timer time1;
  Kokkos::parallel_for(size,contraction<left_type,right_type>(a,l,r));
  Kokkos::fence();  
  double sec1 = time1.seconds();

  double sum1 = 0;
  Kokkos::parallel_reduce(size,dot(a),sum1);
  Kokkos::fence();
  
  Kokkos::Impl::Timer time2;
  Kokkos::parallel_for(size,contraction<right_type,left_type>(a,r,l));
  Kokkos::fence();
  double sec2 = time2.seconds();

  double sum2 = 0;
  Kokkos::parallel_reduce(size,dot(a),sum2);


  printf("Result Left/Rigth %lf Right/Left %lf  (equal result: %i)\n",sec1,sec2,sum2==sum1);  

  Kokkos::finalize();
}

