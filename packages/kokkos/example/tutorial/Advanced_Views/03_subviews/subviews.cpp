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

typedef Kokkos::View<double***, Kokkos::LayoutRight> mesh_type;
typedef Kokkos::View<double**, Kokkos::LayoutStride> xz_plane_type;
typedef Kokkos::View<double**, Kokkos::LayoutRight> yz_plane_type;
typedef Kokkos::View<double**, Kokkos::LayoutStride> xy_plane_type;
typedef Kokkos::View<double***, Kokkos::LayoutStride> inner_mesh_type;



template<class ViewType>
struct set_boundary {
  ViewType a;
  double value;
  set_boundary(ViewType a_, double value_):a(a_),value(value_) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < a.dimension_1(); j++)
      a(i,j) = value;
  }
};

template<class ViewType>
struct set_inner {
  ViewType a;
  double value;
  set_inner(ViewType a_, double value_):a(a_),value(value_) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for(int j = 0; j < a.dimension_1(); j++)
      for(int k = 0; k < a.dimension_2(); k++)
        a(i,j,k) = value;
  }
};

template<class ViewType>
struct update {
  ViewType a;
  double dt;
  update(ViewType a_, double dt_):a(a_),dt(dt_) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    i++;
    for(int j = 1; j < a.dimension_1()-1; j++)
      for(int k = 1; k < a.dimension_2()-1; k++)
        a(i,j,k) += dt*(a(i,j,k+1)-a(i,j,k-1) +
                        a(i,j+1,k)-a(i,j-1,k) +
                        a(i+1,j,k)-a(i-1,j,k));
  }
};


int main(int narg, char* arg[]) {
  Kokkos::initialize(narg,arg);
  
  int size = 100;
  mesh_type A("A",size+2,size+2,size+2);
  inner_mesh_type Ai = 
       Kokkos::subview<inner_mesh_type>
               (A,Kokkos::pair<int,int>(1,size+1),Kokkos::pair<int,int>(1,size+1),Kokkos::pair<int,int>(1,size+1));

  xy_plane_type Zneg_halo = 
       Kokkos::subview<xy_plane_type>
               (A,Kokkos::ALL(),Kokkos::ALL(),0);
  xy_plane_type Zpos_halo = 
       Kokkos::subview<xy_plane_type>
               (A,Kokkos::ALL(),Kokkos::ALL(),101);

  xz_plane_type Yneg_halo = 
       Kokkos::subview<xz_plane_type>
               (A,Kokkos::ALL(),0,Kokkos::ALL());
  xz_plane_type Ypos_halo =
       Kokkos::subview<xz_plane_type>
               (A,Kokkos::ALL(),101,Kokkos::ALL());

  yz_plane_type Xneg_halo = 
       Kokkos::subview<yz_plane_type>
               (A,0,Kokkos::ALL(),Kokkos::ALL());
  yz_plane_type Xpos_halo =
       Kokkos::subview<yz_plane_type>
               (A,101,Kokkos::ALL(),Kokkos::ALL());

  Kokkos::parallel_for(Zneg_halo.dimension_0(),set_boundary<xy_plane_type>(Zneg_halo,1));
  Kokkos::parallel_for(Zpos_halo.dimension_0(),set_boundary<xy_plane_type>(Zpos_halo,-1));
  Kokkos::parallel_for(Yneg_halo.dimension_0(),set_boundary<xz_plane_type>(Yneg_halo,2));
  Kokkos::parallel_for(Ypos_halo.dimension_0(),set_boundary<xz_plane_type>(Ypos_halo,-2));
  Kokkos::parallel_for(Xneg_halo.dimension_0(),set_boundary<yz_plane_type>(Xneg_halo,3));
  Kokkos::parallel_for(Xpos_halo.dimension_0(),set_boundary<yz_plane_type>(Xpos_halo,-3));
  Kokkos::parallel_for(Ai.dimension_0(),set_inner<inner_mesh_type>(Ai,0));
  Kokkos::parallel_for(Ai.dimension_0(),update<mesh_type>(A,0.1));
  


  printf("Done\n");

  Kokkos::finalize();
}

