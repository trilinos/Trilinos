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

typedef Kokkos::DualView<double*> view_type;
typedef Kokkos::DualView<int**> idx_type;


template<class Device>
struct localsum {
  // Define the execution space for the functor (overrides the DefaultExecutionSpace)
  typedef Device device_type;

  // Get the view types on the particular device the functor is instantiated for
  Kokkos::View<idx_type::const_data_type, idx_type::array_layout, Device> idx;
  Kokkos::View<view_type::array_type, view_type::array_layout, Device> dest;
  Kokkos::View<view_type::const_data_type, view_type::array_layout, Device, Kokkos::MemoryRandomAccess > src;

  localsum(idx_type dv_idx, view_type dv_dest,
      view_type dv_src) {
    // Extract the view on the correct Device from the DualView
    idx = dv_idx.view<Device>();
    dest = dv_dest.template view<Device>();
    src = dv_src.template view<Device>();

    // Synchronize the DualView on the correct Device
    dv_idx.sync<Device>();
    dv_dest.template sync<Device>();
    dv_src.template sync<Device>();

    // Mark dest as modified on Device
    dv_dest.template modify<Device>();
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    double tmp = 0.0;
    for(int j = 0; j < idx.dimension_1(); j++) {
      const double val = src(idx(i,j));
      tmp += val*val + 0.5*(idx.dimension_0()*val -idx.dimension_1()*val);
    }
    dest(i) += tmp;
  }
};

int main(int narg, char* arg[]) {
  Kokkos::initialize(narg,arg);

  int size = 1000000;

  // Create DualViews. This will allocate on both the device and its
  // host_mirror_device.
  idx_type idx("Idx",size,64);
  view_type dest("Dest",size);
  view_type src("Src",size);

  srand(134231);

  // Get a reference to the host view of idx directly (equivalent to
  // idx.view<idx_type::host_mirror_space>() )
  idx_type::t_host h_idx = idx.h_view;
  for (int i = 0; i < size; ++i) {
    for (view_type::size_type j=0; j < h_idx.dimension_1 (); ++j) {
      h_idx(i,j) = (size + i + (rand () % 500 - 250)) % size;
    }
  }

  // Mark idx as modified on the host_mirror_space so that a
  // sync to the device will actually move data.  The sync happens in
  // the functor's constructor.
  idx.modify<idx_type::host_mirror_space>();

  // Run on the device.  This will cause a sync of idx to the device,
  // since it was marked as modified on the host.
  Kokkos::Impl::Timer timer;
  Kokkos::parallel_for(size,localsum<view_type::device_type>(idx,dest,src));
  Kokkos::fence();
  double sec1_dev = timer.seconds();

  timer.reset();
  Kokkos::parallel_for(size,localsum<view_type::device_type>(idx,dest,src));
  Kokkos::fence();
  double sec2_dev = timer.seconds();

  // Run on the host's default execution space (could be the same as device).
  // This will cause a sync back to the host of dest.  Note that if the Device is CUDA,
  // the data layout will not be optimal on host, so performance is
  // lower than what it would be for a pure host compilation.
  timer.reset();
  Kokkos::parallel_for(size,localsum<Kokkos::HostSpace::execution_space>(idx,dest,src));
  Kokkos::fence();
  double sec1_host = timer.seconds();

  timer.reset();
  Kokkos::parallel_for(size,localsum<Kokkos::HostSpace::execution_space>(idx,dest,src));
  Kokkos::fence();
  double sec2_host = timer.seconds();

  printf("Device Time with Sync: %lf without Sync: %lf \n",sec1_dev,sec2_dev);
  printf("Host   Time with Sync: %lf without Sync: %lf \n",sec1_host,sec2_host);

  Kokkos::finalize();
}

