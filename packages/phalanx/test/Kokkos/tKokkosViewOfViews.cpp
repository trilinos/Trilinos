// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Kokkos_Core.hpp"
#include "Teuchos_UnitTestHarness.hpp"

// Force this test to always run without UVM.
// For Cuda builds, ignore the default memory space in the Cuda
// execution space since it can be set to UVM at
// configure. For non-Cuda builds, just use the default memory space
// in the execution space.
namespace Kokkos {
  class Cuda;
  class CudaSpace;
}
using exec_t = Kokkos::DefaultExecutionSpace;
using mem_t = std::conditional<std::is_same<exec_t,Kokkos::Cuda>::value,
                               Kokkos::CudaSpace,
                               Kokkos::DefaultExecutionSpace::memory_space>::type;

// *************************************
// This test demonstrates how to create a view of views from separate
// already allocated views for double and FAD types.
// *************************************
TEUCHOS_UNIT_TEST(ViewOfViews,from_separate_views) {

  const int num_cells = 10;
  const int num_pts = 8;
  const int num_equations = 32;

  // Requirement 1: The inner view must be unmanaged on device to
  // prevent double deletion! To initialize correctly, we need to deep
  // copy from host with the inner views propeties matching exactly on
  // host and device.
  using InnerViewUnmanaged = Kokkos::View<double***,mem_t,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using InnerViewManaged = Kokkos::View<double***,mem_t>;
  using OuterViewDeviceUnmanaged = Kokkos::View<InnerViewUnmanaged*,mem_t>;
  using OuterViewHostMirrorUnmanaged = Kokkos::View<InnerViewUnmanaged*,mem_t>::HostMirror; // all inner view args must match for deep_copy to device
  using OuterViewHostMirrorManaged = Kokkos::View<InnerViewManaged*,mem_t>::HostMirror; // used to store views so they don't go out of scope

  // Requirement 2: The host view must exist for the life of the
  // device view! Need to make sure the views are managed on
  // host. However the deep_copy to device needs the matching memory
  // management type on the inner view. So we need an unmanaged
  // version also.
  OuterViewHostMirrorManaged v_host_managed("outer view host mirror managed",3);
  OuterViewHostMirrorUnmanaged v_host_unmanaged("outer view host mirror unmanaged",3);
  {
    // Let the original views go out of scope so that we proves the
    // managed view of views is keeping device allocations alive.
    {
      Kokkos::View<double***,mem_t> a("a",num_cells,num_pts,num_equations);
      Kokkos::View<double***,mem_t> b("b",num_cells,num_pts,num_equations);
      Kokkos::View<double***,mem_t> c("c",num_cells,num_pts,num_equations);

      Kokkos::deep_copy(a,2.0);
      Kokkos::deep_copy(b,3.0);

      // To initialize the inner views on device correctly.
      v_host_unmanaged(0) = a;
      v_host_unmanaged(1) = b;
      v_host_unmanaged(2) = c;

      // To make sure the inner views don't go out of scope and delete.
      v_host_managed(0) = a;
      v_host_managed(1) = b;
      v_host_managed(2) = c;
    }

    // Requirement 3: Need to deep_copy the host view to device to
    // initialize the inner views correctly on device.
    OuterViewDeviceUnmanaged v_dev("outer view device unmanaged",3);
    Kokkos::deep_copy(v_dev,v_host_unmanaged);

    auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0},{num_cells,num_pts,num_equations});
    Kokkos::parallel_for("view of view test",policy,KOKKOS_LAMBDA (const int cell,const int pt, const int eq) {
      v_dev(2)(cell,pt,eq) = v_dev(0)(cell,pt,eq) + v_dev(1)(cell,pt,eq);
    });

    auto c_device = v_host_managed(2);
    auto c_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),c_device);

    const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
    for (int cell=0; cell < num_cells; ++cell)
      for (int pt=0; pt < num_pts; ++pt)
        for (int eq=0; eq < num_equations; ++eq) {
          TEST_FLOATING_EQUALITY(c_host(cell,pt,eq),5.0,tol);
        }
  }

  // Need to manually clean out the inner view memory. The Kokkos view
  // dtor calls a parallel_for in the destructor for the view. A view
  // of views thus calls nested parallel_fors that are not
  // allowed. This results in a hang at the parallel_for mutex.
  {
    for (size_t i=0; i < v_host_managed.size(); ++i)
      v_host_managed(i) = {};
  }
}

// Note: The outer view will call the default ctor on device for first
// touch during allocation. If a default ctor is not avaialble, we can
// create the outer view uninitialized and manually call new/delete:

// Kokkos::View<InnerView*> v_host(Kokkos::view_alloc("c",Kokkos::WithoutInitializing),100,100);
// for (size_t i=0; i < v_host.extent(0); ++i)
//   new (&v_host(i)) InnerView(a);
// for (size_t i=0; i < v_host.extent(0); ++i)
//   v_host(i).~InnerView();

// *************************************
// This test checks if we can use a Kokkos::Array with view of views
// *************************************
TEUCHOS_UNIT_TEST(ViewOfViews,ArrayOfViews) {

  const int num_cells = 10;
  const int num_pts = 8;
  const int num_equations = 32;

  using InnerView = Kokkos::View<double***,mem_t>;

  // For making sure the views are not destroyed.
  // std::vector<InnerView> vec(3);

  Kokkos::Array<InnerView,3> a_of_v;
  {
    InnerView a("a",num_cells,num_pts,num_equations);
    InnerView b("b",num_cells,num_pts,num_equations);
    InnerView c("c",num_cells,num_pts,num_equations);

    a_of_v[0] = a;
    a_of_v[1] = b;
    a_of_v[2] = c;

    Kokkos::deep_copy(a,1.0);
    Kokkos::deep_copy(b,2.0);
    Kokkos::deep_copy(c,3.0);
  }
  
  auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0},{num_cells,num_pts,num_equations});
  Kokkos::parallel_for("view of view test",policy,KOKKOS_LAMBDA (const int cell,const int pt, const int eq) {
    a_of_v[2](cell,pt,eq) += a_of_v[0](cell,pt,eq) + a_of_v[1](cell,pt,eq);
  });

  auto c_host = Kokkos::create_mirror_view(a_of_v[2]);
  Kokkos::deep_copy(c_host,a_of_v[2]);
  const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
  for (int cell=0; cell < num_cells; ++cell)
    for (int pt=0; pt < num_pts; ++pt)
      for (int eq=0; eq < num_equations; ++eq) {
        TEST_FLOATING_EQUALITY(c_host(cell,pt,eq),6.0,tol);
      }

}
  
