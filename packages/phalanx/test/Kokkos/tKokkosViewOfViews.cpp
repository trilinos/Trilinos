#include "Kokkos_Core.hpp"
#include "Teuchos_UnitTestHarness.hpp"

using exec_t = Kokkos::DefaultExecutionSpace;
using mem_t = Kokkos::DefaultExecutionSpace::memory_space;

// *************************************
// This test demonstrates how to create a view of views from separate
// already allocated views for double and FAD types.
// *************************************
TEUCHOS_UNIT_TEST(ViewOfViews,from_separate_views) {

  const int num_cells = 10;
  const int num_pts = 8;
  const int num_equations = 32;

  Kokkos::View<double***,mem_t> a("a",num_cells,num_pts,num_equations);
  Kokkos::View<double***,mem_t> b("b",num_cells,num_pts,num_equations);
  Kokkos::View<double***,mem_t> c("c",num_cells,num_pts,num_equations);

  Kokkos::deep_copy(a,2.0);
  Kokkos::deep_copy(b,3.0);

  // Requirement 1: The inner view must be unmanaged to prevent double deletion!
  using InnerView = Kokkos::View<double***,mem_t,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using OuterView = Kokkos::View<InnerView*,mem_t>;

  // Requirement 2: The host view must exist for the life of the device view!
  OuterView::HostMirror v_host("outer host",3);
  {
    v_host(0) = a;
    v_host(1) = b;
    v_host(2) = c;

    // Requirement 3: Need to deep_copy the host view to device to
    // initialize the inner views correctly.
    OuterView v_dev("outer device",3);
    Kokkos::deep_copy(v_dev,v_host);

    auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0},{num_cells,num_pts,num_equations});
    Kokkos::parallel_for("view of view test",policy,KOKKOS_LAMBDA (const int cell,const int pt, const int eq) {
      v_dev(2)(cell,pt,eq) = v_dev(0)(cell,pt,eq) + v_dev(1)(cell,pt,eq);
    });

    auto c_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),c);

    const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
    for (int cell=0; cell < num_cells; ++cell)
      for (int pt=0; pt < num_pts; ++pt)
        for (int eq=0; eq < num_equations; ++eq) {
          TEST_FLOATING_EQUALITY(c_host(cell,pt,eq),5.0,tol);
        }
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
