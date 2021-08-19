#include "Kokkos_Core.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Phalanx_KokkosViewOfViews.hpp"

// This test demonstrates how to create a view of views for double and FAD types.
using exec_t = Kokkos::DefaultExecutionSpace;
using mem_t = Kokkos::DefaultExecutionSpace::memory_space;

TEUCHOS_UNIT_TEST(PhalanxViewOfViews,double) {

  const int num_cells = 10;
  const int num_pts = 8;
  const int num_equations = 32;

  Kokkos::View<double***,mem_t> a("a",num_cells,num_pts,num_equations);
  Kokkos::View<double***,mem_t> b("b",num_cells,num_pts,num_equations);
  Kokkos::View<double***,mem_t> c("c",num_cells,num_pts,num_equations);
  Kokkos::View<double***,mem_t> d("d",num_cells,num_pts,num_equations);

  Kokkos::deep_copy(a,2.0);
  Kokkos::deep_copy(b,3.0);
  Kokkos::deep_copy(c,4.0);

  {
    using InnerView = Kokkos::View<double***,mem_t,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    constexpr int OuterViewRank = 2;
    PHX::ViewOfViews<OuterViewRank,InnerView,mem_t> v_of_v("outer host",2,2);
    
    v_of_v.addView(a,0,0);
    v_of_v.addView(b,0,1);
    v_of_v.addView(c,1,0);
    v_of_v.addView(d,1,1);
    
    v_of_v.syncHostToDevice();
   
    { 
      auto v_dev = v_of_v.getViewDevice();
      auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0},{num_cells,num_pts,num_equations});
      Kokkos::parallel_for("view of view test",policy,KOKKOS_LAMBDA (const int cell,const int pt, const int eq) {
        v_dev(1,1)(cell,pt,eq) = v_dev(0,0)(cell,pt,eq) + v_dev(0,1)(cell,pt,eq) + v_dev(1,0)(cell,pt,eq);
      });
    }

    // Uncomment the line below to prove the ViewOfViews prevents
    // device views from outliving host view. This line will cause a
    // Kokkos::abort() and error message since v_dev above is still in
    // scope when the ViewOfViews is destoryed.
    // v_of_v = PHX::ViewOfViews<OuterViewRank,InnerView,mem_t>("outer host",2,2);
  }

  auto d_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),d);
  
  const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
  for (int cell=0; cell < num_cells; ++cell)
    for (int pt=0; pt < num_pts; ++pt)
      for (int eq=0; eq < num_equations; ++eq) {
        TEST_FLOATING_EQUALITY(d_host(cell,pt,eq),9.0,tol);
      }
}
