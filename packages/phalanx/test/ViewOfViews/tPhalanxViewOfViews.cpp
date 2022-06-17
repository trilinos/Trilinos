#include "Kokkos_Core.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Phalanx_KokkosViewOfViews.hpp"
#include <vector>

// ********************************
// This test demonstrates how to create a view of views for double and FAD types.
// Original implementation
// ********************************
using exec_t = Kokkos::DefaultExecutionSpace;
using mem_t = Kokkos::DefaultExecutionSpace::memory_space;

TEUCHOS_UNIT_TEST(PhalanxViewOfViews,OldImpl) {

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

// ********************************
// New implementation (automatically adds the unmanaged memory trait to inner view)
// ********************************
TEUCHOS_UNIT_TEST(PhalanxViewOfViews,NewImpl) {

  const int num_cells = 10;
  const int num_pts = 8;
  const int num_equations = 32;

  using InnerView = Kokkos::View<double***,mem_t>;
  constexpr int OuterViewRank = 2;
  PHX::ViewOfViews2<OuterViewRank,InnerView,mem_t> v_of_v("outer host",2,2);

  {

    // Let originals go out of scope to check correct memory management
    {
      Kokkos::View<double***,mem_t> a("a",num_cells,num_pts,num_equations);
      Kokkos::View<double***,mem_t> b("b",num_cells,num_pts,num_equations);
      Kokkos::View<double***,mem_t> c("c",num_cells,num_pts,num_equations);
      Kokkos::View<double***,mem_t> d("d",num_cells,num_pts,num_equations);

      Kokkos::deep_copy(a,2.0);
      Kokkos::deep_copy(b,3.0);
      Kokkos::deep_copy(c,4.0);

      v_of_v.setView(a,0,0);
      v_of_v.setView(b,0,1);
      v_of_v.setView(c,1,0);
      v_of_v.setView(d,1,1);
    }

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

  auto d = v_of_v.getViewHost()(1,1);
  auto d_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),d);

  const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
  for (int cell=0; cell < num_cells; ++cell)
    for (int pt=0; pt < num_pts; ++pt)
      for (int eq=0; eq < num_equations; ++eq) {
        TEST_FLOATING_EQUALITY(d_host(cell,pt,eq),9.0,tol);
      }

}

TEUCHOS_UNIT_TEST(PhalanxViewOfViews,ViewOfView3_EmptyCtor) {

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
    using InnerView = Kokkos::View<double***,mem_t>;
    constexpr int OuterViewRank = 2;
    PHX::ViewOfViews3<OuterViewRank,InnerView,mem_t> v_of_v;

    TEST_ASSERT(!v_of_v.isInitialized());
    v_of_v.initialize("outer host",2,2);
    TEST_ASSERT(v_of_v.isInitialized());

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

    v_of_v.disableSafetyCheck();
    v_of_v.enableSafetyCheck();
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

TEUCHOS_UNIT_TEST(PhalanxViewOfViews,ViewOfView3_TwoArgCtor) {

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
    using InnerView = Kokkos::View<double***,mem_t>;
    constexpr int OuterViewRank = 2;
    PHX::ViewOfViews3<OuterViewRank,InnerView,mem_t> v_of_v("outer host",2,2);

    TEST_ASSERT(v_of_v.isInitialized());

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

// Make sure that an uninitialized ViewOviews3 can be default
// constructed and destoryed. Happens in application unit tests.
struct MeshEvaluationTestStruct {
    using InnerView = Kokkos::View<double***,mem_t>;
    PHX::ViewOfViews3<2,InnerView,mem_t> v_of_v_;
};

TEUCHOS_UNIT_TEST(PhalanxViewOfViews,ViewOfView3_DefaultCtorDtor) {
  auto mesh_eval = std::make_shared<MeshEvaluationTestStruct>();
  TEST_ASSERT(!mesh_eval->v_of_v_.isInitialized());
  mesh_eval = nullptr;
}

// ********************************
// Demonstrates an alternative path for ViewOfViews that uses a user
// defined wrapper and the assignment operator on device to disable
// the reference tracking.
// ********************************

// Force this test to always run without UVM.
// For Cuda builds, ignore the default memory space in the Cuda
// execution space since it can be set to UVM in Trilinos
// configure. For non-Cuda builds, just use the default memory space
// in the execution space.
namespace Kokkos {
  class Cuda;
  class CudaSpace;
}
using DeviceExecutionSpace = Kokkos::DefaultExecutionSpace;
using DeviceMemorySpace = std::conditional<std::is_same<DeviceExecutionSpace,Kokkos::Cuda>::value,
                                           Kokkos::CudaSpace,
                                           Kokkos::DefaultExecutionSpace::memory_space>::type;
using view = Kokkos::View<double*,DeviceMemorySpace>;
using view_host = view::HostMirror;

class Wrapper {
public:
  view a_;
  view b_;
  KOKKOS_INLINE_FUNCTION
  Wrapper() : a_(nullptr,0),b_(nullptr,0) {}

  Wrapper(const view_host& a, const view_host& b)
    : a_(a.label(),a.layout()),b_(b.label(),b.layout())
  {
    Kokkos::deep_copy(a_,a);
    Kokkos::deep_copy(b_,b);
  }

  KOKKOS_DEFAULTED_FUNCTION
  Wrapper& operator=(const Wrapper& src) = default;

  KOKKOS_INLINE_FUNCTION
  double multiply(int i) const {return a_(i) * b_(i);}
};

TEUCHOS_UNIT_TEST(PhalanxViewOfViews,WrapperExample) {

  constexpr int num_objects = 3;
  constexpr int view_size = 20;

  // This object must exist for lifetime if v_of_uo to keep inner
  // views in scope.
  std::vector<Wrapper> uo(num_objects);
  {
    view_host a("a",view_size);
    view_host b("b",view_size);
    view_host c("c",view_size);
    view_host d("d",view_size);
    view_host e("e",view_size);
    view_host f("f",view_size);

    Kokkos::deep_copy(a,0.0);
    Kokkos::deep_copy(b,1.0);
    Kokkos::deep_copy(c,2.0);
    Kokkos::deep_copy(d,3.0);
    Kokkos::deep_copy(e,4.0);
    Kokkos::deep_copy(f,5.0);

    uo[0] = Wrapper(a,b);
    uo[1] = Wrapper(c,d);
    uo[2] = Wrapper(e,f);
  }

  Kokkos::View<Wrapper*,DeviceMemorySpace> v_of_uo("v_of_uo",num_objects);

  for (int i=0; i < 3; ++i) {
    Wrapper tmp = uo[i];
    Kokkos::parallel_for("initialize view_of_uo",1,KOKKOS_LAMBDA(const int ) {
      // reference counting is disabled on device
      v_of_uo(i) = tmp;
    });
    Kokkos::fence();
  }

  Kokkos::View<double**,DeviceMemorySpace> results("results",num_objects,view_size);
  Kokkos::MDRangePolicy<exec_t,Kokkos::Rank<2>> policy({0,0},{num_objects,view_size});
  Kokkos::parallel_for("v_of_uo",policy,KOKKOS_LAMBDA(const int i,const int j) {
    results(i,j) = v_of_uo(i).multiply(j) + static_cast<double>(j);
  });
  Kokkos::fence();

  auto results_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),results);

  constexpr auto tol = std::numeric_limits<double>::epsilon() * 100.0;
  for (int i=0; i < num_objects; ++i) {
    for (int j=0; j < view_size; ++j) {

      double gold_value = static_cast<double>(j);
      if (i == 1)
        gold_value += 6.0;
      else if (i == 2)
        gold_value += 20.0;

      TEST_FLOATING_EQUALITY(results_host(i,j),gold_value,tol);
    }
  }
}

TEUCHOS_UNIT_TEST(PhalanxViewOfViews,CreateHostHost) {

  const int num_cells = 5;

  using InnerView = Kokkos::View<double*,mem_t>;
  InnerView a("a",num_cells);
  InnerView b("b",num_cells);
  InnerView c("c",num_cells);
  InnerView d("d",num_cells);

  Kokkos::deep_copy(a,2.0);
  Kokkos::deep_copy(b,3.0);
  Kokkos::deep_copy(c,4.0);

  // Rank 1 outer view
  {
    PHX::ViewOfViews3<1,InnerView,mem_t> pvov("vov1",4);
    pvov.addView(a,0);
    pvov.addView(b,1);
    pvov.addView(c,2);
    pvov.addView(d,3);
    pvov.syncHostToDevice();

    auto vov = pvov.getViewDevice();

    Kokkos::parallel_for("vov1",num_cells,KOKKOS_LAMBDA(const int cell) {
      vov(3)(cell) = vov(0)(cell) * vov(1)(cell) + vov(2)(cell);
    });

    auto vov_host = PHX::createHostHostViewOfViews(vov);

    const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
    for (int cell=0; cell < num_cells; ++cell) {
      TEST_FLOATING_EQUALITY(vov_host(0)(cell),2.0,tol);
      TEST_FLOATING_EQUALITY(vov_host(1)(cell),3.0,tol);
      TEST_FLOATING_EQUALITY(vov_host(2)(cell),4.0,tol);
      TEST_FLOATING_EQUALITY(vov_host(3)(cell),10.0,tol);
    }
  }

  // Rank 2 outer view
  {
    PHX::ViewOfViews3<2,InnerView,mem_t> pvov("vov1",2,2);
    pvov.addView(a,0,0);
    pvov.addView(b,0,1);
    pvov.addView(c,1,0);
    pvov.addView(d,1,1);
    pvov.syncHostToDevice();

    auto vov = pvov.getViewDevice();

    Kokkos::parallel_for("vov1",num_cells,KOKKOS_LAMBDA(const int cell) {
      vov(1,1)(cell) = vov(0,0)(cell) * vov(0,1)(cell) + vov(1,0)(cell) + 1.0;
    });

    auto vov_host = PHX::createHostHostViewOfViews(vov);

    const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
    for (int cell=0; cell < num_cells; ++cell) {
      TEST_FLOATING_EQUALITY(vov_host(0,0)(cell),2.0,tol);
      TEST_FLOATING_EQUALITY(vov_host(0,1)(cell),3.0,tol);
      TEST_FLOATING_EQUALITY(vov_host(1,0)(cell),4.0,tol);
      TEST_FLOATING_EQUALITY(vov_host(1,1)(cell),11.0,tol);
    }
  }

  // Rank 3 outer view
  {
    PHX::ViewOfViews3<3,InnerView,mem_t> pvov("vov1",3,3,3);
    pvov.addView(a,0,0,0);
    pvov.addView(b,1,1,1);
    pvov.addView(c,2,2,2);
    pvov.addView(d,0,1,2);
    pvov.syncHostToDevice();

    auto vov = pvov.getViewDevice();

    Kokkos::parallel_for("vov1",num_cells,KOKKOS_LAMBDA(const int cell) {
        vov(0,1,2)(cell) = vov(0,0,0)(cell) * vov(1,1,1)(cell) + vov(2,2,2)(cell) + 2.0;
    });

    auto vov_host = PHX::createHostHostViewOfViews(vov);

    const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
    for (int cell=0; cell < num_cells; ++cell) {
      TEST_FLOATING_EQUALITY(vov_host(0,0,0)(cell),2.0,tol);
      TEST_FLOATING_EQUALITY(vov_host(1,1,1)(cell),3.0,tol);
      TEST_FLOATING_EQUALITY(vov_host(2,2,2)(cell),4.0,tol);
      TEST_FLOATING_EQUALITY(vov_host(0,1,2)(cell),12.0,tol);
    }
  }

}
