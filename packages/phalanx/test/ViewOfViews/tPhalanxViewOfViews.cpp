// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Sacado.hpp"
#include "Kokkos_View_Fad.hpp"
#include "Kokkos_Core.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_KokkosViewOfViews.hpp"
#include "Phalanx_Kokkos_Tools_CheckStreams.hpp"
#include <vector>

// ********************************
// This test demonstrates how to create a view of views for double and FAD types.
// Original implementation
// ********************************
using exec_t = Kokkos::DefaultExecutionSpace;
using mem_t = Kokkos::DefaultExecutionSpace::memory_space;

TEUCHOS_UNIT_TEST(PhalanxViewOfViews,ViewOfView_DefaultStreamInitialize) {

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
    PHX::ViewOfViews<OuterViewRank,InnerView,mem_t> v_of_v;

    TEST_ASSERT(!v_of_v.isInitialized());
    v_of_v.initialize("outer host",2,2);
    TEST_ASSERT(v_of_v.isInitialized());

    v_of_v.addView(a,0,0);
    v_of_v.addView(b,0,1);
    v_of_v.addView(c,1,0);
    v_of_v.addView(d,1,1);

    TEST_ASSERT(!v_of_v.deviceViewIsSynced());
    v_of_v.syncHostToDevice();
    TEST_ASSERT(v_of_v.deviceViewIsSynced());

    {
      auto v_dev = v_of_v.getViewDevice();
      auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0},{num_cells,num_pts,num_equations});
      Kokkos::parallel_for("view of view test",policy,KOKKOS_LAMBDA (const int cell,const int pt, const int eq) {
        v_dev(1,1)(cell,pt,eq) = v_dev(0,0)(cell,pt,eq) + v_dev(0,1)(cell,pt,eq) + v_dev(1,0)(cell,pt,eq);
      });
    }

    // Test the query objectes used for VT serialization
    TEST_ASSERT(v_of_v.safetyCheck());
    v_of_v.disableSafetyCheck();
    TEST_ASSERT(!v_of_v.safetyCheck());
    v_of_v.enableSafetyCheck();
    TEST_ASSERT(v_of_v.safetyCheck());

    // Test the const versions of accessors
    {
      const PHX::ViewOfViews<OuterViewRank,InnerView,mem_t>& const_v_of_v = v_of_v;
      const_v_of_v.getViewHost();
      const_v_of_v.getViewDevice();
    }
  }

  auto d_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),d);

  const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
  for (int cell=0; cell < num_cells; ++cell)
    for (int pt=0; pt < num_pts; ++pt)
      for (int eq=0; eq < num_equations; ++eq) {
        TEST_FLOATING_EQUALITY(d_host(cell,pt,eq),9.0,tol);
      }
}

TEUCHOS_UNIT_TEST(PhalanxViewOfViews,ViewOfView_DefaultStreamCtor) {

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
    PHX::ViewOfViews<OuterViewRank,InnerView,mem_t> v_of_v("outer host",2,2);

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
  }

  auto d_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),d);

  const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
  for (int cell=0; cell < num_cells; ++cell)
    for (int pt=0; pt < num_pts; ++pt)
      for (int eq=0; eq < num_equations; ++eq) {
        TEST_FLOATING_EQUALITY(d_host(cell,pt,eq),9.0,tol);
      }
}

TEUCHOS_UNIT_TEST(PhalanxViewOfViews,ViewOfView_UserStreamCtor) {

  const int num_cells = 10;
  const int num_pts = 8;
  const int num_equations = 32;

  std::vector<PHX::Device> streams;
  if (PHX::Device().concurrency() >= 4) {
    std::cout << "Using partition_space, concurrency=" << PHX::Device().concurrency() << std::endl;
    streams = Kokkos::Experimental::partition_space(PHX::Device(),1,1,1,1);
  }
  else {
    std::cout << "NOT using partition_space, concurrency=" << PHX::Device().concurrency() << std::endl;
    for (int i=0; i < 4; ++i)
      streams.push_back(PHX::Device());
  }

  PHX::set_enforce_no_default_stream_use();

  Kokkos::View<double***,mem_t> a(Kokkos::view_alloc(streams[0],"a"),num_cells,num_pts,num_equations);
  Kokkos::View<double***,mem_t> b(Kokkos::view_alloc(streams[1],"b"),num_cells,num_pts,num_equations);
  Kokkos::View<double***,mem_t> c(Kokkos::view_alloc(streams[2],"c"),num_cells,num_pts,num_equations);
  Kokkos::View<double***,mem_t> d(Kokkos::view_alloc(streams[3],"d"),num_cells,num_pts,num_equations);

  Kokkos::deep_copy(streams[0],a,2.0);
  Kokkos::deep_copy(streams[1],b,3.0);
  Kokkos::deep_copy(streams[2],c,4.0);

  streams[0].fence();
  streams[1].fence();
  streams[2].fence();
  streams[3].fence();

  {
    using InnerView = Kokkos::View<double***,mem_t>;
    constexpr int OuterViewRank = 2;
    PHX::ViewOfViews<OuterViewRank,InnerView,mem_t> v_of_v(streams[3],"outer host",2,2);

    TEST_ASSERT(v_of_v.isInitialized());

    // For UVM=ON builds, need to fence here. The outer views are
    // being accessed before the initialization is completed on
    // device. The failure only shows up if you overload the cuda card
    // with multiple tests to slow down the initialization.
    streams[3].fence(); // PHX_UVM_ON_FENCE

    v_of_v.addView(a,0,0);
    v_of_v.addView(b,0,1);
    v_of_v.addView(c,1,0);
    v_of_v.addView(d,1,1);

    v_of_v.syncHostToDevice(streams[3]);

    {
      auto v_dev = v_of_v.getViewDevice();
      auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<3>>(streams[3],{0,0,0},{num_cells,num_pts,num_equations});
      Kokkos::parallel_for("view of view test",policy,KOKKOS_LAMBDA (const int cell,const int pt, const int eq) {
        v_dev(1,1)(cell,pt,eq) = v_dev(0,0)(cell,pt,eq) + v_dev(0,1)(cell,pt,eq) + v_dev(1,0)(cell,pt,eq);
      });
    }

  }

  streams[3].fence();

  PHX::unset_enforce_no_default_stream_use();

  auto d_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),d);

  const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
  for (int cell=0; cell < num_cells; ++cell)
    for (int pt=0; pt < num_pts; ++pt)
      for (int eq=0; eq < num_equations; ++eq) {
        TEST_FLOATING_EQUALITY(d_host(cell,pt,eq),9.0,tol);
      }
}

TEUCHOS_UNIT_TEST(PhalanxViewOfViews,ViewOfView_UserStreamInitialize) {

  const int num_cells = 10;
  const int num_pts = 8;
  const int num_equations = 32;

  std::vector<PHX::Device> streams;
  if (PHX::Device().concurrency() >= 4) {
    std::cout << "Using partition_space, concurrency=" << PHX::Device().concurrency() << std::endl;
    streams = Kokkos::Experimental::partition_space(PHX::Device(),1,1,1,1);
  }
  else {
    std::cout << "NOT using partition_space, concurrency=" << PHX::Device().concurrency() << std::endl;
    for (int i=0; i < 4; ++i)
      streams.push_back(PHX::Device());
  }

  PHX::set_enforce_no_default_stream_use();

  Kokkos::View<double***,mem_t> a(Kokkos::view_alloc(streams[0],"a"),num_cells,num_pts,num_equations);
  Kokkos::View<double***,mem_t> b(Kokkos::view_alloc(streams[1],"b"),num_cells,num_pts,num_equations);
  Kokkos::View<double***,mem_t> c(Kokkos::view_alloc(streams[2],"c"),num_cells,num_pts,num_equations);
  Kokkos::View<double***,mem_t> d(Kokkos::view_alloc(streams[3],"d"),num_cells,num_pts,num_equations);

  Kokkos::deep_copy(streams[0],a,2.0);
  Kokkos::deep_copy(streams[1],b,3.0);
  Kokkos::deep_copy(streams[2],c,4.0);

  streams[0].fence();
  streams[1].fence();
  streams[2].fence();
  streams[3].fence();

  {
    using InnerView = Kokkos::View<double***,mem_t>;
    constexpr int OuterViewRank = 2;
    PHX::ViewOfViews<OuterViewRank,InnerView,mem_t> v_of_v;

    TEST_ASSERT(!v_of_v.isInitialized());
    v_of_v.initialize(streams[3],"outer host",2,2);
    TEST_ASSERT(v_of_v.isInitialized());

    // For UVM=ON builds, need to fence here. The outer views are
    // being accessed before the initialization is completed on
    // device. The failure only shows up if you overload the cuda card
    // with multiple tests to slow down the initialization.
    streams[3].fence(); // PHX_UVM_ON_FENCE

    v_of_v.addView(a,0,0);
    v_of_v.addView(b,0,1);
    v_of_v.addView(c,1,0);
    v_of_v.addView(d,1,1);

    v_of_v.syncHostToDevice(streams[3]);

    {
      auto v_dev = v_of_v.getViewDevice();
      auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<3>>(streams[3],{0,0,0},{num_cells,num_pts,num_equations});
      Kokkos::parallel_for("view of view test",policy,KOKKOS_LAMBDA (const int cell,const int pt, const int eq) {
        v_dev(1,1)(cell,pt,eq) = v_dev(0,0)(cell,pt,eq) + v_dev(0,1)(cell,pt,eq) + v_dev(1,0)(cell,pt,eq);
      });
    }
  }

  streams[3].fence();

  PHX::unset_enforce_no_default_stream_use();

  auto d_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),d);

  const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
  for (int cell=0; cell < num_cells; ++cell)
    for (int pt=0; pt < num_pts; ++pt)
      for (int eq=0; eq < num_equations; ++eq) {
        TEST_FLOATING_EQUALITY(d_host(cell,pt,eq),9.0,tol);
      }
}

/* Temporarily disable:
   https://github.com/kokkos/kokkos-tools/issues/224
   https://github.com/trilinos/Trilinos/pull/11391

TEUCHOS_UNIT_TEST(PhalanxViewOfViews,KokkosToolsDefaultStreamCheck) {
  PHX::set_enforce_no_default_stream_use();
  // Checks are only active for CUDA and HIP backends
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  TEST_THROW(PHX::Device().fence(),std::runtime_error);
#endif
  PHX::unset_enforce_no_default_stream_use();
}

*/

// Make sure that an uninitialized ViewOviews3 can be default
// constructed and destroyed. Happens in application unit tests.
struct MeshEvaluationTestStruct {
    using InnerView = Kokkos::View<double***,mem_t>;
    PHX::ViewOfViews<2,InnerView,mem_t> v_of_v_;
};

TEUCHOS_UNIT_TEST(PhalanxViewOfViews,ViewOfView_DefaultCtorDtor) {
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
    PHX::ViewOfViews<1,InnerView,mem_t> pvov("vov1",4);
    pvov.addView(a,0);
    pvov.addView(b,1);
    pvov.addView(c,2);
    pvov.addView(d,3);
    pvov.syncHostToDevice();

    auto vov = pvov.getViewDevice();

    Kokkos::parallel_for("vov1",num_cells,KOKKOS_LAMBDA(const int cell) {
      vov(3)(cell) = vov(0)(cell) * vov(1)(cell) + vov(2)(cell);
    });

    auto pvov_host = PHX::createHostHostViewOfViews(pvov);
    {
      auto vov_host = pvov_host.getViewHost();

      const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
      for (int cell=0; cell < num_cells; ++cell) {
        TEST_FLOATING_EQUALITY(vov_host(0)(cell),2.0,tol);
        TEST_FLOATING_EQUALITY(vov_host(1)(cell),3.0,tol);
        TEST_FLOATING_EQUALITY(vov_host(2)(cell),4.0,tol);
        TEST_FLOATING_EQUALITY(vov_host(3)(cell),10.0,tol);
      }
    }
  }

  // Rank 2 outer view
  {
    PHX::ViewOfViews<2,InnerView,mem_t> pvov("vov1",2,2);
    pvov.addView(a,0,0);
    pvov.addView(b,0,1);
    pvov.addView(c,1,0);
    pvov.addView(d,1,1);
    pvov.syncHostToDevice();

    auto vov = pvov.getViewDevice();

    Kokkos::parallel_for("vov1",num_cells,KOKKOS_LAMBDA(const int cell) {
      vov(1,1)(cell) = vov(0,0)(cell) * vov(0,1)(cell) + vov(1,0)(cell) + 1.0;
    });

    auto pvov_host = PHX::createHostHostViewOfViews(pvov);
    {
      auto vov_host = pvov_host.getViewHost();

      const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
      for (int cell=0; cell < num_cells; ++cell) {
        TEST_FLOATING_EQUALITY(vov_host(0,0)(cell),2.0,tol);
        TEST_FLOATING_EQUALITY(vov_host(0,1)(cell),3.0,tol);
        TEST_FLOATING_EQUALITY(vov_host(1,0)(cell),4.0,tol);
        TEST_FLOATING_EQUALITY(vov_host(1,1)(cell),11.0,tol);
      }
    }
  }

  // Rank 3 outer view
  {
    PHX::ViewOfViews<3,InnerView,mem_t> pvov("vov1",3,3,3);
    pvov.addView(a,0,0,0);
    pvov.addView(b,1,1,1);
    pvov.addView(c,2,2,2);
    pvov.addView(d,0,1,2);
    pvov.syncHostToDevice();

    auto vov = pvov.getViewDevice();

    Kokkos::parallel_for("vov1",num_cells,KOKKOS_LAMBDA(const int cell) {
        vov(0,1,2)(cell) = vov(0,0,0)(cell) * vov(1,1,1)(cell) + vov(2,2,2)(cell) + 2.0;
    });

    auto pvov_host = PHX::createHostHostViewOfViews(pvov);
    {
      auto vov_host = pvov_host.getViewHost();

      const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
      for (int cell=0; cell < num_cells; ++cell) {
        TEST_FLOATING_EQUALITY(vov_host(0,0,0)(cell),2.0,tol);
        TEST_FLOATING_EQUALITY(vov_host(1,1,1)(cell),3.0,tol);
        TEST_FLOATING_EQUALITY(vov_host(2,2,2)(cell),4.0,tol);
        TEST_FLOATING_EQUALITY(vov_host(0,1,2)(cell),12.0,tol);
      }
    }
  }

}

using ScalarType = Sacado::Fad::DFad<double>;

PHX::ViewOfViews<2,PHX::View<ScalarType**>,PHX::Device> createVoV()
{
  PHX::ViewOfViews<2,PHX::View<ScalarType**>,PHX::Device> tmp;

  tmp.initialize("tmp from createVoV()",2,2);
  const int num_cells = 10;
  const int num_pts = 5;
  const int num_deriv = 2;
  PHX::View<ScalarType**> a("a",num_cells,num_pts,num_deriv+1);
  PHX::View<ScalarType**> b("b",num_cells,num_pts,num_deriv+1);
  PHX::View<ScalarType**> c("c",num_cells,num_pts,num_deriv+1);
  // Don't assign the "d" array. Simulates a jagged outer array.
  // PHX::View<double*> d("d",5,num_deriv+1);
  tmp.addView(a,0,0);
  tmp.addView(b,0,1);
  tmp.addView(c,1,0);
  // tmp.addView(d,1,1);
  tmp.syncHostToDevice();

  return tmp;
}

template<typename VoVType>
void initializeVoV(VoVType& vov)
{
  auto Mat_h = vov.getViewHost();
  auto a = Mat_h(0,0);
  auto b = Mat_h(0,1);
  auto c = Mat_h(1,0);

  // Initialize a, b and c
  Kokkos::MDRangePolicy<PHX::Device::execution_space,Kokkos::Rank<2>> policy({0,0},{a.extent(0),a.extent(1)});
  Kokkos::parallel_for("FadAndAssignment init",policy,KOKKOS_LAMBDA(const int cell, const int pt) {
    a(cell,pt).val() = double(cell) + double(pt);
    a(cell,pt).fastAccessDx(0) = 0.0;
    a(cell,pt).fastAccessDx(1) = 2.0 * double(cell) + double(pt);
    b(cell,pt).val() = double(cell);
    b(cell,pt).fastAccessDx(0) = 0.0;
    b(cell,pt).fastAccessDx(1) = 0.0;
    c(cell,pt).val() = 0.0;
    c(cell,pt).fastAccessDx(0) = 0.0;
    c(cell,pt).fastAccessDx(1) = 0.0;
  });
  PHX::Device::execution_space().fence();
}

template<typename VoVType, typename OstreamType>
void testVoV(VoVType& vov, OstreamType& out, bool& success)
{
  // Compute c using device vov
  auto Mat_h = vov.getViewHost();
  auto c = Mat_h(1,0);
  auto Mat = vov.getViewDevice();
  const bool use_hierarchic = true;
  if (use_hierarchic) {
    Kokkos::TeamPolicy<PHX::exec_space> policy(c.extent(0),Kokkos::AUTO());
    Kokkos::parallel_for("FadAndAssignement compute",policy,KOKKOS_LAMBDA(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) {
      const auto cell = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,c.extent(1)), [&] (const int& pt) {
        Mat(1,0)(cell,pt) = Mat(0,0)(cell,pt) * Mat(0,1)(cell,pt);
      });
    });
  } else {
    Kokkos::MDRangePolicy<PHX::Device::execution_space,Kokkos::Rank<2>> policy({0,0},{c.extent(0),c.extent(1)});
    Kokkos::parallel_for("FadAndAssignement compute",policy,KOKKOS_LAMBDA(const int cell, const int pt) {
      Mat(1,0)(cell,pt) = Mat(0,0)(cell,pt) * Mat(0,1)(cell,pt);
    });
  }
  PHX::Device::execution_space().fence();

  // Check the results
  auto c_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),c);
  const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
  for (size_t cell=0; cell < c.extent(0); ++cell) {
    for (size_t pt=0; pt < c.extent(1); ++pt) {
      // printf("c(%zu,%zu).val=%f, c.dx(0)=%f, c.dx(1)=%f\n",cell,pt,c_h(cell,pt).val(),c_h(cell,pt).fastAccessDx(0),c_h(cell,pt).fastAccessDx(1));
      double gold_val = double(cell) * (double(cell) + double(pt));
      TEST_FLOATING_EQUALITY( c_h(cell,pt).val(),gold_val,tol);
      double gold_dx1 = (2.0 * double(cell) + double(pt)) * double(cell);
      TEST_FLOATING_EQUALITY( c_h(cell,pt).fastAccessDx(1),gold_dx1,tol);
    }
  }
}

// This test checks the copy and assignment objects for a vov. We have
// many use cases where the vov is created in one object and then
// passed around and stored in other objects. The assignment operators
// loop over the test multiple times to prove that both an empty vov
// and a perviously populated vov can both be assigned to a new
// vov. The previously populated version will hang due to nested
// parallel_for in the kokkos view dtors if not handled corectly.
TEUCHOS_UNIT_TEST(PhalanxViewOfViews,CtorAndAssignmentWithFadData) {
  using VoV = PHX::ViewOfViews<2,PHX::View<ScalarType**>,PHX::Device>;

  // Test the move assignment operator
  VoV vov;
  for (int i=0; i < 2; ++i) {
    vov = std::move(createVoV());
    initializeVoV(vov);
    testVoV(vov, out, success);
  }

  // Test the copy constructor
  {
    VoV vov2(vov);
    initializeVoV(vov2);
    testVoV(vov2, out, success);
  }

  // Test the move contructor
  {
    VoV vov3(std::move(createVoV()));
    initializeVoV(vov3);
    testVoV(vov3, out, success);
  }

  // Test the copy assignment operator
  VoV vov4;
  for (int i=0; i < 2; ++i) {
    vov4 = vov;
    initializeVoV(vov4);
    testVoV(vov4, out, success);
  }
}

// This unit test is to check view-of-views with fad data types. There
// is an issue with sacado where MDRange does not index into the
// derivative array correctly if HIERARCHIC parallelism in sacado is
// enabled. Flat and Team kernels work fine. MDRange is using the
// hierarchic parallelism for the calculation and leaves no
// parallelism for sacado to use on the hidden fad dimension. Team
// uses hierarchic, but leaves the vector level parallelism for sacado
// (if we write the kernels correctly).

enum class PHX_KERNEL_TYPE {TEAM,MDRANGE,FLAT};
TEUCHOS_UNIT_TEST(PhalanxViewOfViews,FadHierarchicMDRangeBug) {
  using ST = Sacado::Fad::DFad<double>;
  using VT = PHX::View<ST**>; // PHX::View is a Kokkos::View with Contiguous layout for FAD types.
  const int num_cell = 2;
  const int num_pt = 3;
  const int num_deriv = 2+1;
  VT a("a",num_cell,num_pt,num_deriv);
  VT b("b",num_cell,num_pt,num_deriv);

  {
    Kokkos::MDRangePolicy<PHX::Device::execution_space,Kokkos::Rank<2>> policy({0,0},{a.extent(0),a.extent(1)});
    Kokkos::parallel_for("FadRawPtr init",policy,KOKKOS_LAMBDA(const int cell, const int pt) {
      a(cell,pt).val() = double(cell) + double(pt);
      a(cell,pt).fastAccessDx(1) = 2.0 + double(cell) + double(pt);
    });
    PHX::exec_space().fence(); // don't need this but being safe for debugging
  }

  std::cout << std::endl;
  const PHX_KERNEL_TYPE kernel = PHX_KERNEL_TYPE::TEAM; // This passes
  // const PHX_KERNEL_TYPE kernel = PHX_KERNEL_TYPE::MDRANGE; // This fails for FAD scalar types
  // const PHX_KERNEL_TYPE kernel = PHX_KERNEL_TYPE::FLAT; // This passes
  if (kernel == PHX_KERNEL_TYPE::TEAM) {
    Kokkos::TeamPolicy<PHX::exec_space> policy(num_cell,Kokkos::AUTO());
    Kokkos::parallel_for("FadRawPtr b=a*a",policy,KOKKOS_LAMBDA(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) {
      const auto cell = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,a.extent(1)), [&] (const int& pt) {
        b(cell,pt) = a(cell,pt) * a(cell,pt);
        // printf("DEVICE: b(%d,%d).val()=%f, b.dx(1)=%f\n",cell,pt,b(cell,pt).val(),b(cell,pt).fastAccessDx(1));
      });
    });
  } else if (kernel == PHX_KERNEL_TYPE::MDRANGE) {
    Kokkos::MDRangePolicy<PHX::Device::execution_space,Kokkos::Rank<2>> policy({0,0},{a.extent(0),a.extent(1)});
    Kokkos::parallel_for("FadRawPtr b=a*a",policy,KOKKOS_LAMBDA(const int cell, const int pt) {
      b(cell,pt) = a(cell,pt) * a(cell,pt);
      // printf("DEVICE: b(%d,%d).val()=%f, b.dx(1)=%f\n",cell,pt,b(cell,pt).val(),b(cell,pt).fastAccessDx(1));
    });
  } else {
    Kokkos::parallel_for("FadRawPtr b=a*a",a.extent(0),KOKKOS_LAMBDA(const int cell) {
      for (size_t pt=0; pt < a.extent(1); ++pt) {
        b(cell,pt) = a(cell,pt) * a(cell,pt);
        // printf("DEVICE: b(%d,%d).val()=%f, b.dx(1)=%f\n",cell,pt,b(cell,pt).val(),b(cell,pt).fastAccessDx(1));
      }
    });
  }

  PHX::exec_space().fence();

  auto b_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),b);
  const auto tol = std::numeric_limits<double>::epsilon() * 100.0;
  for (size_t cell=0; cell < a.extent(0); ++cell) {
    for (size_t pt=0; pt < a.extent(1); ++pt) {
      out << "b_h(" << cell << "," << pt << ")\n";
      double gold_val = (double(cell) + double(pt)) * (double(cell) + double(pt));
      TEST_FLOATING_EQUALITY( b_h(cell,pt).val(),gold_val,tol);
      double gold_dx1 = 2.0 * (2.0 + double(cell) + double(pt)) * (double(cell) + double(pt));
      TEST_FLOATING_EQUALITY( b_h(cell,pt).fastAccessDx(1),gold_dx1,tol);
    }
  }
}


// This test checks that the safety check works. Since the safety
// check occurs in the deleter, the error calls abort instead of
// throwning a exception. We would need to switch the test harness to
// gtest to be able to test aborts in unit tests, so the needed code
// check is commented out for now. You can manually run this check by
// uncommenting the line below.

TEUCHOS_UNIT_TEST(PhalanxViewOfViews,ViewOfView_SafetyCheckAbort) {

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
    PHX::ViewOfViews<OuterViewRank,InnerView,mem_t> v_of_v;

    TEST_ASSERT(!v_of_v.isInitialized());
    v_of_v.initialize("outer host",2,2);
    TEST_ASSERT(v_of_v.isInitialized());

    v_of_v.addView(a,0,0);
    v_of_v.addView(b,0,1);
    v_of_v.addView(c,1,0);
    v_of_v.addView(d,1,1);

    TEST_ASSERT(!v_of_v.deviceViewIsSynced());
    v_of_v.syncHostToDevice();
    TEST_ASSERT(v_of_v.deviceViewIsSynced());

    auto v_dev = v_of_v.getViewDevice();
    auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0},{num_cells,num_pts,num_equations});
    Kokkos::parallel_for("view of view test",policy,KOKKOS_LAMBDA (const int cell,const int pt, const int eq) {
      v_dev(1,1)(cell,pt,eq) = v_dev(0,0)(cell,pt,eq) + v_dev(0,1)(cell,pt,eq) + v_dev(1,0)(cell,pt,eq);
    });

    // Uncomment the line below to prove the ViewOfViews prevents
    // device views from outliving host view. This line will cause a
    // Kokkos::abort() and error message since v_dev above is still in
    // scope when the ViewOfViews is destroyed.
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
