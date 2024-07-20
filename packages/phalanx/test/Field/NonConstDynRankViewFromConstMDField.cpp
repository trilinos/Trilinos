// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Phalanx_ExtentTraits.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Phalanx_GetNonConstDynRankViewFromConstMDField.hpp"

// This test demonstrates how to get a nonconst DynRankView from a
// const MDField for double and FAD scalar types.
PHX_EXTENT(CELL)
PHX_EXTENT(QP)
PHX_EXTENT(EQ)

using exec_t = Kokkos::DefaultExecutionSpace;
using mem_t = Kokkos::DefaultExecutionSpace::memory_space;

template<typename Scalar>
using non_const_mdfield = PHX::MDField<Scalar,CELL,QP,EQ>;

template<typename Scalar>
using const_mdfield = const PHX::MDField<const Scalar,CELL,QP,EQ>;

const int num_cells = 10;
const int num_pts = 8;
const int num_equations = 32;
const int num_derivatives = 4;

namespace {
  // function input for a,b,c are all CONST
  template<typename Scalar>
  void runTest(std::ostream& out,
               bool& success,
               const_mdfield<Scalar>& a,
               const_mdfield<Scalar>& b,
               const_mdfield<Scalar>& c)
  {
    auto tmp_a = a.get_static_view();
    auto tmp_b = b.get_static_view();
    // Demonstrate getting a nonconst view from a const view
    auto tmp_c = PHX::getNonConstDynRankViewFromConstMDField<Scalar>(c);

    Kokkos::MDRangePolicy<exec_t,Kokkos::Rank<3>> policy({0,0,0},{num_cells,num_pts,num_equations});
    Kokkos::parallel_for("use non-const DynRankView from const View",policy,KOKKOS_LAMBDA (const int cell,const int pt,const int eq) {
      tmp_c(cell,pt,eq) = tmp_a(cell,pt,eq) + tmp_b(cell,pt,eq);
    });
    exec_t().fence();
  }
}

// Check double values to make sure layout hasn't been switched.
TEUCHOS_UNIT_TEST(NonConstDynRankViewFromView,double) {
  using ScalarType = double; 
  non_const_mdfield<ScalarType> a("a","layout",num_cells,num_pts,num_equations);
  non_const_mdfield<ScalarType> b("b","layout",num_cells,num_pts,num_equations);
  non_const_mdfield<ScalarType> c("c","layout",num_cells,num_pts,num_equations);
  Kokkos::deep_copy(a.get_static_view(),2.0);
  Kokkos::deep_copy(b.get_static_view(),3.0);
  runTest<ScalarType>(out,success,a,b,c);

  auto c_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),c.get_static_view());
  auto tol = 100.0 * Teuchos::ScalarTraits<ScalarType>::eps();
  for (int cell=0; cell < num_cells; ++cell)
    for (int pt=0; pt < num_pts; ++pt)
      for (int eq=0; eq < num_equations; ++eq) {
        TEST_FLOATING_EQUALITY(c_host(cell,pt,eq),5.0,tol);
      }
}

// Check DFAD values to make sure layout hasn't been switched.
TEUCHOS_UNIT_TEST(NonConstDynRankViewFromView,FAD) {
  using RealType = double;
  using ScalarType = Sacado::Fad::DFad<RealType>;
  non_const_mdfield<ScalarType> a("a","layout",num_cells,num_pts,num_equations,num_derivatives);
  non_const_mdfield<ScalarType> b("b","layout",num_cells,num_pts,num_equations,num_derivatives);
  non_const_mdfield<ScalarType> c("c","layout",num_cells,num_pts,num_equations,num_derivatives);
  auto tmp_a = a.get_static_view();
  auto tmp_b = b.get_static_view();
  Kokkos::MDRangePolicy<exec_t,Kokkos::Rank<3>> policy({0,0,0},{num_cells,num_pts,num_equations});
  Kokkos::parallel_for("initialize fads",policy,KOKKOS_LAMBDA (const int cell,const int pt,const int eq) {
    tmp_a(cell,pt,eq).val() = 2.0;
    tmp_b(cell,pt,eq).val() = 3.0;
  });
  exec_t().fence();

  runTest<ScalarType>(out,success,a,b,c);

  auto c_host = Kokkos::create_mirror_view(c.get_static_view());
  Kokkos::deep_copy(c_host,c.get_static_view());
  auto tol = 100.0 * Teuchos::ScalarTraits<RealType>::eps();
  for (int cell=0; cell < num_cells; ++cell)
    for (int pt=0; pt < num_pts; ++pt)
      for (int eq=0; eq < num_equations; ++eq) {
        TEST_FLOATING_EQUALITY(c_host(cell,pt,eq).val(),5.0,tol);
      }
}
