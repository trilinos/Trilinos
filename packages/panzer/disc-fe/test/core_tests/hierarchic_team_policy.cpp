// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>

#include "KokkosExp_View_Fad.hpp"
#include "Kokkos_DynRankView_Fad.hpp"
#include "Panzer_HierarchicParallelism.hpp"
#include "Sacado.hpp"

namespace panzer_test {

  const int M = 100;
  const int N = 16;

  template<typename Scalar, typename VectorType,typename OutputStream>
  void checkPolicy(bool use_stream_instance,
                   VectorType& a, VectorType& b, VectorType& c,
                   bool& success, OutputStream& out)
  {
    Kokkos::deep_copy(a,0.0);
    Kokkos::deep_copy(b,1.0);
    Kokkos::deep_copy(c,2.0);

    if (use_stream_instance) {
      PHX::ExecSpace exec_space;
      auto policy = panzer::HP::inst().teamPolicy<Scalar>(exec_space,M);
      Kokkos::parallel_for("test 0",policy,KOKKOS_LAMBDA (const Kokkos::TeamPolicy<PHX::ExecSpace>::member_type team){
        const int i = team.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,N), [&] (const int j) {
          a(i,j) += b(i,j) + c(i,j);
        });
      });
    }
    else {
      auto policy = panzer::HP::inst().teamPolicy<Scalar>(M);
      Kokkos::parallel_for("test 0",policy,KOKKOS_LAMBDA (const Kokkos::TeamPolicy<PHX::ExecSpace>::member_type team){
        const int i = team.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,N), [&] (const int j) {
          a(i,j) += b(i,j) + c(i,j);
        });
      });
    }
    Kokkos::fence();

    auto a_host = Kokkos::create_mirror_view(a);
    Kokkos::deep_copy(a_host,a);
    auto tol = 1000.0 * std::numeric_limits<double>::epsilon();

    for (int i=0; i < M; ++i) {
      for (int j=0; j < N; ++j) {
        TEST_FLOATING_EQUALITY(Sacado::scalarValue(a_host(i,j)),3.0,tol);
      }
    }
  }

  TEUCHOS_UNIT_TEST(HierarchicTeamPolicy, StreamsDouble)
  {
    using Scalar = double;
    PHX::View<Scalar**> a("a",M,N);
    PHX::View<Scalar**> b("b",M,N);
    PHX::View<Scalar**> c("c",M,N);
    panzer_test::checkPolicy<Scalar>(false,a,b,c,success,out); // default exec space
    panzer_test::checkPolicy<Scalar>(true,a,b,c,success,out);  // specify exec space
  }

  TEUCHOS_UNIT_TEST(HierarchicTeamPolicy, StreamsDFAD)
  {
    using Scalar = Sacado::Fad::DFad<double>;
    const int deriv_dim = 8;
    PHX::View<Scalar**> a("a",M,N,deriv_dim);
    PHX::View<Scalar**> b("b",M,N,deriv_dim);
    PHX::View<Scalar**> c("c",M,N,deriv_dim);
    panzer_test::checkPolicy<Scalar>(false,a,b,c,success,out); // default exec space
    panzer_test::checkPolicy<Scalar>(true,a,b,c,success,out);  // specify exec space
  }

}
