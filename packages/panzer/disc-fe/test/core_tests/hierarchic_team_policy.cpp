// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>

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
        TEST_FLOATING_EQUALITY(Sacado::ScalarValue<Scalar>::eval(a(i,j)),3.0,tol);
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
