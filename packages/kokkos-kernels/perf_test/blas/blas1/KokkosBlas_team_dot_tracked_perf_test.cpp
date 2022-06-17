/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Core.hpp>
#include <KokkosBlas1_team_dot.hpp>
#include <Kokkos_Random.hpp>
// For RPS implementation
#include "KokkosBlas_team_dot_perf_test.hpp"
#ifdef KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE
#include <PerfTestUtilities.hpp>
#endif  // KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE

template <class ExecSpace, class Layout>
testData_rps_team_dot<ExecSpace, Layout> setup_test(int m, int repeat,
                                                    const int numberOfTeams) {
  // use constructor to generate testData_object
  testData_rps_team_dot<ExecSpace, Layout> testData_rps_team_dot_obj(m);

  // set a field in the struct
  testData_rps_team_dot_obj.m             = m;
  testData_rps_team_dot_obj.repeat        = repeat;
  testData_rps_team_dot_obj.numberOfTeams = numberOfTeams;

  return testData_rps_team_dot_obj;
}

test_list construct_team_dot_kernel_base(const rajaperf::RunParams& run_params)

{
  // instantiate test_list (the type) as kernel_base_vector
  test_list kernel_base_vector;

  /////////////////////////////////////////////////////////////////////////////
  // https://github.com/kokkos/kokkos-kernels/wiki/BLAS-1::team-dot
  /////////////////////////////////////////////////////////////////////////////

  using test_data_type = decltype(
      setup_test<Kokkos::DefaultExecutionSpace,
                 Kokkos::DefaultExecutionSpace::array_layout>(1, 1, 1));

  kernel_base_vector.push_back(rajaperf::make_kernel_base(
      "BLAS_TEAM_DOT ", run_params,
      [=](const int m, const int repeat) {
        // returns a tuple of testData_obj
        return std::make_tuple(
            // TODO: Discuss decltype
            // TODO: Ask KK what values they want tested?
            setup_test<Kokkos::DefaultExecutionSpace,
                       Kokkos::DefaultExecutionSpace::array_layout>(m, repeat,
                                                                    1));
      },
      [&](const int, const int, test_data_type& data) {
        Kokkos::parallel_for(
            "TeamDotUsage_RPS",
            test_data_type::policy(data.numberOfTeams, Kokkos::AUTO),
            KOKKOS_LAMBDA(const test_data_type::member_type& team) {
              // loop body
              KokkosBlas::Experimental::dot(team, data.x, data.y);
            });
      }));

  // Overall return - a vector of kernel base objects
  // containing data for setup and run lambdas of type test_list

  return kernel_base_vector;
}
