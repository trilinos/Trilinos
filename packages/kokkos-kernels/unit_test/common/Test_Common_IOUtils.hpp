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

/// \file Test_Common_IOUtils.hpp
/// \brief Tests for IO and print routines

#ifndef KOKKOSKERNELS_IOTEST_HPP
#define KOKKOSKERNELS_IOTEST_HPP

#include <Kokkos_Core.hpp>
#include <KokkosKernels_PrintUtils.hpp>

template <typename OutStream, typename ExecSpace>
class ViewPrintHelper {
 public:
  ViewPrintHelper(OutStream &out_, const char *sep_) : out(out_), sep(sep_) {}

  template <typename T>
  void operator()(T view, int limit = 0) {
    const auto v = Kokkos::create_mirror_view_and_copy(space, view);
    kk_print_1Dview(out, v, limit < 1, sep, limit);
  }

 private:
  OutStream &out;
  ExecSpace space;
  const char *sep;
};

template <typename exec_space>
void testPrintView() {
  using scalar_t   = default_scalar;
  using Unmanaged  = Kokkos::MemoryTraits<Kokkos::Unmanaged>;
  using rank0_view = Kokkos::View<scalar_t, Kokkos::HostSpace, Unmanaged>;
  using rank1_view = Kokkos::View<scalar_t *, Kokkos::HostSpace, Unmanaged>;
  using rank2_view = Kokkos::View<scalar_t **, Kokkos::HostSpace, Unmanaged>;

  // Note: try custom separator
  std::stringstream out;
  auto test = ViewPrintHelper<decltype(out), exec_space>(out, "|");

  std::vector<scalar_t> vals = {4, 5, 6, 7};
  const rank1_view hv1(vals.data(), vals.size());
  test(hv1);
  test(hv1, 4);
  test(hv1, 3);
  test(rank0_view(vals.data()));
  test(rank2_view(vals.data(), vals.size(), 1));
  test(rank2_view(vals.data(), vals.size() / 2, 2));

  EXPECT_EQ(out.str(),
            "4|5|6|7|\n"
            "4|5|6|7|\n"
            "4|... ... ...|7|\n"
            "4|\n"
            "4|5|6|7|\n"
            "[2x2 multi-vector]\n");
}

TEST_F(TestCategory, common_print_view) { testPrintView<TestExecSpace>(); }

#endif  // KOKKOSKERNELS_IOTEST_HPP
