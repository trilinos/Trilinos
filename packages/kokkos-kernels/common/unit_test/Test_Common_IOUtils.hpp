// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
  using scalar_t   = KokkosKernels::default_scalar;
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

TEST_F(TestCategory, common_print_view) { testPrintView<TestDevice>(); }

#endif  // KOKKOSKERNELS_IOTEST_HPP
