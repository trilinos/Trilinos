// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

//
// TEMPORARY diagnostic for the SLFad failures on SYCL.  Delete once they are
// understood.
//
// Every failing test in both SYCL suites performs one operation in common:
// assigning a scalar into a *device* View of Fad.
//
//   Fad_KokkosTests      ValueAssign          Kokkos::deep_copy(a, 2.3456)
//   Fad_KokkosTests      AtomicAdd            Kokkos::deep_copy(v, 2.3456)
//   Fad_KokkosTests      LocalDeepCopy[Team]  Kokkos::Experimental::local_deep_copy
//   Fad_KokkosAtomicTests  (all of them)      Kokkos::deep_copy(s, tag.init())
//
// Those are also the ONLY tests in Fad_KokkosTests that touch a Fad view that
// way -- every other deep_copy there targets Sacado::as_scalar_view(), a plain
// double view -- which is why 16 tests fail in one suite and all 48 SLFad
// tests fail in the other.
//
// This strips the Teuchos harness away and does only that, for SFad (which
// passes) and SLFad (which does not), labelling every step so the failing one
// is obvious.  It also reports the two build settings that turn a Sacado debug
// check into a hard device failure:
//
//   SACADO_DEBUG  -- enables the size checks in ExprAssign/StaticStorage.  If
//                    they are reachable from device code they emit a
//                    Kokkos::abort(), which on SYCL is a device-side printf.
//   NDEBUG        -- when NOT set, Kokkos::abort() on SYCL calls __assert_fail
//                    (Kokkos_SYCL_Abort.hpp), which surfaces as an opaque
//                    UR_RESULT_ERROR_UNKNOWN rather than a message
//
// Run plain, then with SYCL_UR_TRACE=1 to see the underlying Level Zero error.
//

#include "Sacado.hpp"
#include "Sacado_Fad_Kokkos.hpp"

#include <cstdio>
#include <exception>

typedef Sacado::Fad::SFad<double, 5> SFadType;    // static_size = 5
typedef Sacado::Fad::SLFad<double, 10> SLFadType; // static_size = 0, max 10

#if defined(KOKKOS_ENABLE_SYCL)
typedef Kokkos::SYCL exec_space;
#else
typedef Kokkos::DefaultExecutionSpace exec_space;
#endif

namespace {

int g_failures = 0;

// Run one labelled step, reporting whether it threw rather than letting the
// exception escape, so later steps still run.
template <typename F> void step(const char *what, F f) {
  std::printf("    %-42s ", what);
  std::fflush(stdout);
  try {
    f();
    Kokkos::fence();
    std::printf("ok\n");
  } catch (const std::exception &e) {
    std::printf("THREW\n      %s\n", e.what());
    ++g_failures;
  }
  std::fflush(stdout);
}

template <typename FadType, typename Layout>
void probe(const char *fad_name, const char *layout_name) {
  typedef Kokkos::View<FadType *, Layout, exec_space> ViewType;
  typedef Kokkos::View<FadType **, Layout, exec_space> View2Type;
  typedef Kokkos::View<FadType, Layout, exec_space> ScalarViewType;

  const int num_rows = 11;
  const int num_cols = 7;
  const int fad_size = 5;

  std::printf("\n  %s / %s\n", fad_name, layout_name);
  std::printf("    sizeof(%s) = %zu bytes; view allocates %d doubles = %zu "
              "bytes per element\n",
              fad_name, sizeof(FadType), fad_size + 1,
              (fad_size + 1) * sizeof(double));

  ViewType v;
  View2Type v2;
  ScalarViewType s;

  step("allocate views", [&] {
    v = ViewType("v", num_rows, fad_size + 1);
    v2 = View2Type("v2", num_rows, num_cols, fad_size + 1);
    s = ScalarViewType("s", fad_size + 1);
  });

  // The operation shared by every failing test.
  step("deep_copy(rank-0 Fad view, scalar)",
       [&] { Kokkos::deep_copy(s, 2.3456); });

  step("deep_copy(rank-1 Fad view, scalar)",
       [&] { Kokkos::deep_copy(v, 2.3456); });

  // For contrast: a plain kernel touching the same views.  If this succeeds
  // while the deep_copy above fails, the Fad view itself is sound and the
  // fault is specific to the scalar fill.
  step("parallel_for writing v(i).val()", [&] {
    ViewType vl = v;
    Kokkos::parallel_for(
        Kokkos::RangePolicy<exec_space>(0, num_rows),
        KOKKOS_LAMBDA(const int i) { vl(i).val() = 1.0 * i; });
  });

  step("parallel_for writing v(i).fastAccessDx(j)", [&] {
    ViewType vl = v;
    Kokkos::parallel_for(
        Kokkos::RangePolicy<exec_space>(0, num_rows),
        KOKKOS_LAMBDA(const int i) {
          for (int j = 0; j < fad_size; ++j)
            vl(i).fastAccessDx(j) = 7.89 + j;
        });
  });

  // The device-side sibling used by the LocalDeepCopy tests.
  step("local_deep_copy(subview, Fad value)", [&] {
    View2Type vl = v2;
    ScalarViewType sl = s;
    Kokkos::parallel_for(
        Kokkos::RangePolicy<exec_space>(0, num_rows),
        KOKKOS_LAMBDA(const int i) {
          auto row = Kokkos::subview(vl, i, Kokkos::ALL);
          Kokkos::Experimental::local_deep_copy(row, sl());
        });
  });

  // Fad-to-Fad assignment on device, which is where the SACADO_DEBUG size
  // check in ExprAssign lives.
  step("parallel_for assigning v(i) = s()", [&] {
    ViewType vl = v;
    ScalarViewType sl = s;
    Kokkos::parallel_for(
        Kokkos::RangePolicy<exec_space>(0, num_rows),
        KOKKOS_LAMBDA(const int i) { vl(i) = sl(); });
  });
}

} // namespace

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    std::printf("\n======== Sacado Fad deep_copy repro ========\n\n");

    std::printf("Build configuration:\n");
#if defined(SACADO_DEBUG)
    std::printf("  SACADO_DEBUG : DEFINED -- the size checks in ExprAssign and\n"
                "                 Static{,Fixed}Storage are active on the host;\n"
                "                 they are guarded out of device code\n");
#else
    std::printf("  SACADO_DEBUG : not defined -- Sacado's debug size checks "
                "are compiled out\n");
#endif
#if defined(NDEBUG)
    std::printf("  NDEBUG       : defined -- Kokkos::abort() on SYCL prints a "
                "message and continues\n");
#else
    std::printf("  NDEBUG       : NOT DEFINED -- Kokkos::abort() on SYCL calls "
                "__assert_fail(),\n                 which shows up as an opaque "
                "backend error\n");
#endif
    std::printf("  exec space   : %s\n", exec_space::name());

    using Kokkos::LayoutLeft;
    typedef Sacado::LayoutContiguous<Kokkos::LayoutLeft> LeftContiguous;

    probe<SFadType, LayoutLeft>("SFad<double,5>", "LayoutLeft");
    probe<SLFadType, LayoutLeft>("SLFad<double,10>", "LayoutLeft");
    probe<SFadType, LeftContiguous>("SFad<double,5>", "LayoutContiguous");
    probe<SLFadType, LeftContiguous>("SLFad<double,10>", "LayoutContiguous");

    std::printf("\n%d step(s) threw\n", g_failures);
  }
  Kokkos::finalize();

  if (g_failures == 0)
    std::printf("End Result: TEST PASSED\n");
  else
    std::printf("End Result: TEST FAILED\n");

  return g_failures == 0 ? 0 : 1;
}
