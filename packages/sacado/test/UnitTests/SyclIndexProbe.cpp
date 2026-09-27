// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

//
// Diagnostic probe for Sacado's SYCL port.  This uses only Kokkos and raw
// SYCL -- no Sacado headers -- and answers the questions that decide how the
// hierarchical (vector-partitioned Fad) code should obtain thread indices on
// SYCL:
//
//   1. What is the ACTUAL vector width of a Kokkos TeamPolicy on this device?
//      Kokkos clamps the requested vector size to the device's maximum
//      sub-group size (Kokkos_SYCL_ParallelFor_Team.hpp), so a request for 32
//      may silently become 16.  Sacado's compile-time partition stride comes
//      from the layout (Sacado::LayoutContiguous<Layout,Stride>) and must
//      agree with it, or partitioned views hand out wrong strides.
//
//   2. Is `get_nd_item<2>().get_local_id(1)` really the vector lane -- i.e.
//      the portable equivalent of CUDA's threadIdx.x?  Kokkos maps
//      team_rank() to get_local_id(0) and iterates ThreadVectorRange over
//      get_local_id(1) (Kokkos_SYCL_Team.hpp), so this should hold, but it is
//      the assumption the whole hierarchical port rests on.
//
//   3. Do the oneAPI free-function work-item queries even compile and run in
//      the kernel shapes Sacado reaches them from -- team kernels (nd_range)
//      and flat kernels (which Kokkos may launch as a basic sycl::range)?
//
//   4. Do get_sub_group() and group_ballot() work inside a flat RangePolicy
//      kernel?  That is exactly what the Phase 1 Fad atomics need.
//
// The probe is deliberately built out of two independent measurements of the
// vector width: a portable work-item count that needs no SYCL calls at all,
// and the nd_item query.  If the query turns out to be unusable, the portable
// measurement still answers question 1.
//
// Compile-time switches (all default to the useful setting):
//
//   -DSACADO_PROBE_ND_ITEM=0        disable the nd_item query in team kernels
//                                   (use if it fails to compile -- that
//                                   failure is itself the answer)
//   -DSACADO_PROBE_RANGE_SUBGROUP=0 disable the sub-group/ballot probe in flat
//                                   kernels
//   -DSACADO_PROBE_RANGE_ND_ITEM=1  opt in to calling get_nd_item<2>() from a
//                                   flat RangePolicy kernel.  This is
//                                   expected to be invalid; enable it only to
//                                   find out HOW it fails (compile error,
//                                   runtime error, or silent garbage).
//

#include "Teuchos_GlobalMPISession.hpp"

#include "Kokkos_Core.hpp"

#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#if !defined(KOKKOS_ENABLE_SYCL)

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  std::printf("SyclIndexProbe: Kokkos was not built with the SYCL backend; "
              "nothing to probe.\n");
  std::printf("End Result: TEST PASSED\n");
  return 0;
}

#else

// The spelling of the oneAPI free-function work-item queries changed in the
// 2025.0 compilers.  Mirror the switch desul itself uses, see
// tpls/desul/include/desul/atomics/Lock_Based_Fetch_Op_SYCL.hpp.
#if defined(__INTEL_LLVM_COMPILER) && __INTEL_LLVM_COMPILER >= 20250000
#define SACADO_PROBE_ND_ITEM_2()                                               \
  sycl::ext::oneapi::this_work_item::get_nd_item<2>()
#define SACADO_PROBE_SUB_GROUP()                                               \
  sycl::ext::oneapi::this_work_item::get_sub_group()
#define SACADO_PROBE_QUERY_SPELLING "sycl::ext::oneapi::this_work_item::get_*"
#else
#define SACADO_PROBE_ND_ITEM_2() sycl::ext::oneapi::experimental::this_nd_item<2>()
#define SACADO_PROBE_SUB_GROUP() sycl::ext::oneapi::experimental::this_sub_group()
#define SACADO_PROBE_QUERY_SPELLING "sycl::ext::oneapi::experimental::this_*"
#endif

#ifndef SACADO_PROBE_ND_ITEM
#define SACADO_PROBE_ND_ITEM 1
#endif

#ifndef SACADO_PROBE_RANGE_SUBGROUP
#define SACADO_PROBE_RANGE_SUBGROUP 1
#endif

#ifndef SACADO_PROBE_RANGE_ND_ITEM
#define SACADO_PROBE_RANGE_ND_ITEM 0
#endif

typedef Kokkos::SYCL exec_space;
typedef Kokkos::View<int *, exec_space> result_view;
typedef Kokkos::View<int *, exec_space>::host_mirror_type host_result_view;

enum ResultSlot {
  R_WORK_ITEMS = 0,      // work items that executed the team body
  R_TEAM_SIZE,           // team.team_size()
  R_ND_RANGE_0,          // nd_item.get_local_range(0)
  R_ND_RANGE_1,          // nd_item.get_local_range(1)
  R_MAX_LOCAL_ID_0,      // max nd_item.get_local_id(0)
  R_MAX_LOCAL_ID_1,      // max nd_item.get_local_id(1)
  R_ID0_TEAM_RANK_BAD,   // count of work items where local_id(0) != team_rank()
  R_LANE_ID1_BAD,        // count where ThreadVectorRange lane != local_id(1)
  R_SUB_GROUP_SIZE,      // sub_group.get_local_range()[0]
  R_MAX_SG_LOCAL_ID,     // max sub_group.get_local_id()[0]
  R_ND_QUERY_RAN,        // 1 if the nd_item query executed on device
  R_BALLOT_RAN,          // 1 if group_ballot() executed and behaved sanely
  R_CHECKSUM,            // keeps the register-pressure work alive
  R_COUNT
};

namespace {

// ---------------------------------------------------------------------------
// Team kernel probe
// ---------------------------------------------------------------------------
host_result_view probe_team(int league_size, int team_size, int vector_size) {
  result_view r("sacado_sycl_probe", R_COUNT);
  Kokkos::deep_copy(r, 0);

  typedef Kokkos::TeamPolicy<exec_space> policy_type;
  policy_type policy(league_size, team_size, vector_size);

  Kokkos::parallel_for(
      "sacado_sycl_probe_team", policy,
      KOKKOS_LAMBDA(const policy_type::member_type &team) {
        // Kokkos executes the team body redundantly on every work item of the
        // team, vector lanes included.  Counting them measures the real
        // work-group size without any SYCL-specific call, which is why this
        // number is trustworthy even if the queries below are not.
        Kokkos::atomic_add(&r(R_WORK_ITEMS), 1);
        Kokkos::atomic_max(&r(R_TEAM_SIZE), static_cast<int>(team.team_size()));

#if SACADO_PROBE_ND_ITEM
#if defined(__SYCL_DEVICE_ONLY__)
        auto item = SACADO_PROBE_ND_ITEM_2();
        const int id0 = static_cast<int>(item.get_local_id(0));
        const int id1 = static_cast<int>(item.get_local_id(1));
        const int range0 = static_cast<int>(item.get_local_range(0));
        const int range1 = static_cast<int>(item.get_local_range(1));

        Kokkos::atomic_max(&r(R_ND_RANGE_0), range0);
        Kokkos::atomic_max(&r(R_ND_RANGE_1), range1);
        Kokkos::atomic_max(&r(R_MAX_LOCAL_ID_0), id0);
        Kokkos::atomic_max(&r(R_MAX_LOCAL_ID_1), id1);
        Kokkos::atomic_max(&r(R_ND_QUERY_RAN), 1);

        if (id0 != static_cast<int>(team.team_rank()))
          Kokkos::atomic_add(&r(R_ID0_TEAM_RANK_BAD), 1);

        auto sg = SACADO_PROBE_SUB_GROUP();
        Kokkos::atomic_max(&r(R_SUB_GROUP_SIZE),
                           static_cast<int>(sg.get_local_range()[0]));
        Kokkos::atomic_max(&r(R_MAX_SG_LOCAL_ID),
                           static_cast<int>(sg.get_local_id()[0]));

        // The decisive check for the port: ThreadVectorRange over exactly the
        // vector width gives each lane one iteration, and that iteration index
        // is what Sacado would use as threadIdx.x.  It must equal local_id(1).
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(team, static_cast<int>(range1)),
            [&](const int lane) {
              if (lane != id1)
                Kokkos::atomic_add(&r(R_LANE_ID1_BAD), 1);
            });
#endif // __SYCL_DEVICE_ONLY__
#endif // SACADO_PROBE_ND_ITEM
      });
  Kokkos::fence();

  host_result_view h = Kokkos::create_mirror_view(r);
  Kokkos::deep_copy(h, r);
  return h;
}

// ---------------------------------------------------------------------------
// Register-pressure probe
// ---------------------------------------------------------------------------
//
// Kokkos does not ask the device for the vector width it will deliver, it asks
// the *kernel*:
//
//   auto max_sg_size = kernel.get_info<
//       sycl::info::kernel_device_specific::max_sub_group_size>(q.get_device());
//   auto final_vector_size = std::min<int>(m_vector_size, max_sg_size);
//                                    -- Kokkos_SYCL_ParallelFor_Team.hpp
//
// On Intel GPUs the compiler picks a kernel's SIMD width from its resource
// usage, so a register-hungry kernel can be compiled narrower than the device
// maximum and Kokkos then silently reduces the vector dimension to match.  That
// matters for Sacado because the partitioned Fad type is sized at compile time
// from the layout stride while the loop bounds come from the runtime width: if
// they disagree, each thread writes more derivative components than its local
// Fad has room for.
//
// This sweeps the size of a live per-thread array and reports the vector width
// actually delivered, so the effect can be measured rather than assumed.
template <int N>
host_result_view probe_register_pressure(int league_size, int team_size,
                                         int vector_size) {
  result_view r("sacado_sycl_probe", R_COUNT);
  Kokkos::deep_copy(r, 0);

  typedef Kokkos::TeamPolicy<exec_space> policy_type;
  policy_type policy(league_size, team_size, vector_size);

  Kokkos::parallel_for(
      "sacado_sycl_probe_registers", policy,
      KOKKOS_LAMBDA(const policy_type::member_type &team) {
        Kokkos::atomic_add(&r(R_WORK_ITEMS), 1);
        Kokkos::atomic_max(&r(R_TEAM_SIZE), static_cast<int>(team.team_size()));

        // A per-thread array, written and read back so it cannot be optimized
        // away -- the rolling dependency keeps every element live at once,
        // which is what drives the register pressure.
        double a[N];
        for (int i = 0; i < N; ++i)
          a[i] = 1.0 + i + team.league_rank();
        for (int k = 0; k < 4; ++k)
          for (int i = 0; i < N; ++i)
            a[i] = a[i] * 1.000001 + a[(i + 1) % N];
        double sum = 0.0;
        for (int i = 0; i < N; ++i)
          sum += a[i];
        Kokkos::atomic_max(&r(R_CHECKSUM), static_cast<int>(sum) & 0xffff);

#if SACADO_PROBE_ND_ITEM
#if defined(__SYCL_DEVICE_ONLY__)
        auto item = SACADO_PROBE_ND_ITEM_2();
        Kokkos::atomic_max(&r(R_ND_RANGE_1),
                           static_cast<int>(item.get_local_range(1)));
        auto sg = SACADO_PROBE_SUB_GROUP();
        Kokkos::atomic_max(&r(R_SUB_GROUP_SIZE),
                           static_cast<int>(sg.get_local_range()[0]));
        Kokkos::atomic_max(&r(R_ND_QUERY_RAN), 1);
#endif
#endif
      });
  Kokkos::fence();

  host_result_view h = Kokkos::create_mirror_view(r);
  Kokkos::deep_copy(h, r);
  return h;
}

// Kokkos throws for a team size the backend cannot launch; fall back to a
// single-thread team so a small backend still reports something.
template <int N>
bool run_register_pressure(int league_size, int team_size, int vector_size,
                           host_result_view &out) {
  const int team_sizes[2] = {team_size, 1};
  for (int t = 0; t < 2; ++t) {
    try {
      out = probe_register_pressure<N>(league_size, team_sizes[t], vector_size);
      return true;
    } catch (const std::exception &) {
    }
  }
  return false;
}

// ---------------------------------------------------------------------------
// Flat (RangePolicy) kernel probe
// ---------------------------------------------------------------------------
host_result_view probe_range(int n, int chunk_size) {
  result_view r("sacado_sycl_probe", R_COUNT);
  Kokkos::deep_copy(r, 0);

  Kokkos::RangePolicy<exec_space> policy(0, n);
  if (chunk_size > 0)
    policy.set_chunk_size(chunk_size);

  Kokkos::parallel_for(
      "sacado_sycl_probe_range", policy, KOKKOS_LAMBDA(const int) {
        // Touch `r` unconditionally, OUTSIDE the __SYCL_DEVICE_ONLY__ guard.
        // KOKKOS_LAMBDA captures by [=], so implicit capture only picks up
        // what the body names.  If every use of `r` sat behind the guard, the
        // host pass would capture nothing while the device pass captured the
        // View, the two passes would disagree on sizeof(lambda), and DPC++
        // rejects that with "Unexpected kernel lambda size".  This also
        // doubles as a count of the work items that actually ran.
        Kokkos::atomic_add(&r(R_WORK_ITEMS), 1);

#if defined(__SYCL_DEVICE_ONLY__)

#if SACADO_PROBE_RANGE_SUBGROUP
        // This is what the Phase 1 Fad atomics do: acquire a sub-group and
        // ballot over it to avoid deadlocking lock-stepped work items.
        auto sg = SACADO_PROBE_SUB_GROUP();
        Kokkos::atomic_max(&r(R_SUB_GROUP_SIZE),
                           static_cast<int>(sg.get_local_range()[0]));
        Kokkos::atomic_max(&r(R_MAX_SG_LOCAL_ID),
                           static_cast<int>(sg.get_local_id()[0]));

        using sycl::ext::oneapi::group_ballot;
        using sycl::ext::oneapi::sub_group_mask;
        sub_group_mask active = group_ballot(sg, 1);
        sub_group_mask none = group_ballot(sg, 0);
        // A sane ballot has at least this work item set in `active` and
        // nothing set in `none`.
        if (active != none)
          Kokkos::atomic_max(&r(R_BALLOT_RAN), 1);
#endif

#if SACADO_PROBE_RANGE_ND_ITEM
        // Expected to be invalid -- Kokkos may launch a RangePolicy as a basic
        // sycl::range<1> kernel, and even the nd_range<1> form is the wrong
        // dimensionality for a <2> query.
        auto item = SACADO_PROBE_ND_ITEM_2();
        Kokkos::atomic_max(&r(R_ND_RANGE_0),
                           static_cast<int>(item.get_local_range(0)));
        Kokkos::atomic_max(&r(R_ND_RANGE_1),
                           static_cast<int>(item.get_local_range(1)));
        Kokkos::atomic_max(&r(R_MAX_LOCAL_ID_1),
                           static_cast<int>(item.get_local_id(1)));
        Kokkos::atomic_max(&r(R_ND_QUERY_RAN), 1);
#endif

#endif // __SYCL_DEVICE_ONLY__
      });
  Kokkos::fence();

  host_result_view h = Kokkos::create_mirror_view(r);
  Kokkos::deep_copy(h, r);
  return h;
}

void print_device_info() {
  sycl::device dev = exec_space().sycl_queue().get_device();

  std::printf("  name                  : %s\n",
              dev.get_info<sycl::info::device::name>().c_str());
  std::printf("  max_work_group_size   : %zu\n",
              dev.get_info<sycl::info::device::max_work_group_size>());
  std::printf("  max_num_sub_groups    : %u\n",
              dev.get_info<sycl::info::device::max_num_sub_groups>());

  std::vector<size_t> sg_sizes =
      dev.get_info<sycl::info::device::sub_group_sizes>();
  std::printf("  sub_group_sizes       : [");
  for (size_t i = 0; i < sg_sizes.size(); ++i)
    std::printf("%s%zu", i ? ", " : "", sg_sizes[i]);
  std::printf("]\n");
}

} // namespace

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::InitializationSettings init_args;
  init_args.set_device_id(0);
  Kokkos::initialize(init_args);

  bool success = true;
  {
    std::printf("\n======== Sacado SYCL index probe ========\n\n");

    std::printf("Build configuration:\n");
#if defined(__INTEL_LLVM_COMPILER)
    std::printf("  __INTEL_LLVM_COMPILER : %d\n", __INTEL_LLVM_COMPILER);
#else
    std::printf("  __INTEL_LLVM_COMPILER : not defined\n");
#endif
    std::printf("  KOKKOS_VERSION        : %d\n", KOKKOS_VERSION);
    std::printf("  query spelling        : %s\n", SACADO_PROBE_QUERY_SPELLING);
#if defined(DESUL_SYCL_DEVICE_GLOBAL_SUPPORTED)
    std::printf("  DESUL_SYCL_DEVICE_GLOBAL_SUPPORTED : defined  (Fad atomics "
                "can use desul's SYCL locks)\n");
#else
    std::printf("  DESUL_SYCL_DEVICE_GLOBAL_SUPPORTED : NOT DEFINED  -- "
                "desul's SYCL lock_address is a stub that returns true, so "
                "Fad atomics would race silently.  Configure with "
                "Kokkos_ARCH_INTEL_PVC (or similar) before trusting them.\n");
#endif
    std::printf("  probes enabled        : nd_item(team)=%d, "
                "subgroup(range)=%d, nd_item(range)=%d\n",
                SACADO_PROBE_ND_ITEM, SACADO_PROBE_RANGE_SUBGROUP,
                SACADO_PROBE_RANGE_ND_ITEM);

    std::printf("\nDevice:\n");
    print_device_info();

    // ---------------------------------------------------------------------
    // Vector width sweep
    // ---------------------------------------------------------------------
    const int league_size = 8;
    const int team_size = 4;
    const int requested[] = {1, 2, 4, 8, 16, 32, 64};
    const int n_requested = sizeof(requested) / sizeof(requested[0]);

    std::printf("\nTeamPolicy(league=%d, team_size=%d, vector_size=V):\n",
                league_size, team_size);
    std::printf("  V requested | work items/team | team_size() | "
                "local_range(0) | local_range(1) | max local_id(1) | sg size\n");

    int width_for_32 = -1;
    int range1_for_32 = -1;
    int lane_bad_total = 0;
    int id0_bad_total = 0;
    bool nd_query_ran = false;

    for (int i = 0; i < n_requested; ++i) {
      const int v = requested[i];

      // Kokkos throws for a team_size/vector_size combination the backend
      // cannot launch (too large a work group, for instance).  Fall back to a
      // single-thread team, and if that is rejected too just report the row as
      // unavailable -- losing one row is better than losing the whole sweep.
      host_result_view h;
      bool launched = false;
      std::string launch_error;
      const int team_sizes[2] = {team_size, 1};
      for (int t = 0; t < 2 && !launched; ++t) {
        try {
          h = probe_team(league_size, team_sizes[t], v);
          launched = true;
        } catch (const std::exception &e) {
          launch_error = e.what();
        }
      }
      if (!launched) {
        std::printf("  %11d | launch rejected: %s\n", v, launch_error.c_str());
        continue;
      }

      // Work items per team, measured without any SYCL call.
      const int per_team = h(R_WORK_ITEMS) / (league_size * h(R_TEAM_SIZE));
      const int actual_vector_width = per_team;

      std::printf("  %11d | %15d | %11d | %14d | %14d | %15d | %7d%s\n",
                  v, actual_vector_width, h(R_TEAM_SIZE), h(R_ND_RANGE_0),
                  h(R_ND_RANGE_1), h(R_MAX_LOCAL_ID_1), h(R_SUB_GROUP_SIZE),
                  (actual_vector_width != v) ? "   <-- CLAMPED" : "");

      lane_bad_total += h(R_LANE_ID1_BAD);
      id0_bad_total += h(R_ID0_TEAM_RANK_BAD);
      if (h(R_ND_QUERY_RAN))
        nd_query_ran = true;

      if (v == 32) {
        width_for_32 = actual_vector_width;
        range1_for_32 = h(R_ND_RANGE_1);
      }
    }

    // ---------------------------------------------------------------------
    // Register pressure vs delivered vector width
    // ---------------------------------------------------------------------
    std::printf("\nTeamPolicy(vector_size=32) with a live per-thread array:\n");
    std::printf("  doubles/thread | work items/team | local_range(1) | sg size\n");

    int narrowed_at = 0;
    int measured = 0;
    {
      const int ts = 8;              // 8 x 32 = 256 work items per group
      const int ls = 64;
      int sizes[6] = {4, 8, 16, 32, 64, 128};
      host_result_view hs[6];
      bool ok[6] = {false, false, false, false, false, false};
      ok[0] = run_register_pressure<4>  (ls, ts, 32, hs[0]);
      ok[1] = run_register_pressure<8>  (ls, ts, 32, hs[1]);
      ok[2] = run_register_pressure<16> (ls, ts, 32, hs[2]);
      ok[3] = run_register_pressure<32> (ls, ts, 32, hs[3]);
      ok[4] = run_register_pressure<64> (ls, ts, 32, hs[4]);
      ok[5] = run_register_pressure<128>(ls, ts, 32, hs[5]);

      for (int i = 0; i < 6; ++i) {
        if (!ok[i]) { std::printf("  %14d | launch rejected\n", sizes[i]); continue; }
        const int per_team =
            hs[i](R_WORK_ITEMS) / (ls * (hs[i](R_TEAM_SIZE) ? hs[i](R_TEAM_SIZE) : 1));
        std::printf("  %14d | %15d | %14d | %7d%s\n", sizes[i], per_team,
                    hs[i](R_ND_RANGE_1), hs[i](R_SUB_GROUP_SIZE),
                    (per_team != 32) ? "   <-- NARROWED" : "");
        ++measured;
        if (per_team != 32 && narrowed_at == 0)
          narrowed_at = sizes[i];
      }
    }

    // ---------------------------------------------------------------------
    // Flat kernel probes
    // ---------------------------------------------------------------------
    std::printf("\nRangePolicy (flat) kernels:\n");
    {
      host_result_view h = probe_range(1024, 0);
      std::printf("  default chunk size : work items=%d, sub_group size=%d, "
                  "max sg local_id=%d, ballot ok=%s\n",
                  h(R_WORK_ITEMS), h(R_SUB_GROUP_SIZE), h(R_MAX_SG_LOCAL_ID),
                  h(R_BALLOT_RAN) ? "yes" : "NO");
#if SACADO_PROBE_RANGE_ND_ITEM
      std::printf("  default chunk size : nd_item query ran=%s, "
                  "local_range=(%d,%d), max local_id(1)=%d\n",
                  h(R_ND_QUERY_RAN) ? "yes" : "no", h(R_ND_RANGE_0),
                  h(R_ND_RANGE_1), h(R_MAX_LOCAL_ID_1));
#endif
    }
    {
      host_result_view h = probe_range(1024, 64);
      std::printf("  chunk_size = 64    : work items=%d, sub_group size=%d, "
                  "max sg local_id=%d, ballot ok=%s\n",
                  h(R_WORK_ITEMS), h(R_SUB_GROUP_SIZE), h(R_MAX_SG_LOCAL_ID),
                  h(R_BALLOT_RAN) ? "yes" : "NO");
#if SACADO_PROBE_RANGE_ND_ITEM
      std::printf("  chunk_size = 64    : nd_item query ran=%s, "
                  "local_range=(%d,%d), max local_id(1)=%d\n",
                  h(R_ND_QUERY_RAN) ? "yes" : "no", h(R_ND_RANGE_0),
                  h(R_ND_RANGE_1), h(R_MAX_LOCAL_ID_1));
#endif
    }

    // ---------------------------------------------------------------------
    // Verdict
    // ---------------------------------------------------------------------
    std::printf("\n---- findings ----\n");

    if (measured == 0)
      std::printf("0. Register pressure sweep did not run -- no team launch was "
                  "accepted.\n   Nothing measured about vector-width "
                  "narrowing.\n");
    else if (narrowed_at == 0)
      std::printf("0. Vector width stayed 32 at every register load tested.  "
                  "Kokkos did not narrow\n   the vector dimension for these "
                  "kernels on this device.\n");
    else {
      std::printf("0. Vector width DROPPED below 32 once a thread held %d "
                  "doubles.\n", narrowed_at);
      std::printf("   Sacado's partitioned Fad type is sized at compile time "
                  "from the layout\n   stride, so a kernel narrowed this way "
                  "writes more derivative components\n   than its local Fad "
                  "holds.  See SacadoSyclStrideIssue.txt.\n");
    }


    std::printf("1. Actual vector width for a requested 32: %d\n",
                width_for_32);
    if (width_for_32 != 32)
      std::printf("   Kokkos clamped it.  Sacado's LayoutContiguous stride for "
                  "SYCL hierarchical tests must be %d, not 32.\n",
                  width_for_32);

#if SACADO_PROBE_ND_ITEM
    if (!nd_query_ran) {
      std::printf("2. The nd_item free-function query compiled but never "
                  "reported a value -- treat as unusable.\n");
      success = false;
    } else {
      std::printf("2. nd_item query works in team kernels.  "
                  "local_range(1) for requested 32 = %d "
                  "(work-item count said %d -- %s).\n",
                  range1_for_32, width_for_32,
                  (range1_for_32 == width_for_32) ? "agree" : "DISAGREE");
      if (range1_for_32 != width_for_32)
        success = false;

      std::printf("3. local_id(0) == team_rank() : %s (%d mismatches)\n",
                  id0_bad_total == 0 ? "yes" : "NO", id0_bad_total);
      std::printf("4. ThreadVectorRange lane == local_id(1) : %s (%d "
                  "mismatches)\n",
                  lane_bad_total == 0 ? "yes" : "NO", lane_bad_total);
      if (id0_bad_total != 0 || lane_bad_total != 0)
        success = false;

      if (success)
        std::printf("\n   => get_nd_item<2>().get_local_id(1) is the vector "
                    "lane (threadIdx.x) and get_local_range(1) is the vector "
                    "width (blockDim.x).\n");
    }
#else
    std::printf("2-4. nd_item probes were disabled at compile time.\n");
#endif

    std::printf("\n");
  }
  Kokkos::finalize();

  // The probe is diagnostic: it fails only if an assumption the port depends
  // on is actually violated.  A clamped vector width is information, not a
  // failure.
  if (success)
    std::printf("End Result: TEST PASSED\n");
  else
    std::printf("End Result: TEST FAILED\n");

  return success ? 0 : 1;
}

#endif // KOKKOS_ENABLE_SYCL
