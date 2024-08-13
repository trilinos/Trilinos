#include "gtest/gtest.h"
#include "stk_search/morton_lbvh/MortonLBVH_ParallelConsistencyUtils.hpp"
#include "stk_search/Box.hpp"
#include "stk_search/BoxIdent.hpp"
#include "stk_search/IdentProc.hpp"
#include "stk_util/parallel/Parallel.hpp"

namespace {
using BoxType = stk::search::Box<double>;
using IdentProcType = stk::search::IdentProc<int, int>;
using BoxIdentProcType = stk::search::BoxIdentProc<BoxType, IdentProcType>;
using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using BoxIdentProcView = Kokkos::View<BoxIdentProcType*, ExecutionSpace>;

using HostExecutionSpace = Kokkos::DefaultHostExecutionSpace;
using BoxIdentProcViewHost = Kokkos::View<BoxIdentProcType*, HostExecutionSpace>;
}

TEST(ParallelConsistencyUtils, ProcBoundingBoxView1Proc)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 1)
  {
    GTEST_SKIP();
  }

  ExecutionSpace execSpace{};

  int myrank = stk::parallel_machine_rank(stk::parallel_machine_world());
  Kokkos::View<BoxIdentProcType*, ExecutionSpace> boxIdentProcs("box_ident_procs", 2);
  double delta = 0.5;
  double x0 = delta * myrank;
  boxIdentProcs(0) = BoxIdentProcType{BoxType(x0, 0,     0, x0 + delta,   delta, delta), IdentProcType(0, myrank)};
  boxIdentProcs(1) = BoxIdentProcType{BoxType(x0, delta, 0, x0 + delta, 2*delta, delta), IdentProcType(1, myrank)};

  Kokkos::View<BoxType*, ExecutionSpace> procBoxes = stk::search::gather_all_processor_superset_domain_boxes(boxIdentProcs, execSpace, stk::parallel_machine_world());

  auto procBoxesHost = Kokkos::create_mirror_view_and_copy(execSpace, procBoxes);

  EXPECT_EQ(procBoxesHost.extent(0), 1u);
  EXPECT_EQ(procBoxesHost(0), BoxType(0, 0, 0,   delta, 2*delta, delta));
}

TEST(ParallelConsistencyUtils, ProcBoundingBoxView1ProcSphere)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 1)
  {
    GTEST_SKIP();
  }

  using BoxType = stk::search::Sphere<double>;
  using BoxIdentProcType = stk::search::BoxIdentProc<BoxType, IdentProcType>;

  ExecutionSpace execSpace{};

  int myrank = stk::parallel_machine_rank(stk::parallel_machine_world());
  Kokkos::View<BoxIdentProcType*, ExecutionSpace> boxIdentProcs("box_ident_procs", 2);
  boxIdentProcs(0) = BoxIdentProcType{BoxType({0, 0, 0}, 1), IdentProcType(0, myrank)};
  boxIdentProcs(1) = BoxIdentProcType{BoxType({1, 0, 0}, 1), IdentProcType(1, myrank)};

  Kokkos::View<stk::search::Box<double>*, ExecutionSpace> procBoxes = stk::search::gather_all_processor_superset_domain_boxes(boxIdentProcs, execSpace, stk::parallel_machine_world());

  auto procBoxesHost = Kokkos::create_mirror_view_and_copy(execSpace, procBoxes);

  EXPECT_EQ(procBoxesHost.extent(0), 1u);
  EXPECT_EQ(procBoxesHost(0), stk::search::Box<double>(-1, -1, -1, 2, 1, 1));
}

TEST(ParallelConsistencyUtils, ProcBoundingBoxView2Proc)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2)
  {
    GTEST_SKIP();
  }

  ExecutionSpace execSpace{};

  int myrank = stk::parallel_machine_rank(stk::parallel_machine_world());
  Kokkos::View<BoxIdentProcType*, ExecutionSpace> boxIdentProcs("box_ident_procs", 2);
  double delta = 0.5;
  double x0 = delta * myrank;
  boxIdentProcs(0) = BoxIdentProcType{BoxType(x0, 0,     0, x0 + delta,   delta, delta), IdentProcType(0, myrank)};
  boxIdentProcs(1) = BoxIdentProcType{BoxType(x0, delta, 0, x0 + delta, 2*delta, delta), IdentProcType(1, myrank)};

  Kokkos::View<BoxType*, ExecutionSpace> procBoxes = stk::search::gather_all_processor_superset_domain_boxes(boxIdentProcs, execSpace, stk::parallel_machine_world());

  auto procBoxesHost = Kokkos::create_mirror_view_and_copy(execSpace, procBoxes);

  EXPECT_EQ(procBoxesHost.extent(0), 2u);
  EXPECT_EQ(procBoxesHost(0), BoxType(0,     0, 0,   delta, 2*delta, delta));
  EXPECT_EQ(procBoxesHost(1), BoxType(delta, 0, 0, 2*delta, 2*delta, delta));
}

TEST(ParallelConsistencyUtils, ProcBoundingBoxView4Proc)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 4)
  {
    GTEST_SKIP();
  }

  ExecutionSpace execSpace{};

  int myrank = stk::parallel_machine_rank(stk::parallel_machine_world());
  int i = myrank % 2;
  int j = myrank / 2;
  Kokkos::View<BoxIdentProcType*, ExecutionSpace> boxIdentProcs("box_ident_procs", 1);

  double delta = 0.5;
  double x0 = delta * i;
  double y0 = delta * j;
  boxIdentProcs(0) = BoxIdentProcType{BoxType(x0, y0, 0, x0 + delta, y0 + delta, delta), IdentProcType(0, myrank)};

  Kokkos::View<BoxType*, ExecutionSpace> procBoxes = stk::search::gather_all_processor_superset_domain_boxes(boxIdentProcs, execSpace, stk::parallel_machine_world());

  auto procBoxesHost = Kokkos::create_mirror_view_and_copy(execSpace, procBoxes);

  EXPECT_EQ(procBoxesHost.extent(0), 4u);
  EXPECT_EQ(procBoxesHost(0), BoxType(0,     0,     0, delta,   delta,   delta));
  EXPECT_EQ(procBoxesHost(1), BoxType(delta, 0,     0, 2*delta, delta,   delta));
  EXPECT_EQ(procBoxesHost(2), BoxType(0,     delta, 0, delta,   2*delta, delta));
  EXPECT_EQ(procBoxesHost(3), BoxType(delta, delta, 0, 2*delta, 2*delta, delta));
}

TEST(ParallelConsistencyUtils, ExtendRangeWithRemoteBoxesLocal)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 1)
  {
    GTEST_SKIP();
  }

  BoxIdentProcView domainBoxes("domain_boxes", 1);
  BoxIdentProcView rangeBoxes("range_boxes", 2);

  auto domainBoxesHost = Kokkos::create_mirror_view(Kokkos::DefaultHostExecutionSpace{}, domainBoxes);
  auto rangeBoxesHost = Kokkos::create_mirror_view(Kokkos::DefaultHostExecutionSpace{}, rangeBoxes);

  domainBoxesHost(0) = BoxIdentProcType{BoxType(0, 0, 0, 1, 1, 1), IdentProcType(0, 0)};

  rangeBoxesHost(0) = BoxIdentProcType{BoxType(0, 0, 0, 0.5, 1, 1), IdentProcType(0, 0)};
  rangeBoxesHost(1) = BoxIdentProcType{BoxType(0.5, 0, 0, 1, 1, 1), IdentProcType(1, 0)};

  Kokkos::deep_copy(domainBoxes, domainBoxesHost);
  Kokkos::deep_copy(rangeBoxes, rangeBoxesHost);

  auto [extendedRange, remoteIdents] = stk::search::morton_extend_local_range_with_remote_boxes_that_might_intersect(
                                        domainBoxes, rangeBoxes, ExecutionSpace{}, stk::parallel_machine_world());

  EXPECT_EQ(extendedRange.extent(0), 2U);
  EXPECT_EQ(remoteIdents.extent(0), 0U);
}

TEST(ParallelConsistencyUtils, ExtendRangeWithRemoteBoxes2Procs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2)
  {
    GTEST_SKIP();
  }

  BoxIdentProcView domainBoxes("domain_boxes", 1);
  BoxIdentProcView rangeBoxes("range_boxes", 2);

  auto domainBoxesHost = Kokkos::create_mirror_view(Kokkos::DefaultHostExecutionSpace{}, domainBoxes);
  auto rangeBoxesHost = Kokkos::create_mirror_view(Kokkos::DefaultHostExecutionSpace{}, rangeBoxes);

  int myrank =  stk::parallel_machine_rank(stk::parallel_machine_world());;
  double delta_x = 0.5;
  double x0 = delta_x * myrank;
  domainBoxesHost(0) = BoxIdentProcType{BoxType(x0, 0, 0, x0 + delta_x, 1, 1), IdentProcType(0, 0)};

  if (myrank == 0)
  {
    rangeBoxesHost(0) = BoxIdentProcType{BoxType(0,   0, 0, 0.4, 1, 1), IdentProcType(0, myrank)};
    rangeBoxesHost(1) = BoxIdentProcType{BoxType(0.4, 0, 0, 0.6, 1, 1), IdentProcType(1, myrank)};
  } else
  {
    rangeBoxesHost(0) = BoxIdentProcType{BoxType(0.6,   0, 0, 1, 1, 1), IdentProcType(0, myrank)};
    rangeBoxesHost(1) = BoxIdentProcType{BoxType(0.4, 0, 0, 0.7, 1, 1), IdentProcType(1, myrank)};
  }

  Kokkos::deep_copy(domainBoxes, domainBoxesHost);
  Kokkos::deep_copy(rangeBoxes, rangeBoxesHost);

  auto [extendedRange, remoteIdentProcs] = stk::search::morton_extend_local_range_with_remote_boxes_that_might_intersect(
                                        domainBoxes, rangeBoxes, ExecutionSpace{}, stk::parallel_machine_world());

  EXPECT_EQ(extendedRange.extent(0), 3U);
  EXPECT_EQ(remoteIdentProcs.extent(0), 1U);

  auto extendedRangeHost    = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace{}, extendedRange);
  auto remoteIdentProcsHost = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace{}, remoteIdentProcs);

  if (myrank == 0)
  {
    EXPECT_EQ(extendedRangeHost(0), rangeBoxesHost(0).box);

    EXPECT_EQ(extendedRangeHost(1), rangeBoxesHost(1).box);

    EXPECT_EQ(extendedRangeHost(2),    BoxType(0.4, 0, 0, 0.7, 1, 1));
    EXPECT_EQ(remoteIdentProcsHost(0), IdentProcType(1, 1));
  } else
  {
    EXPECT_EQ(extendedRangeHost(0), rangeBoxesHost(0).box);

    EXPECT_EQ(extendedRangeHost(1), rangeBoxesHost(1).box);

    EXPECT_EQ(extendedRangeHost(2), BoxType(0.4, 0, 0, 0.6, 1, 1));
    EXPECT_EQ(remoteIdentProcsHost(0), IdentProcType(1, 0));
  }
}

TEST(ParallelConsistencyUtils, ExtendRangeWithRemoteBoxes4Procs)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 4)
  {
    GTEST_SKIP();
  }

  BoxIdentProcView domainBoxes("domain_boxes", 1);
  BoxIdentProcView rangeBoxes("range_boxes", 2);

  auto domainBoxesHost = Kokkos::create_mirror_view(Kokkos::DefaultHostExecutionSpace{}, domainBoxes);
  auto rangeBoxesHost = Kokkos::create_mirror_view(Kokkos::DefaultHostExecutionSpace{}, rangeBoxes);

  int myrank =  stk::parallel_machine_rank(stk::parallel_machine_world());;

  if (myrank == 0)
  {
    domainBoxesHost(0) = BoxIdentProcType{BoxType(0, 0, 0, 0.5, 0.5, 1), IdentProcType(0, myrank)};

    rangeBoxesHost(0) = BoxIdentProcType{BoxType(0,   0,   0, 0.4, 0.4, 1), IdentProcType(0, myrank)};
    rangeBoxesHost(1) = BoxIdentProcType{BoxType(0.3, 0.4, 0, 0.6, 0.6, 1), IdentProcType(1, myrank)};
  } else if (myrank == 1)
  {
    domainBoxesHost(0) = BoxIdentProcType{BoxType(0.5, 0, 0, 1.0, 0.5, 1), IdentProcType(0, myrank)};

    rangeBoxesHost(0) = BoxIdentProcType{BoxType(0.6, 0.0, 0, 1.0, 0.4, 1), IdentProcType(0, myrank)};
    rangeBoxesHost(1) = BoxIdentProcType{BoxType(0.2, 0.4, 0, 0.6, 0.6, 1), IdentProcType(1, myrank)};
  } else if (myrank == 2)
  {
    domainBoxesHost(0) = BoxIdentProcType{BoxType(0.0, 0.5, 0, 0.5, 1.0, 1), IdentProcType(0, myrank)};

    rangeBoxesHost(0) = BoxIdentProcType{BoxType(0.0, 0.6, 0, 0.4, 1.0, 1), IdentProcType(0, myrank)};
    rangeBoxesHost(1) = BoxIdentProcType{BoxType(0.4, 0.4, 0, 0.7, 0.6, 1), IdentProcType(1, myrank)};
  } else
  {
    domainBoxesHost(0) = BoxIdentProcType{BoxType(0.5, 0.5, 0, 1.0, 1.0, 1), IdentProcType(0, myrank)};

    rangeBoxesHost(0) = BoxIdentProcType{BoxType(0.6, 0.6, 0, 1.0, 1.0, 1), IdentProcType(0, myrank)};
    rangeBoxesHost(1) = BoxIdentProcType{BoxType(0.4, 0.4, 0, 0.8, 0.6, 1), IdentProcType(1, myrank)};
  }

  Kokkos::deep_copy(domainBoxes, domainBoxesHost);
  Kokkos::deep_copy(rangeBoxes, rangeBoxesHost);

  auto [extendedRange, remoteIdentProcs] = stk::search::morton_extend_local_range_with_remote_boxes_that_might_intersect(
                                        domainBoxes, rangeBoxes, ExecutionSpace{}, stk::parallel_machine_world());

  EXPECT_EQ(extendedRange.extent(0), 5U);
  EXPECT_EQ(remoteIdentProcs.extent(0), 3U);

  auto extendedRangeHost    = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace{}, extendedRange);
  auto remoteIdentProcsHost = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace{}, remoteIdentProcs);

  EXPECT_EQ(extendedRangeHost(0), rangeBoxesHost(0).box);
  EXPECT_EQ(extendedRangeHost(1), rangeBoxesHost(1).box);
  if (myrank == 0)
  {
    EXPECT_EQ(extendedRangeHost(2),    BoxType(0.2, 0.4, 0, 0.6, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(0), IdentProcType(1, 1));

    EXPECT_EQ(extendedRangeHost(3),    BoxType(0.4, 0.4, 0, 0.7, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(1), IdentProcType(1, 2));

    EXPECT_EQ(extendedRangeHost(4),    BoxType(0.4, 0.4, 0, 0.8, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(2), IdentProcType(1, 3));
  } else if (myrank == 1)
  {
    EXPECT_EQ(extendedRangeHost(2),    BoxType(0.3, 0.4, 0, 0.6, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(0), IdentProcType(1, 0));

    EXPECT_EQ(extendedRangeHost(3),    BoxType(0.4, 0.4, 0, 0.7, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(1), IdentProcType(1, 2));

    EXPECT_EQ(extendedRangeHost(4),    BoxType(0.4, 0.4, 0, 0.8, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(2), IdentProcType(1, 3));
  } else if (myrank == 2)
  {
    EXPECT_EQ(extendedRangeHost(2),    BoxType(0.3, 0.4, 0, 0.6, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(0), IdentProcType(1, 0));

    EXPECT_EQ(extendedRangeHost(3),    BoxType(0.2, 0.4, 0, 0.6, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(1), IdentProcType(1, 1));

    EXPECT_EQ(extendedRangeHost(4),    BoxType(0.4, 0.4, 0, 0.8, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(2), IdentProcType(1, 3));
  } else
  {
    EXPECT_EQ(extendedRangeHost(2),    BoxType(0.3, 0.4, 0, 0.6, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(0), IdentProcType(1, 0));

    EXPECT_EQ(extendedRangeHost(3),    BoxType(0.2, 0.4, 0, 0.6, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(1), IdentProcType(1, 1));

    EXPECT_EQ(extendedRangeHost(4),    BoxType(0.4, 0.4, 0, 0.7, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(2), IdentProcType(1, 2));
  }
}

TEST(ParallelConsistencyUtils, ExtendRangeWithRemoteBoxes4ProcsHostSpace)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 4)
  {
    GTEST_SKIP();
  }

  BoxIdentProcViewHost domainBoxesHost("domain_boxes", 1);
  BoxIdentProcViewHost rangeBoxesHost("range_boxes", 2);

  int myrank =  stk::parallel_machine_rank(stk::parallel_machine_world());;

  if (myrank == 0)
  {
    domainBoxesHost(0) = BoxIdentProcType{BoxType(0, 0, 0, 0.5, 0.5, 1), IdentProcType(0, myrank)};

    rangeBoxesHost(0) = BoxIdentProcType{BoxType(0,   0,   0, 0.4, 0.4, 1), IdentProcType(0, myrank)};
    rangeBoxesHost(1) = BoxIdentProcType{BoxType(0.3, 0.4, 0, 0.6, 0.6, 1), IdentProcType(1, myrank)};
  } else if (myrank == 1)
  {
    domainBoxesHost(0) = BoxIdentProcType{BoxType(0.5, 0, 0, 1.0, 0.5, 1), IdentProcType(0, myrank)};

    rangeBoxesHost(0) = BoxIdentProcType{BoxType(0.6, 0.0, 0, 1.0, 0.4, 1), IdentProcType(0, myrank)};
    rangeBoxesHost(1) = BoxIdentProcType{BoxType(0.2, 0.4, 0, 0.6, 0.6, 1), IdentProcType(1, myrank)};
  } else if (myrank == 2)
  {
    domainBoxesHost(0) = BoxIdentProcType{BoxType(0.0, 0.5, 0, 0.5, 1.0, 1), IdentProcType(0, myrank)};

    rangeBoxesHost(0) = BoxIdentProcType{BoxType(0.0, 0.6, 0, 0.4, 1.0, 1), IdentProcType(0, myrank)};
    rangeBoxesHost(1) = BoxIdentProcType{BoxType(0.4, 0.4, 0, 0.7, 0.6, 1), IdentProcType(1, myrank)};
  } else
  {
    domainBoxesHost(0) = BoxIdentProcType{BoxType(0.5, 0.5, 0, 1.0, 1.0, 1), IdentProcType(0, myrank)};

    rangeBoxesHost(0) = BoxIdentProcType{BoxType(0.6, 0.6, 0, 1.0, 1.0, 1), IdentProcType(0, myrank)};
    rangeBoxesHost(1) = BoxIdentProcType{BoxType(0.4, 0.4, 0, 0.8, 0.6, 1), IdentProcType(1, myrank)};
  }

  auto [extendedRangeHost, remoteIdentProcsHost] = stk::search::morton_extend_local_range_with_remote_boxes_that_might_intersect(
                                        domainBoxesHost, rangeBoxesHost, HostExecutionSpace{}, stk::parallel_machine_world());

  EXPECT_EQ(extendedRangeHost.extent(0), 5U);
  EXPECT_EQ(remoteIdentProcsHost.extent(0), 3U);

  EXPECT_EQ(extendedRangeHost(0), rangeBoxesHost(0).box);
  EXPECT_EQ(extendedRangeHost(1), rangeBoxesHost(1).box);
  if (myrank == 0)
  {
    EXPECT_EQ(extendedRangeHost(2),    BoxType(0.2, 0.4, 0, 0.6, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(0), IdentProcType(1, 1));

    EXPECT_EQ(extendedRangeHost(3),    BoxType(0.4, 0.4, 0, 0.7, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(1), IdentProcType(1, 2));

    EXPECT_EQ(extendedRangeHost(4),    BoxType(0.4, 0.4, 0, 0.8, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(2), IdentProcType(1, 3));
  } else if (myrank == 1)
  {
    EXPECT_EQ(extendedRangeHost(2),    BoxType(0.3, 0.4, 0, 0.6, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(0), IdentProcType(1, 0));

    EXPECT_EQ(extendedRangeHost(3),    BoxType(0.4, 0.4, 0, 0.7, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(1), IdentProcType(1, 2));

    EXPECT_EQ(extendedRangeHost(4),    BoxType(0.4, 0.4, 0, 0.8, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(2), IdentProcType(1, 3));
  } else if (myrank == 2)
  {
    EXPECT_EQ(extendedRangeHost(2),    BoxType(0.3, 0.4, 0, 0.6, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(0), IdentProcType(1, 0));

    EXPECT_EQ(extendedRangeHost(3),    BoxType(0.2, 0.4, 0, 0.6, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(1), IdentProcType(1, 1));

    EXPECT_EQ(extendedRangeHost(4),    BoxType(0.4, 0.4, 0, 0.8, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(2), IdentProcType(1, 3));
  } else
  {
    EXPECT_EQ(extendedRangeHost(2),    BoxType(0.3, 0.4, 0, 0.6, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(0), IdentProcType(1, 0));

    EXPECT_EQ(extendedRangeHost(3),    BoxType(0.2, 0.4, 0, 0.6, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(1), IdentProcType(1, 1));

    EXPECT_EQ(extendedRangeHost(4),    BoxType(0.4, 0.4, 0, 0.7, 0.6, 1));
    EXPECT_EQ(remoteIdentProcsHost(2), IdentProcType(1, 2));
  }
}
