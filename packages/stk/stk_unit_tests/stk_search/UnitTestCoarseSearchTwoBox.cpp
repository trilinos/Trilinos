// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/LocalCoarseSearch.hpp>
#include <stk_unit_test_utils/Search_UnitTestUtils.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/ngp/NgpSpaces.hpp>

#include <gtest/gtest.h>
#include <vector>
#include "stk_search/SearchMethod.hpp"

namespace
{

void runTwoBoxTest(stk::search::SearchMethod searchMethod, const double distanceBetweenBoxCenters,
                   const double boxSize, const unsigned expectedNumOverlap)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(comm);
  int procId = stk::parallel_machine_rank(comm);

  std::vector<std::pair<StkBox, IdentProc>> boxVector1;
  if (procId == 0) {
    boxVector1.push_back(stk::unit_test_util::generateBoundingVolume<StkBox>(0, 0, 0, boxSize/2, 1, procId));
  }

  std::vector<std::pair<StkBox, IdentProc>> boxVector2;
  if (procId == numProcs-1) {
    boxVector2.push_back(stk::unit_test_util::generateBoundingVolume<StkBox>(distanceBetweenBoxCenters, 0, 0, boxSize/2, 2, procId));
  }

  SearchResults boxIdPairResults;
  stk::search::coarse_search(boxVector1, boxVector2, searchMethod, comm, boxIdPairResults);

  if (procId == 0 || (procId == numProcs-1)) {
    EXPECT_EQ(expectedNumOverlap, boxIdPairResults.size());
  }
  else {
    EXPECT_EQ(0u, boxIdPairResults.size());
  }
}

void device_runTwoBoxTest(stk::search::SearchMethod searchMethod, const double distanceBetweenBoxCenters,
                                const double boxSize, const unsigned expectedNumOverlap)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(comm);
  int procId = stk::parallel_machine_rank(comm);

  auto domain = Kokkos::View<StkBoxIdentProc*, stk::ngp::ExecSpace>("domain box-ident", 1);
  auto range = Kokkos::View<StkBoxIdentProc*, stk::ngp::ExecSpace>("range box-ident", 1);

  if (procId == 0) {
    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
      KOKKOS_LAMBDA(const unsigned & i) {
        domain[0] =
            stk::unit_test_util::device_generateBoxIdentProc<StkBox, IdentProc>(0, 0, 0, boxSize/2, 1, procId);
      });
  }

  if (procId == numProcs - 1) {
    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
      KOKKOS_LAMBDA(const unsigned & i) {
        range[0] =
            stk::unit_test_util::device_generateBoxIdentProc<StkBox, IdentProc>(distanceBetweenBoxCenters, 0, 0,
                                                                                          boxSize/2, 2, procId);
      });
  }

  auto intersections = Kokkos::View<IdentProcIntersection*, stk::ngp::ExecSpace>("intersections", 0);

  stk::search::coarse_search(domain, range, searchMethod, comm, intersections);

  Kokkos::View<IdentProcIntersection*>::HostMirror hostIntersections = Kokkos::create_mirror_view(intersections);
  Kokkos::deep_copy(hostIntersections, intersections);

  if (procId == 0 || (procId == numProcs-1)) {
    EXPECT_EQ(expectedNumOverlap, intersections.extent(0));
  }
  else {
    EXPECT_EQ(0u, intersections.extent(0));
  }
}

constexpr double boxSize = 1.0;

TEST(CoarseSearchCorrectness, OverlappingBoxes_KDTREE)
{
  const double distanceBetweenBoxCenters = 0.5;
  const unsigned expectedNumOverlap = 1;
  runTwoBoxTest(stk::search::KDTREE, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, NonOverlappingBoxes_KDTREE)
{
  const double distanceBetweenBoxCenters = 2.0;
  const unsigned expectedNumOverlap = 0;
  runTwoBoxTest(stk::search::KDTREE, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, JustEdgeOverlappingBoxes_KDTREE)
{
  double distanceBetweenBoxCenters = 0.999999999;
  const unsigned expectedNumOverlap = 1;
  runTwoBoxTest(stk::search::KDTREE, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, NotQuiteEdgeOverlappingBoxes_KDTREE)
{
  double distanceBetweenBoxCenters = 1.0000000001;
  const unsigned expectedNumOverlap = 0;
  runTwoBoxTest(stk::search::KDTREE, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}


TEST(CoarseSearchCorrectness, OverlappingBoxes_MORTON_LBVH)
{
  const double distanceBetweenBoxCenters = 0.5;
  const unsigned expectedNumOverlap = 1;
  runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, NonOverlappingBoxes_MORTON_LBVH)
{
  const double distanceBetweenBoxCenters = 2.0;
  const unsigned expectedNumOverlap = 0;
  runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, JustEdgeOverlappingBoxes_MORTON_LBVH)
{
  double distanceBetweenBoxCenters = 0.999999999;
  const unsigned expectedNumOverlap = 1;
  runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, NotQuiteEdgeOverlappingBoxes_MORTON_LBVH)
{
  double distanceBetweenBoxCenters = 1.0000000001;
  const unsigned expectedNumOverlap = 0;
  runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, OverlappingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double distanceBetweenBoxCenters = 0.5;
  const unsigned expectedNumOverlap = 1;
  runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, NonOverlappingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double distanceBetweenBoxCenters = 2.0;
  const unsigned expectedNumOverlap = 0;
  runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, JustEdgeOverlappingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  double distanceBetweenBoxCenters = 0.999999999;
  const unsigned expectedNumOverlap = 1;
  runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, NotQuiteEdgeOverlappingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  double distanceBetweenBoxCenters = 1.00001;
  const unsigned expectedNumOverlap = 0;
  runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, NotQuiteEdgeOverlappingBoxes_FloatTruncation_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  double distanceBetweenBoxCenters = 1.0000000001;
  const unsigned expectedNumOverlap = 1;
  runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

void host_local_runTwoBoxTest(stk::search::SearchMethod searchMethod,
    const double distanceBetweenBoxCenters,
    const double boxSize,
    const unsigned expectedNumOverlap)
{
  StkBoxIdentVector domain;
  StkBoxIdentVector range;

  domain.emplace_back(stk::unit_test_util::box_ident_to_pair(
      stk::unit_test_util::device_generateBoxIdent<StkBox, Ident>(0, 0, 0, boxSize / 2, 1)));

  range.emplace_back(stk::unit_test_util::box_ident_to_pair(
      stk::unit_test_util::device_generateBoxIdent<StkBox, Ident>(
          distanceBetweenBoxCenters, 0, 0, boxSize / 2, 2)));

  LocalSearchResults intersections;
  
  bool sortSearchResults = true;
  stk::search::local_coarse_search(domain, range, searchMethod, intersections, sortSearchResults);

  ASSERT_EQ(intersections.size(), expectedNumOverlap);

  for (unsigned i = 0; i < expectedNumOverlap; ++i) {
    EXPECT_EQ(intersections[i].first, 1);
    EXPECT_EQ(intersections[i].second, 2);
  }
}

void device_local_runTwoBoxTest(stk::search::SearchMethod searchMethod, const double distanceBetweenBoxCenters,
                                const double boxSize, const unsigned expectedNumOverlap)
{
  auto domain = Kokkos::View<StkBoxIdent*, stk::ngp::ExecSpace>("domain box-ident", 1);
  auto range = Kokkos::View<StkBoxIdent*, stk::ngp::ExecSpace>("range box-ident", 1);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
    KOKKOS_LAMBDA(const unsigned & i) {
      domain[0] =
          stk::unit_test_util::device_generateBoxIdent<StkBox, Ident>(0, 0, 0, boxSize/2, 1);

      range[0] =
          stk::unit_test_util::device_generateBoxIdent<StkBox, Ident>(distanceBetweenBoxCenters, 0, 0,
                                                                                          boxSize/2, 2);
    });

  auto intersections = Kokkos::View<IdentIntersection*, stk::ngp::ExecSpace>("intersections", 0);

  auto execSpace = stk::ngp::ExecSpace{};
  bool sortSearchResults = true;
  stk::search::local_coarse_search(domain, range, searchMethod, intersections, execSpace, sortSearchResults);

  Kokkos::View<IdentIntersection*>::HostMirror hostIntersections = Kokkos::create_mirror_view(intersections);
  Kokkos::deep_copy(hostIntersections, intersections);

  ASSERT_EQ(intersections.extent(0), expectedNumOverlap);

  for (unsigned i = 0; i < expectedNumOverlap; ++i) {
    EXPECT_EQ(intersections(i).domainIdent, 1);
    EXPECT_EQ(intersections(i).rangeIdent,  2);
  }
}


TEST(CoarseSearchCorrectness, Local_OverlappingBoxes_KDTREE)
{
  const double distanceBetweenBoxCenters = 0.5;
  const unsigned expectedNumOverlap = 1;
  host_local_runTwoBoxTest(stk::search::KDTREE, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Local_NonOverlappingBoxes_KDTREE)
{
  const double distanceBetweenBoxCenters = 2.0;
  const unsigned expectedNumOverlap = 0;
  host_local_runTwoBoxTest(stk::search::KDTREE, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Local_JustEdgeOverlappingBoxes_KDTREE)
{
  double distanceBetweenBoxCenters = 0.999999999;
  const unsigned expectedNumOverlap = 1;
  host_local_runTwoBoxTest(stk::search::KDTREE, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Local_NotQuiteEdgeOverlappingBoxes_KDTREE)
{
  double distanceBetweenBoxCenters = 1.0000000001;
  const unsigned expectedNumOverlap = 0;
  host_local_runTwoBoxTest(stk::search::KDTREE, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Ngp_Local_OverlappingBoxes_MORTON_LBVH)
{
  const double distanceBetweenBoxCenters = 0.5;
  const unsigned expectedNumOverlap = 1;
  host_local_runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_local_runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Ngp_Local_NonOverlappingBoxes_MORTON_LBVH)
{
  const double distanceBetweenBoxCenters = 2.0;
  const unsigned expectedNumOverlap = 0;
  host_local_runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_local_runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Ngp_Local_JustEdgeOverlappingBoxes_MORTON_LBVH)
{
  double distanceBetweenBoxCenters = 0.999999999;
  const unsigned expectedNumOverlap = 1;
  host_local_runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_local_runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Ngp_Local_NotQuiteEdgeOverlappingBoxes_MORTON_LBVH)
{
  double distanceBetweenBoxCenters = 1.0000000001;
  const unsigned expectedNumOverlap = 0;
  host_local_runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_local_runTwoBoxTest(stk::search::MORTON_LBVH, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Ngp_Local_OverlappingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double distanceBetweenBoxCenters = 0.5;
  const unsigned expectedNumOverlap = 1;
  host_local_runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_local_runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Ngp_Local_NonOverlappingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double distanceBetweenBoxCenters = 2.0;
  const unsigned expectedNumOverlap = 0;
  host_local_runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_local_runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Ngp_Local_JustEdgeOverlappingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  double distanceBetweenBoxCenters = 0.999999999;
  const unsigned expectedNumOverlap = 1;
  host_local_runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_local_runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Ngp_Local_NotQuiteEdgeOverlappingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  double distanceBetweenBoxCenters = 1.0000000001;
  const unsigned expectedNumOverlap = 1;
  host_local_runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
  device_local_runTwoBoxTest(stk::search::ARBORX, distanceBetweenBoxCenters, boxSize, expectedNumOverlap);
}

}
