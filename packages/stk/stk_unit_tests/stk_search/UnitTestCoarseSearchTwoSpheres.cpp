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

void runTwoSpheresTest(stk::search::SearchMethod searchMethod, const double distanceBetweenSphereCenters,
                       const double radius, const unsigned expectedNumOverlap)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(comm);
  int procId = stk::parallel_machine_rank(comm);

  std::vector<std::pair<Sphere, IdentProc>> boxVector1;
  if (procId == 0) {
    boxVector1.push_back(stk::unit_test_util::generateBoundingVolume<Sphere>(0, 0, 0, radius, 1, procId));
  }

  std::vector<std::pair<Sphere, IdentProc>> boxVector2;
  if (procId == numProcs-1) {
    boxVector2.push_back(stk::unit_test_util::generateBoundingVolume<Sphere>(distanceBetweenSphereCenters, 0, 0, radius, 2, procId));
  }

  SearchResults boxIdPairResults;
  stk::search::coarse_search(boxVector1, boxVector2, searchMethod, comm, boxIdPairResults);

  if (procId == 0 || (procId == numProcs-1)) {
    ASSERT_EQ(expectedNumOverlap, boxIdPairResults.size());
  }
  else {
    ASSERT_EQ(0u, boxIdPairResults.size());
  }
}

void device_runTwoSpheresTest(stk::search::SearchMethod searchMethod, const double distanceBetweenSphereCenters,
                                    const double radius, const unsigned expectedNumOverlap)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(comm);
  int procId = stk::parallel_machine_rank(comm);

  auto domain = Kokkos::View<SphereIdentProc*, stk::ngp::ExecSpace>("domain box-ident-proc", 1);
  auto range = Kokkos::View<SphereIdentProc*, stk::ngp::ExecSpace>("range box-ident-proc", 1);

  if (procId == 0) {
    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
      KOKKOS_LAMBDA(const unsigned & i) {
        domain[0] =
            stk::unit_test_util::device_generateBoxIdentProc<Sphere, IdentProc>(0, 0, 0, radius, 1, procId);
    });
  }

  if (procId == numProcs - 1) {
    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
      KOKKOS_LAMBDA(const unsigned & i) {
        const double axisOffset = distanceBetweenSphereCenters / sqrt(2.0);
        range[0] =
            stk::unit_test_util::device_generateBoxIdentProc<Sphere, IdentProc>(axisOffset, axisOffset, 0,
                                                                                                    radius, 2, procId);
    });
  }

  auto intersections = Kokkos::View<IdentProcIntersection*, stk::ngp::ExecSpace>("intersections", 0);

  stk::search::coarse_search(domain, range, searchMethod, comm, intersections);

  Kokkos::View<IdentProcIntersection*>::HostMirror hostIntersections = Kokkos::create_mirror_view(intersections);
  Kokkos::deep_copy(hostIntersections, intersections);

  if (procId == 0 || (procId == numProcs-1)) {
    ASSERT_EQ(expectedNumOverlap, intersections.size());
  }
  else {
    ASSERT_EQ(0u, intersections.size());
  }

}

const double radiusOfOneHalf = 0.5;

TEST(CoarseSearchCorrectness, OverlappingSpheres_KDTREE)
{
  double distanceBetweenSphereCenters = 0.5;
  const unsigned expectedNumOverlap = 1;
  runTwoSpheresTest(stk::search::KDTREE, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, NonOverlappingSpheres_KDTREE)
{
  double distanceBetweenSphereCenters = 2.0;
  const unsigned expectedNumOverlap = 0;
  runTwoSpheresTest(stk::search::KDTREE, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, JustEdgeOverlappingSpheres_KDTREE)
{
  double distanceBetweenSphereCenters = 0.999999999;
  const unsigned expectedNumOverlap = 1;
  runTwoSpheresTest(stk::search::KDTREE, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, NotQuiteEdgeOverlappingSpheres_KDTREE)
{
  double distanceBetweenSphereCenters = 1.0000000001;
  const unsigned expectedNumOverlap = 0;
  runTwoSpheresTest(stk::search::KDTREE, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, OverlappingSpheres_MORTON_LBVH)
{
  double distanceBetweenSphereCenters = 0.5;
  const unsigned expectedNumOverlap = 1;
  runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, NonOverlappingSpheres_MORTON_LBVH)
{
  double distanceBetweenSphereCenters = 2.0;
  const unsigned expectedNumOverlap = 0;
  runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, JustEdgeOverlappingSpheres_MORTON_LBVH)
{
  double distanceBetweenSphereCenters = 0.999999999;
  const unsigned expectedNumOverlap = 1;
  runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, NotQuiteEdgeOverlappingSpheres_MORTON_LBVH)
{
  double distanceBetweenSphereCenters = 1.0000000001;
  const unsigned expectedNumOverlap = 0;
  runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, OverlappingSpheres_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  double distanceBetweenSphereCenters = 0.5;
  const unsigned expectedNumOverlap = 1;
  runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, NonOverlappingSpheres_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  double distanceBetweenSphereCenters = 2.0;
  const unsigned expectedNumOverlap = 0;
  runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, JustEdgeOverlappingSpheres_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  double distanceBetweenSphereCenters = 0.999999999;
  const unsigned expectedNumOverlap = 1;
  runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, NotQuiteEdgeOverlappingSpheres_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  double distanceBetweenSphereCenters = 1.00001;
  const unsigned expectedNumOverlap = 0;
  runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, NotQuiteEdgeOverlappingSpheres_FloatTruncation_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  double distanceBetweenSphereCenters = 1.0000000001;
  const unsigned expectedNumOverlap = 0;
  runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

void host_local_runTwoSpheresTest(stk::search::SearchMethod searchMethod,
    const double distanceBetweenSphereCenters,
    const double radius,
    const unsigned expectedNumOverlap)
{
  std::vector<SphereIdentPair> domain;
  std::vector<SphereIdentPair> range;

  domain.emplace_back(stk::unit_test_util::box_ident_to_pair(
      stk::unit_test_util::device_generateBoxIdent<Sphere, Ident>(0, 0, 0, radius, 1)));

  const double axisOffset = distanceBetweenSphereCenters / sqrt(2.0);
  range.emplace_back(stk::unit_test_util::box_ident_to_pair(
      stk::unit_test_util::device_generateBoxIdent<Sphere, Ident>(
          axisOffset, axisOffset, 0, radius, 2)));

  LocalSearchResults intersections;

  bool sortSearchResults = true;
  stk::search::local_coarse_search(domain, range, searchMethod, intersections, sortSearchResults);

  ASSERT_EQ(intersections.size(), expectedNumOverlap);

  for (unsigned i = 0; i < expectedNumOverlap; ++i) {
    EXPECT_EQ(intersections[i].first, 1);
    EXPECT_EQ(intersections[i].second, 2);
  }
}

void device_local_runTwoSpheresTest(stk::search::SearchMethod searchMethod, const double distanceBetweenSphereCenters,
                                    const double radius, const unsigned expectedNumOverlap)
{
  auto domain = Kokkos::View<SphereIdent*, stk::ngp::ExecSpace>("domain box-ident", 1);
  auto range = Kokkos::View<SphereIdent*, stk::ngp::ExecSpace>("range box-ident", 1);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
    KOKKOS_LAMBDA(const unsigned & i) {
      domain[0] =
          stk::unit_test_util::device_generateBoxIdent<Sphere, Ident>(0, 0, 0, radius, 1);

      const double axisOffset = distanceBetweenSphereCenters / sqrt(2.0);
      range[0] =
          stk::unit_test_util::device_generateBoxIdent<Sphere, Ident>(axisOffset, axisOffset, 0,
                                                                                          radius, 2);
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

TEST(CoarseSearchCorrectness, Local_OverlappingSpheres_KDTREE)
{
  double distanceBetweenSphereCenters = 0.5;
  const unsigned expectedNumOverlap = 1;
  host_local_runTwoSpheresTest(stk::search::KDTREE, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Local_NonOverlappingSpheres_KDTREE)
{
  double distanceBetweenSphereCenters = 2.0;
  const unsigned expectedNumOverlap = 0;
  host_local_runTwoSpheresTest(stk::search::KDTREE, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Local_JustEdgeOverlappingSpheres_KDTREE)
{
  double distanceBetweenSphereCenters = 0.999999999;
  const unsigned expectedNumOverlap = 1;
  host_local_runTwoSpheresTest(stk::search::KDTREE, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Local_NotQuiteEdgeOverlappingSpheres_KDTREE)
{
  double distanceBetweenSphereCenters = 1.0000000001;
  const unsigned expectedNumOverlap = 0;
  host_local_runTwoSpheresTest(stk::search::KDTREE, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}


TEST(CoarseSearchCorrectness, Ngp_Local_OverlappingSpheres_MORTON_LBVH)
{
  double distanceBetweenSphereCenters = 0.5;
  const unsigned expectedNumOverlap = 1;
  host_local_runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_local_runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Ngp_Local_NonOverlappingSpheres_MORTON_LBVH)
{
  double distanceBetweenSphereCenters = 2.0;
  const unsigned expectedNumOverlap = 0;
  host_local_runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_local_runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Ngp_Local_JustEdgeOverlappingSpheres_MORTON_LBVH)
{
  double distanceBetweenSphereCenters = 0.999999999;
  const unsigned expectedNumOverlap = 1;
  host_local_runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_local_runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Ngp_Local_NotQuiteEdgeOverlappingSpheres_MORTON_LBVH)
{
  double distanceBetweenSphereCenters = 1.0000000001;
  const unsigned expectedNumOverlap = 0;
  host_local_runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_local_runTwoSpheresTest(stk::search::MORTON_LBVH, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Ngp_Local_OverlappingSpheres_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  double distanceBetweenSphereCenters = 0.5;
  const unsigned expectedNumOverlap = 1;
  host_local_runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_local_runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Ngp_Local_NonOverlappingSpheres_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  double distanceBetweenSphereCenters = 2.0;
  const unsigned expectedNumOverlap = 0;
  host_local_runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_local_runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Ngp_Local_JustEdgeOverlappingSpheres_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  double distanceBetweenSphereCenters = 0.999999999;
  const unsigned expectedNumOverlap = 1;
  host_local_runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_local_runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

TEST(CoarseSearchCorrectness, Ngp_Local_NotQuiteEdgeOverlappingSpheres_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  double distanceBetweenSphereCenters = 1.0000000001;
  const unsigned expectedNumOverlap = 0;
  host_local_runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
  device_local_runTwoSpheresTest(stk::search::ARBORX, distanceBetweenSphereCenters, radiusOfOneHalf, expectedNumOverlap);
}

stk::search::Box<float> make_naive_arbor_x_box_from_double_sphere(stk::search::Sphere<double> sphere)
{

    auto center = sphere.center();
    float radius = static_cast<float>(sphere.radius());
    float x_min = static_cast<float>(center[0]) - radius;
    float y_min = static_cast<float>(center[1]) - radius;
    float z_min = static_cast<float>(center[2]) - radius;
    float x_max = static_cast<float>(center[0]) + radius;
    float y_max = static_cast<float>(center[1]) + radius;
    float z_max = static_cast<float>(center[2]) + radius;

    return stk::search::Box(x_min, y_min, z_min, x_max, y_max, z_max);
}

TEST(CoarseSearchCorrectness, ARBORX_FLOAT_PRECISION_SPHERES)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif

  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  typedef double SpherePrecision;

  SpherePrecision domainRadius = 1.400000000001;
  stk::search::Point<SpherePrecision> domainCenter(0.0, 0.0, 0.0);
  SpherePrecision rangeRadius = 1.0;
  stk::search::Point<SpherePrecision> rangeCenter(2.400000000001, 0.0, 0.0);

  stk::search::Sphere<SpherePrecision> domain(domainCenter, domainRadius);
  std::vector<std::pair<stk::search::Sphere<SpherePrecision>, IdentProc>> domainVector;
  domainVector.push_back({domain, IdentProc(0,0)});

  stk::search::Sphere<SpherePrecision> range(rangeCenter, rangeRadius);
  std::vector<std::pair<stk::search::Sphere<SpherePrecision>, IdentProc>> rangeVector;
  rangeVector.push_back({range, IdentProc(0,0)});

  std::vector<std::pair<IdentProc, IdentProc>> searchResultsMorton;
  stk::search::coarse_search(domainVector, rangeVector, stk::search::MORTON_LBVH, MPI_COMM_WORLD, searchResultsMorton);
  EXPECT_EQ(1u, searchResultsMorton.size());

  /* ArborX only uses single precision values; as a result, any double-precision StkBox will be cast down
     to floats. To prevent the case where an intersection between two boxes is missed due to truncation
     by ArborX, the boxes are expanded outward when converting from StkBox to ArborX. The two boxes
     constructed below are representative of the "naive" StkBox to ArborX conversion that would happen
     without expanding the boxes. */
  stk::search::Box<float> naiveArborXDomain = make_naive_arbor_x_box_from_double_sphere(domain);
  stk::search::Box<float> naiveArborXRange  = make_naive_arbor_x_box_from_double_sphere(range);
  std::vector<std::pair<stk::search::Box<float>, IdentProc>> naiveArborXDomainVector;
  naiveArborXDomainVector.push_back({naiveArborXDomain, IdentProc(0,0)});
  std::vector<std::pair<stk::search::Box<float>, IdentProc>> naiveArborXRangeVector;
  naiveArborXRangeVector.push_back({naiveArborXRange, IdentProc(0,0)});

  std::vector<std::pair<IdentProc, IdentProc>> searchResultsNaiveArborX;
  stk::search::coarse_search(naiveArborXDomainVector, naiveArborXRangeVector, stk::search::ARBORX, MPI_COMM_WORLD, searchResultsNaiveArborX);
  EXPECT_EQ(0u, searchResultsNaiveArborX.size());

  /*This case captures the "expanded" box behavior that is now the default when converting from StkBoxes to ArborX.*/
  stk::search::Sphere<SpherePrecision> domainExpanded(domainCenter, domainRadius);
  std::vector<std::pair<stk::search::Sphere<SpherePrecision>, IdentProc>> domainExpandedVector;
  domainExpandedVector.push_back({domainExpanded, IdentProc(0,0)});

  stk::search::Sphere<SpherePrecision> rangeExpanded(rangeCenter, rangeRadius);
  std::vector<std::pair<stk::search::Sphere<SpherePrecision>, IdentProc>> rangeExpandedVector;
  rangeExpandedVector.push_back({rangeExpanded, IdentProc(0,0)});

  std::vector<std::pair<IdentProc, IdentProc>> searchResultsExpanded;
  stk::search::coarse_search(domainExpandedVector, rangeExpandedVector, stk::search::ARBORX, MPI_COMM_WORLD, searchResultsExpanded);
  EXPECT_EQ(1u, searchResultsExpanded.size());

}

}
