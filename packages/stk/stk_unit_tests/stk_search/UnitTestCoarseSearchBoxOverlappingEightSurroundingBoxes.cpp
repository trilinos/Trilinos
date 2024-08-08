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

template<class InnerBoundingBoxType, class OuterBoundingBoxType>
void runBoxOverlappingEightSurroundingBoxes(stk::search::SearchMethod searchMethod, const double radius,
                                            const unsigned numExpectedResults)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProc = stk::parallel_machine_size(comm);
  int procId = stk::parallel_machine_rank(comm);

  std::vector<std::pair<OuterBoundingBoxType, IdentProc>> boxVector1;
  if (procId == 0) {
    boxVector1.push_back(stk::unit_test_util::generateBoundingVolume<OuterBoundingBoxType>(0, 0, 0, radius, 1, procId));
    boxVector1.push_back(stk::unit_test_util::generateBoundingVolume<OuterBoundingBoxType>(1, 0, 0, radius, 2, procId));
    boxVector1.push_back(stk::unit_test_util::generateBoundingVolume<OuterBoundingBoxType>(2, 0, 0, radius, 3, procId));
    boxVector1.push_back(stk::unit_test_util::generateBoundingVolume<OuterBoundingBoxType>(0, 1, 0, radius, 4, procId));
    //skip middle one
    boxVector1.push_back(stk::unit_test_util::generateBoundingVolume<OuterBoundingBoxType>(2, 1, 0, radius, 6, procId));
    boxVector1.push_back(stk::unit_test_util::generateBoundingVolume<OuterBoundingBoxType>(0, 2, 0, radius, 7, procId));
    boxVector1.push_back(stk::unit_test_util::generateBoundingVolume<OuterBoundingBoxType>(1, 2, 0, radius, 8, procId));
    boxVector1.push_back(stk::unit_test_util::generateBoundingVolume<OuterBoundingBoxType>(2, 2, 0, radius, 9, procId));
  }

  std::vector<std::pair<InnerBoundingBoxType, IdentProc>> boxVector2;
  if (procId == numProc-1) {
    boxVector2.push_back(stk::unit_test_util::generateBoundingVolume<InnerBoundingBoxType>(1, 1, 0, radius, 5, procId));
  }

  std::vector< std::pair<IdentProc, IdentProc> > boxIdPairResults;
  stk::search::coarse_search(boxVector1, boxVector2, searchMethod, comm, boxIdPairResults);

  if (procId == 0 || procId == numProc-1) {
    if (numExpectedResults != boxIdPairResults.size()) {
      for (size_t i = 0; i < boxIdPairResults.size(); ++i) {
        std::cerr << boxIdPairResults[i].first << ", " << boxIdPairResults[i].second << std::endl;
      }
    }
    ASSERT_EQ(numExpectedResults, boxIdPairResults.size());
  }
}

template<class InnerBoxType, class OuterBoxType>
void device_runBoxOverlappingEightSurroundingBoxes(stk::search::SearchMethod searchMethod, const double radius,
                                                         const unsigned numExpectedResults)
{

  MPI_Comm comm = MPI_COMM_WORLD;
  int numProc = stk::parallel_machine_size(comm);
  int procId = stk::parallel_machine_rank(comm);

  using InnerBoxIdentProcType = stk::search::BoxIdentProc<InnerBoxType, IdentProc>;
  using OuterBoxIdentProcType = stk::search::BoxIdentProc<OuterBoxType, IdentProc>;

  auto domain = Kokkos::View<OuterBoxIdentProcType*, stk::ngp::ExecSpace>("domain box-ident", 0);
  auto range = Kokkos::View<InnerBoxIdentProcType*, stk::ngp::ExecSpace>("range box-ident", 0);

  if (procId == 0) {
    Kokkos::resize(Kokkos::WithoutInitializing, domain, 8);
    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
      KOKKOS_LAMBDA(const unsigned & i) {
        domain[0] = stk::unit_test_util::device_generateBoxIdentProc<OuterBoxType, IdentProc>(0, 0, 0, radius, 1, procId);
        domain[1] = stk::unit_test_util::device_generateBoxIdentProc<OuterBoxType, IdentProc>(1, 0, 0, radius, 2, procId);
        domain[2] = stk::unit_test_util::device_generateBoxIdentProc<OuterBoxType, IdentProc>(2, 0, 0, radius, 3, procId);
        domain[3] = stk::unit_test_util::device_generateBoxIdentProc<OuterBoxType, IdentProc>(0, 1, 0, radius, 4, procId);
        // Skip middle box
        domain[4] = stk::unit_test_util::device_generateBoxIdentProc<OuterBoxType, IdentProc>(2, 1, 0, radius, 6, procId);
        domain[5] = stk::unit_test_util::device_generateBoxIdentProc<OuterBoxType, IdentProc>(1, 2, 0, radius, 7, procId);
        domain[6] = stk::unit_test_util::device_generateBoxIdentProc<OuterBoxType, IdentProc>(2, 2, 0, radius, 8, procId);
        domain[7] = stk::unit_test_util::device_generateBoxIdentProc<OuterBoxType, IdentProc>(0, 2, 0, radius, 9, procId);

      });
  }

  if (procId == numProc-1) {
    Kokkos::resize(Kokkos::WithoutInitializing, range, 1);
    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
      KOKKOS_LAMBDA(const unsigned & i) {
        range[0] = stk::unit_test_util::device_generateBoxIdentProc<InnerBoxType, IdentProc>(1, 1, 0, radius, 5, procId);
      });
  }

  auto intersections = Kokkos::View<IdentProcIntersection*, stk::ngp::ExecSpace>("intersections", 0);

  stk::search::coarse_search(domain, range, searchMethod, comm, intersections);

  Kokkos::View<IdentProcIntersection*>::HostMirror hostIntersections = Kokkos::create_mirror_view(intersections);
  Kokkos::deep_copy(hostIntersections, intersections);

  if (procId == 0 || procId == numProc - 1) {
    ASSERT_EQ(intersections.extent(0), numExpectedResults);
  }
  else {
    ASSERT_EQ(intersections.extent(0), 0u);
  }
}

TEST(CoarseSearchCorrectness, SphereOverlappingEightSurroundingSpheres_KDTREE)
{
  const double radius = 0.708;
  const unsigned numExpectedResults = 8;
  runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, SphereOverlappingFourOfEightSurroundingSpheres_KDTREE)
{
  const double radius = 0.706;
  const unsigned numExpectedResults = 4;
  runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, SphereOverlappingEightSurroundingSpheres_MORTON_LBVH)
{
  const double radius = 0.708;
  const unsigned numExpectedResults = 8;
  runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, SphereOverlappingFourOfEightSurroundingSpheres_MORTON_LBVH)
{
  const double radius = 0.706;
  const unsigned numExpectedResults = 4;
  runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}


TEST(CoarseSearchCorrectness, SphereOverlappingNoSurroundingPoints_KDTREE)
{
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, SphereOverlappingFourSurroundingPoints_KDTREE)
{
  const double radius = 1.41;
  const unsigned numExpectedResults = 4;
  runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, SphereOverlappingEightSurroundingPoints_KDTREE)
{
  const double radius = 1.42;
  const unsigned numExpectedResults = 8;
  runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, SphereOverlappingNoSurroundingPoints_MORTON_LBVH)
{
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, SphereOverlappingFourSurroundingPoints_MORTON_LBVH)
{
  const double radius = 1.41;
  const unsigned numExpectedResults = 4;
  runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, SphereOverlappingEightSurroundingPoints_MORTON_LBVH)
{
  const double radius = 1.42;
  const unsigned numExpectedResults = 8;
  runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}


TEST(CoarseSearchCorrectness, BoxOverlappingNoSurroundingPoints_KDTREE)
{
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, BoxOverlappingEightSurroundingPoints_KDTREE)
{
  const double radius = 1.01;
  const unsigned numExpectedResults = 8;
  runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, BoxOverlappingNoSurroundingPoints_MORTON_LBVH)
{
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, BoxOverlappingEightSurroundingPoints_MORTON_LBVH)
{
  const double radius = 1.01;
  const unsigned numExpectedResults = 8;
  runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}


TEST(CoarseSearchCorrectness, PointOverlappingNoSurroundingBoxes_KDTREE)
{
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, PointOverlappingEightSurroundingBoxes_KDTREE)
{
  const double radius = 1.01;
  const unsigned numExpectedResults = 8;
  runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, PointOverlappingNoSurroundingBoxes_MORTON_LBVH)
{
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, PointOverlappingEightSurroundingBoxes_MORTON_LBVH)
{
  const double radius = 1.01;
  const unsigned numExpectedResults = 8;
  runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}


TEST(CoarseSearchCorrectness, BoxOverlappingNoSurroundingBoxes_KDTREE)
{
  const double radius = 0.49;
  const unsigned numExpectedResults = 0;
  runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, BoxOverlappingEightSurroundingBoxes_KDTREE)
{
  const double radius = 0.51;
  const unsigned numExpectedResults = 8;
  runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, BoxOverlappingNoSurroundingBoxes_MORTON_LBVH)
{
  const double radius = 0.49;
  const unsigned numExpectedResults = 0;
  runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, BoxOverlappingEightSurroundingBoxes_MORTON_LBVH)
{
  const double radius = 0.51;
  const unsigned numExpectedResults = 8;
  runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}


TEST(CoarseSearchCorrectness, SphereOverlappingEightSurroundingSpheres_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 0.708;
  const unsigned numExpectedResults = 8;
  runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::ARBORX, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, SphereOverlappingNoSurroundingPoints_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::ARBORX, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, SphereOverlappingFourSurroundingPoints_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 1.41;
  const unsigned numExpectedResults = 4;
  runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::ARBORX, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, SphereOverlappingEightSurroundingPoints_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 1.42;
  const unsigned numExpectedResults = 8;
  runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::ARBORX, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, SphereOverlappingFourOfEightSurroundingSpheres_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 0.706;
  const unsigned numExpectedResults = 4;
  runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::ARBORX, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, BoxOverlappingNoSurroundingPoints_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::ARBORX, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, BoxOverlappingEightSurroundingPoints_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 1.01;
  const unsigned numExpectedResults = 8;
  runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::ARBORX, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, PointOverlappingNoSurroundingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, PointOverlappingEightSurroundingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 1.01;
  const unsigned numExpectedResults = 8;
  runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, BoxOverlappingNoSurroundingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 0.49;
  const unsigned numExpectedResults = 0;
  runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, BoxOverlappingEightSurroundingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 0.51;
  const unsigned numExpectedResults = 8;
  runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
  device_runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
}

template <class InnerBoxType, class OuterBoxType>
void host_local_runBoxOverlappingEightSurroundingBoxes(
    stk::search::SearchMethod searchMethod, const double radius, const unsigned numExpectedResults)
{
  using InnerBoxIdentType = std::pair<InnerBoxType, Ident>;
  using OuterBoxIdentType = std::pair<OuterBoxType, Ident>;

  std::vector<OuterBoxIdentType> domain;
  std::vector<InnerBoxIdentType> range;

  domain.push_back(stk::unit_test_util::box_ident_to_pair(
      stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(0, 0, 0, radius, 1)));
  domain.push_back(stk::unit_test_util::box_ident_to_pair(
      stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(1, 0, 0, radius, 2)));
  domain.push_back(stk::unit_test_util::box_ident_to_pair(
      stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(2, 0, 0, radius, 3)));
  domain.push_back(stk::unit_test_util::box_ident_to_pair(
      stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(0, 1, 0, radius, 4)));
  // Skip middle box
  domain.push_back(stk::unit_test_util::box_ident_to_pair(
      stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(2, 1, 0, radius, 6)));
  domain.push_back(stk::unit_test_util::box_ident_to_pair(
      stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(1, 2, 0, radius, 7)));
  domain.push_back(stk::unit_test_util::box_ident_to_pair(
      stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(2, 2, 0, radius, 8)));
  domain.push_back(stk::unit_test_util::box_ident_to_pair(
      stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(0, 2, 0, radius, 9)));

  range.push_back(stk::unit_test_util::box_ident_to_pair(
      stk::unit_test_util::device_generateBoxIdent<InnerBoxType, Ident>(1, 1, 0, radius, 5)));

  LocalSearchResults intersections;

  stk::search::local_coarse_search(domain, range, searchMethod, intersections);

  ASSERT_EQ(intersections.size(), numExpectedResults);
}

template<class InnerBoxType, class OuterBoxType>
void device_local_runBoxOverlappingEightSurroundingBoxes(stk::search::SearchMethod searchMethod, const double radius,
                                                         const unsigned numExpectedResults)
{
  using InnerBoxIdentType = stk::search::BoxIdent<InnerBoxType, Ident>;
  using OuterBoxIdentType = stk::search::BoxIdent<OuterBoxType, Ident>;

  auto domain = Kokkos::View<OuterBoxIdentType*, stk::ngp::ExecSpace>("domain box-ident", 8);
  auto range = Kokkos::View<InnerBoxIdentType*, stk::ngp::ExecSpace>("range box-ident", 1);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
    KOKKOS_LAMBDA(const unsigned & i) {
      domain[0] = stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(0, 0, 0, radius, 1);
      domain[1] = stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(1, 0, 0, radius, 2);
      domain[2] = stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(2, 0, 0, radius, 3);
      domain[3] = stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(0, 1, 0, radius, 4);
      // Skip middle box
      domain[4] = stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(2, 1, 0, radius, 6);
      domain[5] = stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(1, 2, 0, radius, 7);
      domain[6] = stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(2, 2, 0, radius, 8);
      domain[7] = stk::unit_test_util::device_generateBoxIdent<OuterBoxType, Ident>(0, 2, 0, radius, 9);

      range[0] = stk::unit_test_util::device_generateBoxIdent<InnerBoxType, Ident>(1, 1, 0, radius, 5);
    });

  auto intersections = Kokkos::View<IdentIntersection*, stk::ngp::ExecSpace>("intersections", 0);

  stk::search::local_coarse_search(domain, range, searchMethod, intersections);

  Kokkos::View<IdentIntersection*>::HostMirror hostIntersections = Kokkos::create_mirror_view(intersections);
  Kokkos::deep_copy(hostIntersections, intersections);

  ASSERT_EQ(intersections.extent(0), numExpectedResults);
}





TEST(CoarseSearchCorrectness, Local_SphereOverlappingEightSurroundingSpheres_KDTREE)
{
  const double radius = 0.708;
  const unsigned numExpectedResults = 8;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Local_SphereOverlappingFourOfEightSurroundingSpheres_KDTREE)
{
  const double radius = 0.706;
  const unsigned numExpectedResults = 4;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Local_SphereOverlappingNoSurroundingPoints_KDTREE)
{
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Local_SphereOverlappingFourSurroundingPoints_KDTREE)
{
  const double radius = 1.41;
  const unsigned numExpectedResults = 4;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Local_SphereOverlappingEightSurroundingPoints_KDTREE)
{
  const double radius = 1.42;
  const unsigned numExpectedResults = 8;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::KDTREE, radius, numExpectedResults);
}


TEST(CoarseSearchCorrectness, Local_BoxOverlappingNoSurroundingPoints_KDTREE)
{
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  host_local_runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Local_BoxOverlappingEightSurroundingPoints_KDTREE)
{
  const double radius = 1.01;
  const unsigned numExpectedResults = 8;
  host_local_runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::KDTREE, radius, numExpectedResults);
}


TEST(CoarseSearchCorrectness, Local_PointOverlappingNoSurroundingBoxes_KDTREE)
{
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  host_local_runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Local_PointOverlappingEightSurroundingBoxes_KDTREE)
{
  const double radius = 1.01;
  const unsigned numExpectedResults = 8;
  host_local_runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::KDTREE, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Ngp_Local_SphereOverlappingEightSurroundingSpheres_MORTON_LBVH)
{
  const double radius = 0.708;
  const unsigned numExpectedResults = 8;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Ngp_Local_SphereOverlappingFourOfEightSurroundingSpheres_MORTON_LBVH)
{
  const double radius = 0.706;
  const unsigned numExpectedResults = 4;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}


TEST(CoarseSearchCorrectness, Ngp_Local_SphereOverlappingNoSurroundingPoints_MORTON_LBVH)
{
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Ngp_Local_SphereOverlappingFourSurroundingPoints_MORTON_LBVH)
{
  const double radius = 1.41;
  const unsigned numExpectedResults = 4;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Ngp_Local_SphereOverlappingEightSurroundingPoints_MORTON_LBVH)
{
  const double radius = 1.42;
  const unsigned numExpectedResults = 8;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}


TEST(CoarseSearchCorrectness, Ngp_Local_BoxOverlappingNoSurroundingPoints_MORTON_LBVH)
{
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  host_local_runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Ngp_Local_BoxOverlappingEightSurroundingPoints_MORTON_LBVH)
{
  const double radius = 1.01;
  const unsigned numExpectedResults = 8;
  host_local_runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}


TEST(CoarseSearchCorrectness, Ngp_Local_PointOverlappingNoSurroundingBoxes_MORTON_LBVH)
{
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  host_local_runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Ngp_Local_PointOverlappingEightSurroundingBoxes_MORTON_LBVH)
{
  const double radius = 1.01;
  const unsigned numExpectedResults = 8;
  host_local_runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}


TEST(CoarseSearchCorrectness, Ngp_Local_BoxOverlappingNoSurroundingBoxes_MORTON_LBVH)
{
  const double radius = 0.49;
  const unsigned numExpectedResults = 0;
  host_local_runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Ngp_Local_BoxOverlappingEightSurroundingBoxes_MORTON_LBVH)
{
  const double radius = 0.51;
  const unsigned numExpectedResults = 8;
  host_local_runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::MORTON_LBVH, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Ngp_Local_SphereOverlappingEightSurroundingSpheres_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 0.708;
  const unsigned numExpectedResults = 8;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::ARBORX, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Ngp_Local_SphereOverlappingFourOfEightSurroundingSpheres_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 0.706;
  const unsigned numExpectedResults = 4;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::ARBORX, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::ARBORX, radius, numExpectedResults);
}


TEST(CoarseSearchCorrectness, Ngp_Local_SphereOverlappingNoSurroundingPoints_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::ARBORX, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Ngp_Local_SphereOverlappingFourSurroundingPoints_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 1.41;
  const unsigned numExpectedResults = 4;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::ARBORX, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Ngp_Local_SphereOverlappingEightSurroundingPoints_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 1.42;
  const unsigned numExpectedResults = 8;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::ARBORX, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::ARBORX, radius, numExpectedResults);
}


TEST(CoarseSearchCorrectness, Ngp_Local_BoxOverlappingNoSurroundingPoints_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  host_local_runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::ARBORX, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Ngp_Local_BoxOverlappingEightSurroundingPoints_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 1.01;
  const unsigned numExpectedResults = 8;
  host_local_runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::ARBORX, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::ARBORX, radius, numExpectedResults);
}


TEST(CoarseSearchCorrectness, Ngp_Local_PointOverlappingNoSurroundingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 0.99;
  const unsigned numExpectedResults = 0;
  host_local_runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Ngp_Local_PointOverlappingEightSurroundingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 1.01;
  const unsigned numExpectedResults = 8;
  host_local_runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
}


TEST(CoarseSearchCorrectness, Ngp_Local_BoxOverlappingNoSurroundingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 0.49;
  const unsigned numExpectedResults = 0;
  host_local_runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
}

TEST(CoarseSearchCorrectness, Ngp_Local_BoxOverlappingEightSurroundingBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  const double radius = 0.51;
  const unsigned numExpectedResults = 8;
  host_local_runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
  device_local_runBoxOverlappingEightSurroundingBoxes<StkBox,StkBox>(stk::search::ARBORX, radius, numExpectedResults);
}

}
