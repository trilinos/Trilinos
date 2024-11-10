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

enum Axis
{
  xDim, yDim, zDim
};

template<class BoundingBoxType>
void runLineOfBoundingBoxes(stk::search::SearchMethod searchMethod, enum Axis axis)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int procId = stk::parallel_machine_rank(comm);

  const double radius = 0.708;
  const double distanceBetweenCenters = 1.0;

  const double paramCoord = procId * distanceBetweenCenters;

  std::vector<std::pair<BoundingBoxType, IdentProc>> boxVector1;
  std::vector<std::pair<BoundingBoxType, IdentProc>> boxVector2;
  if (procId % 2 == 0) {
    switch(axis)
    {
    case xDim:
      boxVector1.push_back(stk::unit_test_util::generateBoundingVolume<BoundingBoxType>(paramCoord, 0, 0, radius, 1, procId));
      break;
    case yDim:
      boxVector1.push_back(stk::unit_test_util::generateBoundingVolume<BoundingBoxType>(0, paramCoord, 0, radius, 1, procId));
      break;
    case zDim:
      boxVector1.push_back(stk::unit_test_util::generateBoundingVolume<BoundingBoxType>(0, 0, paramCoord, radius, 1, procId));
      break;
    }
  }
  else {
    switch(axis)
    {
    case xDim:
      boxVector2.push_back(stk::unit_test_util::generateBoundingVolume<BoundingBoxType>(paramCoord, 0, 0, radius, 1, procId));
      break;
    case yDim:
      boxVector2.push_back(stk::unit_test_util::generateBoundingVolume<BoundingBoxType>(0, paramCoord, 0, radius, 1, procId));
      break;
    case zDim:
      boxVector2.push_back(stk::unit_test_util::generateBoundingVolume<BoundingBoxType>(0, 0, paramCoord, radius, 1, procId));
      break;
    }
  }

  SearchResults boxIdPairResults;
  stk::search::coarse_search(boxVector1, boxVector2, searchMethod, comm, boxIdPairResults);

  int numProcs = stk::parallel_machine_size(comm);

  unsigned numExpectedResults = 2;
  bool doOwnFirstOrLastSphereInLine = procId == 0 || procId == numProcs-1;

  if (doOwnFirstOrLastSphereInLine) {
    numExpectedResults = 1;
  }

  if (numProcs == 1) {
    numExpectedResults = 0;
  }

  EXPECT_EQ(numExpectedResults, boxIdPairResults.size()) << "on proc id " << procId;
}

template<class BoundingBoxType>
void device_runLineOfBoundingBoxes(stk::search::SearchMethod searchMethod, enum Axis axis)
{

  using InnerBoxIdentProcType = stk::search::BoxIdentProc<BoundingBoxType, IdentProc>;
  using OuterBoxIdentProcType = stk::search::BoxIdentProc<BoundingBoxType, IdentProc>;

  MPI_Comm comm = MPI_COMM_WORLD;
  int procId = stk::parallel_machine_rank(comm);

  const double radius = 0.708;
  const double distanceBetweenCenters = 1.0;

  const double paramCoord = procId * distanceBetweenCenters;

  auto domain = Kokkos::View<OuterBoxIdentProcType*, stk::ngp::ExecSpace>("domain box-ident", 1);
  auto range = Kokkos::View<InnerBoxIdentProcType*, stk::ngp::ExecSpace>("range box-ident", 1);

  if (procId % 2 == 0) {
    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
      KOKKOS_LAMBDA(const unsigned & i) {
        switch(axis) 
        {
        case xDim:
          domain[0] = stk::unit_test_util::device_generateBoxIdentProc<BoundingBoxType, IdentProc>(paramCoord, 0, 0, radius, 1, procId);
          break;
        case yDim:
          domain[0] = stk::unit_test_util::device_generateBoxIdentProc<BoundingBoxType, IdentProc>(0, paramCoord, 0, radius, 1, procId);
          break;
        case zDim:
          domain[0] = stk::unit_test_util::device_generateBoxIdentProc<BoundingBoxType, IdentProc>(0, 0, paramCoord, radius, 1, procId);
          break;
        }
      });
  }

  else {
    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
      KOKKOS_LAMBDA(const unsigned & i) {
        switch(axis) 
        {
        case xDim:
          range[0] = stk::unit_test_util::device_generateBoxIdentProc<BoundingBoxType, IdentProc>(paramCoord, 0, 0, radius, 1, procId);
          break;
        case yDim:
          range[0] = stk::unit_test_util::device_generateBoxIdentProc<BoundingBoxType, IdentProc>(0, paramCoord, 0, radius, 1, procId);
          break;
        case zDim:
          range[0] = stk::unit_test_util::device_generateBoxIdentProc<BoundingBoxType, IdentProc>(0, 0, paramCoord, radius, 1, procId);
          break;
        }
      });
  }

  auto searchResults = Kokkos::View<IdentProcIntersection*, stk::ngp::ExecSpace>("searchResults", 0);
  stk::search::coarse_search(domain, range, searchMethod, comm, searchResults);

  int numProcs = stk::parallel_machine_size(comm);

  unsigned numExpectedResults = 2;
  bool doOwnFirstOrLastSphereInLine = procId == 0 || procId == numProcs-1;

  if (doOwnFirstOrLastSphereInLine) {
    numExpectedResults = 1;
  }

  if (numProcs == 1) {
    numExpectedResults = 0;
  }

  EXPECT_EQ(numExpectedResults, searchResults.size()) << "on proc id " << procId;
}

TEST(CoarseSearchCorrectness, LineOfSpheres_KDTREE)
{
  runLineOfBoundingBoxes<Sphere>(stk::search::KDTREE, xDim);
}

TEST(CoarseSearchCorrectness, LineOfBoxes_KDTREE)
{
  runLineOfBoundingBoxes<StkBox>(stk::search::KDTREE, yDim);
}

TEST(CoarseSearchCorrectness, LineOfSpheresZDimension_KDTREE)
{
  runLineOfBoundingBoxes<Sphere>(stk::search::KDTREE, zDim);
}


TEST(CoarseSearchCorrectness, LineOfSpheres_MORTON_LBVH)
{
  runLineOfBoundingBoxes<Sphere>(stk::search::MORTON_LBVH, xDim);
  device_runLineOfBoundingBoxes<Sphere>(stk::search::MORTON_LBVH, xDim);
}

TEST(CoarseSearchCorrectness, LineOfBoxes_MORTON_LBVH)
{
  runLineOfBoundingBoxes<StkBox>(stk::search::MORTON_LBVH, yDim);
  device_runLineOfBoundingBoxes<StkBox>(stk::search::MORTON_LBVH, yDim);
}

TEST(CoarseSearchCorrectness, LineOfSpheresZDimension_MORTON_LBVH)
{
  runLineOfBoundingBoxes<Sphere>(stk::search::MORTON_LBVH, zDim);
  device_runLineOfBoundingBoxes<Sphere>(stk::search::MORTON_LBVH, zDim);
}

TEST(CoarseSearchCorrectness, LineOfSpheres_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  runLineOfBoundingBoxes<Sphere>(stk::search::ARBORX, xDim);
  device_runLineOfBoundingBoxes<Sphere>(stk::search::ARBORX, xDim);
}

TEST(CoarseSearchCorrectness, LineOfBoxes_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  runLineOfBoundingBoxes<StkBox>(stk::search::ARBORX, yDim);
  device_runLineOfBoundingBoxes<StkBox>(stk::search::ARBORX, yDim);
}

TEST(CoarseSearchCorrectness, LineOfSpheresZDimension_ARBORX)
{
#ifndef STK_HAS_ARBORX
  GTEST_SKIP();
#endif
  runLineOfBoundingBoxes<Sphere>(stk::search::ARBORX, zDim);
  device_runLineOfBoundingBoxes<Sphere>(stk::search::ARBORX, zDim);
}

}
