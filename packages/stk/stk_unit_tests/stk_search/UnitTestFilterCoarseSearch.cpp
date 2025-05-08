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
// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <memory>
#include <math.h>                                     // for sqrt
#include <ostream>                                    // for basic_ostream::...
#include <stdexcept>                                  // for runtime_error
#include <string>                                     // for operator<<, string
#include <utility>                                    // for move
#include <vector>                                     // for vector, swap

#include <gtest/gtest.h>
#include "mpi.h"                                      // for MPI_COMM_WORLD

#include "stk_search/FilterCoarseSearch.hpp"
#include "stk_search/IdentProc.hpp"                   // for IdentProc
#include "stk_util/parallel/Parallel.hpp"             // for parallel_machin...
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_util/parallel/ParallelReduceBool.hpp"
#include <stk_unit_test_utils/getOption.h>            // for get_option

namespace {

class SearchFilterTester : public ::testing::Test {
 public:
  class SendMesh {
   public:
    SendMesh(SearchFilterTester& owner)
      : m_owner(owner)
    {
    }
    using EntityKey = int;
    using EntityProc = stk::search::IdentProc<EntityKey, unsigned>;

    stk::search::ObjectOutsideDomainPolicy get_extrapolate_option() const
    {
      return stk::search::ObjectOutsideDomainPolicy::EXTRAPOLATE;
    }

    void find_parametric_coords(const EntityKey k, const std::vector<double>& /*coords*/, std::vector<double>& /*parametricCoords*/,
                                double& parametricDistance, bool& isWithinParametricTolerance) const
    {
      // parametric tolerance applied here
      parametricDistance = m_owner.parametric_dist.at(k);
      isWithinParametricTolerance = (parametricDistance <= 1.0 + m_owner.parametricTolerance);
    }

    bool modify_search_outside_parametric_tolerance(const EntityKey /*k*/, const std::vector<double>& /*toCoords*/,
                                                    std::vector<double>& /*parametricCoords*/,
                                                    double& /*geometricDistanceSquared*/,
                                                    bool& /*isWithinGeometricTolerance*/) const
    {
      return false;
    }

    void coordinates(const EntityKey /*k*/, std::vector<double>& /*coords*/) const { }

    double get_closest_geometric_distance_squared(const EntityKey k, const std::vector<double>& /*coords*/) const { return m_owner.geometric_dist.at(k); }

    double get_distance_squared_from_centroid(const EntityKey k, const std::vector<double>& /*coords*/) const { return m_owner.geometric_dist.at(k); }

    SearchFilterTester& m_owner;
  };

  class RecvMesh {
   public:
    RecvMesh(SearchFilterTester& owner)
      : m_owner(owner)
    {
    }

    using EntityKey = int;
    using EntityProc = stk::search::IdentProc<EntityKey, unsigned>;

    stk::ParallelMachine comm() const { return MPI_COMM_WORLD; }

    double get_search_tolerance() const { return m_owner.geometricTolerance; }

    void coordinates(int /*id*/, std::vector<double>& /*coords*/) const { }

    double get_distance_from_nearest_node(int /*nodeId*/, const std::vector<double>& /*coords*/) { return 0.0; }

    SearchFilterTester& m_owner;
  };

  using Relation = std::pair<RecvMesh::EntityProc, SendMesh::EntityProc>;
  using RelationVec = std::vector<Relation> ;

  SearchFilterTester()
    : sendMesh(*this)
    , recvMesh(*this)
  {
  }

  RelationVec coarse_results;
  SendMesh sendMesh;
  RecvMesh recvMesh;
  stk::search::FilterCoarseSearchResultMap<RecvMesh> filteredSearchResults;

  void add_coarse_results(const std::vector<double>& par_dist, const std::vector<double>& geo_dist)
  {
    int num_matches = (int)par_dist.size();
    EXPECT_EQ(par_dist.size(), geo_dist.size());

    // Single recv node
    auto node = RecvMesh::EntityProc(1, 0);

    // Send element list
    for(int i = 0; i < num_matches; ++i) {
      auto elem = SendMesh::EntityProc(i + 1, 0);
      coarse_results.emplace_back(node, elem);
      parametric_dist[i + 1] = par_dist[i];
      geometric_dist[i + 1] = geo_dist[i];
    }
  }

  std::map<int, double> parametric_dist;
  std::map<int, double> geometric_dist;
  double geometricTolerance = 0.0;
  double parametricTolerance = 0.0;
};

TEST_F(SearchFilterTester, checkEarlyGeometricMatch)
{
  if(stk::parallel_machine_size(recvMesh.comm()) != 1) {
    return;
  }

  geometricTolerance = 0.5;
  parametricTolerance = 0.5;

  // mix of points within and out of parametric and geo tol
  add_coarse_results({ 2.5, 1.1, 2.0, 0.5 }, { 2.0, 0.2, 2.5, 0.0 });

  auto extrapolateOption = stk::search::ObjectOutsideDomainPolicy::ABORT;

  stk::search::impl::FilterCoarseSearchStats stats = stk::search::impl::filter_coarse_search_by_range(coarse_results,
                                                                                                      sendMesh, recvMesh,
                                                                                                      false, false,
                                                                                                      extrapolateOption,
                                                                                                      filteredSearchResults);

  EXPECT_EQ(1u, stats.numEntitiesWithinTolerance);
  EXPECT_EQ(0u, stats.numEntitiesOutsideTolerance);
  EXPECT_EQ(1u, coarse_results.size());
  EXPECT_EQ(4, coarse_results[0].second.id());
}

TEST_F(SearchFilterTester, checkAllGeometric)
{
  if(stk::parallel_machine_size(recvMesh.comm()) != 1) {
    return;
  }

  geometricTolerance = 0.5;
  parametricTolerance = 0.5;

  // all points outside parametric tol - pick closest geo
  add_coarse_results({ 2.5, 1.6, 2.0, 3.5 }, { 2.0, 0.2, 2.5, 0.7 });

  auto extrapolateOption = stk::search::ObjectOutsideDomainPolicy::ABORT;

  stk::search::impl::FilterCoarseSearchStats stats = stk::search::impl::filter_coarse_search_by_range(coarse_results,
                                                                                                      sendMesh, recvMesh,
                                                                                                      false, false,
                                                                                                      extrapolateOption,
                                                                                                      filteredSearchResults);

  EXPECT_EQ(1u, stats.numEntitiesWithinTolerance);
  EXPECT_EQ(0u, stats.numEntitiesOutsideTolerance);
  EXPECT_EQ(1u, coarse_results.size());
  EXPECT_EQ(2, coarse_results[0].second.id());
}

TEST_F(SearchFilterTester, checkAllOutsideTolerance)
{
  if(stk::parallel_machine_size(recvMesh.comm()) != 1) {
    return;
  }

  geometricTolerance = 0.5;
  parametricTolerance = 0.5;

  // all points outside both parametric and geometric tolerance
  add_coarse_results({ 2.5, 1.6, 2.0, 3.5 }, { 2.0, 0.75, 2.5, 0.7 });

  auto extrapolateOption = stk::search::ObjectOutsideDomainPolicy::ABORT;

  stk::search::impl::FilterCoarseSearchStats stats = stk::search::impl::filter_coarse_search_by_range(coarse_results,
                                                                                                      sendMesh, recvMesh,
                                                                                                      false, false,
                                                                                                      extrapolateOption,
                                                                                                      filteredSearchResults);

  EXPECT_EQ(0u, stats.numEntitiesWithinTolerance);
  EXPECT_EQ(1u, stats.numEntitiesOutsideTolerance);
  EXPECT_EQ(1u, coarse_results.size());
  EXPECT_EQ(4, coarse_results[0].second.id());
}


}
