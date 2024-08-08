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

#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_middle_mesh/create_mesh.hpp>
#include <stk_middle_mesh_util/create_stk_mesh.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/timer.hpp>

namespace {

static constexpr int VERT_RANK = 0;
static constexpr int EDGE_RANK = 1;
static constexpr int ELEM_RANK = 2;

TEST(MiddleMeshOps, GetAndDeleteDown)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::middle_mesh::utils::impl::comm_size(comm) > 1)
    GTEST_SKIP();

  const int NUM_RUNS  = 40;
  const int NUM_ITERS = 1;

  stk::middle_mesh::mesh::impl::MeshSpec spec{1000, 1000, 0, 1, 0, 1};
  auto func = [&](stk::middle_mesh::utils::Point const& pt) { return pt; };

  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();
  for (unsigned i = 0; i < NUM_RUNS; ++i) {
    auto mesh = stk::middle_mesh::mesh::impl::create_mesh(spec, func);

    batchTimer.start_batch_timer();

    auto& elements = mesh->get_mesh_entities(ELEM_RANK);

    for (auto elem : elements) {
      if (elem == nullptr) { continue; }
      [[maybe_unused]] auto elemId = elem->get_id();
      auto numEdges = elem->count_down();

      for (auto e = numEdges-1; e > 0; --e) {
        auto edge = elem->get_down(e);
        if (edge == nullptr) { continue; }
        [[maybe_unused]] auto eId = edge->get_id();
        auto numVerts = edge->count_down();

        for (auto v = numVerts-1; v > 0; --v) {
          auto vert = edge->get_down(v);
          if (vert == nullptr) { continue; }
          [[maybe_unused]] auto vId = vert->get_id();
          edge->delete_down(v);
        }
        elem->delete_down(e);
      }
    }
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(MiddleMeshOps, GetAndDeleteUp)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::middle_mesh::utils::impl::comm_size(comm) > 1)
    GTEST_SKIP();

  const int NUM_RUNS  = 40;
  const int NUM_ITERS = 1;

  stk::middle_mesh::mesh::impl::MeshSpec spec{1000, 1000, 0, 1, 0, 1};
  auto func = [&](stk::middle_mesh::utils::Point const& pt) { return pt; };

  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();
  for (unsigned i = 0; i < NUM_RUNS; ++i) {
    auto mesh = stk::middle_mesh::mesh::impl::create_mesh(spec, func);

    batchTimer.start_batch_timer();

    auto& verts = mesh->get_mesh_entities(VERT_RANK);

    for (auto vert : verts) {
      if (vert == nullptr) { continue; }
      [[maybe_unused]] auto vId = vert->get_id();
      auto numEdges = vert->count_up();

      for (auto e = numEdges-1; e > 0; --e) {
        auto edge = vert->get_up(e);
        if (edge == nullptr) { continue; }
        [[maybe_unused]] auto eId = edge->get_id();
        auto numElems = edge->count_up();

        for (auto el = numElems-1; el > 0; --el) {
          auto elem = edge->get_up(el);
          if (elem == nullptr) { continue; }
          [[maybe_unused]] auto elemId = elem->get_id();
          edge->delete_up(el);
        }
        vert->delete_up(e);
      }
    }
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(MiddleMeshOps, GetUpward)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::middle_mesh::utils::impl::comm_size(comm) > 1)
    GTEST_SKIP();

  const int NUM_RUNS  = 25;
  const int NUM_ITERS = 1;

  stk::middle_mesh::mesh::impl::MeshSpec spec{1000, 1000, 0, 1, 0, 1};
  auto func = [&](stk::middle_mesh::utils::Point const& pt) { return pt; };

  std::vector<stk::middle_mesh::mesh::MeshEntityPtr> edges, elems;

  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  auto mesh = stk::middle_mesh::mesh::impl::create_mesh(spec, func);

  for (unsigned i = 0; i < NUM_RUNS; ++i) {
    batchTimer.start_batch_timer();
    auto& verts = mesh->get_mesh_entities(VERT_RANK);

    for (auto vert : verts) {
      get_upward(vert, EDGE_RANK, edges);
      get_upward(vert, ELEM_RANK, elems);

      for (auto edge : edges) {
        get_upward(edge, ELEM_RANK, elems);
      }
    }
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(MiddleMeshOps, GetDownward)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::middle_mesh::utils::impl::comm_size(comm) > 1)
    GTEST_SKIP();

  const int NUM_RUNS  = 50;
  const int NUM_ITERS = 1;

  stk::middle_mesh::mesh::impl::MeshSpec spec{1000, 1000, 0, 1, 0, 1};
  auto func = [&](stk::middle_mesh::utils::Point const& pt) { return pt; };

  std::vector<stk::middle_mesh::mesh::MeshEntityPtr> edges, verts;

  edges.resize(stk::middle_mesh::mesh::MAX_DOWN);
  verts.resize(stk::middle_mesh::mesh::MAX_DOWN);

  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  auto mesh = stk::middle_mesh::mesh::impl::create_mesh(spec, func);

  for (unsigned i = 0; i < NUM_RUNS; ++i) {
    batchTimer.start_batch_timer();
    auto& elems = mesh->get_mesh_entities(ELEM_RANK);

    for (auto elem : elems) {
      get_downward(elem, EDGE_RANK, edges.data());
      get_downward(elem, VERT_RANK, verts.data());

      for (auto edge : edges) {
        get_downward(edge, VERT_RANK, verts.data());
      }
    }
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

}
