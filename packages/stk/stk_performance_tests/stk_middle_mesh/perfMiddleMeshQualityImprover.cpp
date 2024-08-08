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
#include <stk_middle_mesh/mesh_boundary_snapper.hpp>
#include <stk_middle_mesh/incremental_mesh_boundary_snapper.hpp>
#include <stk_middle_mesh_util/create_stk_mesh.hpp>
#include <stk_unit_tests/stk_middle_mesh/util/meshes.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/timer.hpp>

namespace {

auto pi = std::atan(1) * 4;

template <typename T>
std::shared_ptr<stk::middle_mesh::mesh::Mesh> make_annulus_mesh(const int nelemR, const int nelemTheta, const double rIn,
                                                                const double rOut, double dtheta, MPI_Comm comm, T func2)
{
  assert(rIn > 0);
  assert(rOut > rIn);

  stk::middle_mesh::mesh::impl::MeshSpec spec;
  spec.numelX    = nelemR;
  spec.numelY    = nelemTheta;
  spec.xmin      = rIn;
  spec.xmax      = rOut;
  spec.ymin      = 0;
  spec.ymax      = 2 * pi;
  spec.yPeriodic = true;

  auto func = [&](const stk::middle_mesh::utils::Point& pt) {
    double r     = pt.x;
    double theta = pt.y;
    if (std::fmod(theta, 2 * pi) > 1e-12)
      theta += dtheta;

    double x = r * std::cos(theta);
    double y = r * std::sin(theta);
    double z = func2(x, y);
    stk::middle_mesh::utils::Point pt2(x, y, z);
    return pt2;
  };

  return stk::middle_mesh::mesh::impl::create_mesh(spec, func, comm);
}

std::shared_ptr<stk::middle_mesh::mesh::Mesh> make_annulus_mesh(const int nelemR, const int nelemTheta, const double rIn,
                                                                const double rOut, double dtheta, MPI_Comm comm = MPI_COMM_WORLD)
{
  auto f2 = [](const double x, const double y) { return 0.0; };
  return make_annulus_mesh(nelemR, nelemTheta, rIn, rOut, dtheta, comm, f2);
}


TEST(MiddleMeshQualityImprover, IncrementalBoundarySnapperAnnulusRotation)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  const int NUM_RUNS  = 3;
  const int NUM_ITERS = 1;

  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();
  for (unsigned run = 0; run < NUM_RUNS; ++run) {
    std::cout << "run = " << run << std::endl;
    auto mesh1 = make_annulus_mesh(35, 35, 0.5, 1.5, 0);
    auto mesh2 = make_annulus_mesh(71, 71, 0.5, 1.5, pi/32);

    batchTimer.start_batch_timer();
    auto bsnapper = stk::middle_mesh::mesh::impl::make_incremental_boundary_snapper(mesh1, mesh2, comm);
    bsnapper->snap();
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(MiddleMeshQualityImprover, IncrementalBoundarySnapperLargeElemCount)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  const int NUM_RUNS  = 3;
  const int NUM_ITERS = 1;

  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();
  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    stk::middle_mesh::mesh::impl::MeshSpec spec{300, 300, 0, 1, 0, 1};
    stk::middle_mesh::mesh::impl::MeshSpec spec2{30, 30, 0, 1, 0, 1};
    auto func = [&](stk::middle_mesh::utils::Point const& pt) { return pt; };

    auto mesh1 = create_mesh(spec, func);
    auto mesh2 = create_mesh(spec2, func);

    batchTimer.start_batch_timer();
    auto bsnapper = stk::middle_mesh::mesh::impl::make_incremental_boundary_snapper(mesh1, mesh2, comm);
    bsnapper->snap();
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

}
