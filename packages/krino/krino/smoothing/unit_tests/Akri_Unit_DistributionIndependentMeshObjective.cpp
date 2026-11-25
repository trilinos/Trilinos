// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// The value and gradient of a mesh objective, which steer the optimizer that smooths the mesh, do not
// depend on how the mesh is decomposed across processors or on its ids: each element's contribution is
// rounded to float, so that the sums over elements, and over the elements of each node (with the
// processors' parts added), are exact while the contributions span less than the double mantissa (see
// round_sum_for_distribution_independence).

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include <gtest/gtest.h>
#include <mpi.h>

#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Objective_MeanRatioNorm.hpp>
#include <Akri_Objective_Mesh.hpp>
#include <Akri_Smoothing_Utils.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino
{

namespace
{

// A node's coordinates and the objective's gradient there.
using NodeCoordinatesAndGradient = std::array<double, 6>;

struct MeshObjectiveValueAndGradient
{
  double value;
  std::vector<NodeCoordinatesAndGradient> sortedNodesGradients; // of every node, on every processor
};

// The owned entries, from every processor, sorted by the node coordinates.
std::vector<NodeCoordinatesAndGradient> gather_sorted(const std::vector<NodeCoordinatesAndGradient> & owned)
{
  int numProcs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  const int localCount = 6 * int(owned.size());
  std::vector<int> counts(numProcs), offsets(numProcs, 0);
  MPI_Allgather(&localCount, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
  for (int p = 1; p < numProcs; ++p)
    offsets[p] = offsets[p - 1] + counts[p - 1];
  std::vector<double> all(offsets.back() + counts.back());
  MPI_Allgatherv(owned.empty() ? nullptr : owned.front().data(), localCount, MPI_DOUBLE,
      all.data(), counts.data(), offsets.data(), MPI_DOUBLE, MPI_COMM_WORLD);
  std::vector<NodeCoordinatesAndGradient> gathered(all.size() / 6);
  for (size_t i = 0; i < gathered.size(); ++i)
    std::copy(all.begin() + 6 * i, all.begin() + 6 * (i + 1), gathered[i].begin());
  std::sort(gathered.begin(), gathered.end(),
      [](const NodeCoordinatesAndGradient & a, const NodeCoordinatesAndGradient & b)
      { return std::lexicographical_compare(a.begin(), a.begin() + 3, b.begin(), b.begin() + 3); });
  return gathered;
}

// A BCC tet mesh of [-1,1]^3, as generated or with its ids reversed and decomposed across the
// processors differently, with its nodes moved so that its elements range from thin (near z = 0) to
// large, and from well shaped to poor.
// The mean-ratio-norm objective's value and gradient at the node coordinates.
MeshObjectiveValueAndGradient mean_ratio_norm_objective_of_graded_distorted_box(const bool isVaried)
{
  BoundingBoxMesh bboxMesh(stk::topology::TETRAHEDRON_4, MPI_COMM_WORLD);
  bboxMesh.set_mesh_structure_type(BCC_BOUNDING_BOX_MESH);
  bboxMesh.set_reverse_entity_ids(isVaried);
  bboxMesh.set_alternate_processor_decomposition(isVaried);
  populate_bounding_box_mesh_and_activate(bboxMesh, stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.), 0.25);

  stk::mesh::BulkData & mesh = bboxMesh.bulk_data();
  const stk::mesh::MetaData & meta = bboxMesh.meta_data();
  const CoordinatesFieldRef coordsField(meta.coordinate_field(), 3);
  auto coordsData = coordsField.data<double, stk::mesh::ReadWrite>();
  for (auto * bucket : mesh.buckets(stk::topology::NODE_RANK))
    for (auto node : *bucket)
    {
      auto x = coordsData.entity_values(node);
      const stk::mesh::ComponentIdx ix(0), iy(1), iz(2);
      const stk::math::Vector3d orig(x(ix), x(iy), x(iz));
      const double distortion = 0.08 * std::sin(7.3 * orig[0] + 3.1 * orig[1] + 5.7 * orig[2]);
      x(ix) = orig[0] + distortion;
      x(iy) = orig[1] - distortion;
      x(iz) = std::copysign(std::pow(std::abs(orig[2]), 4.), orig[2]);
    }

  const MeanRatioNormObjective elementObjective;
  const MeshObjective objective(mesh, coordsField, AuxMetaData::get(meta).active_part(), elementObjective);
  const DistributedEntityVector & nodes = objective.get_nodes();
  const DistributedVector nodeCoords = gather_node_coordinates(mesh, coordsField, nodes);

  MeshObjectiveValueAndGradient result;
  result.value = objective.compute_value(nodeCoords);
  DistributedVector gradient;
  objective.fill_gradient(nodeCoords, gradient);

  std::vector<NodeCoordinatesAndGradient> owned;
  for (size_t i = 0; i < nodes.local_size(); ++i)
    owned.push_back({{nodeCoords[3 * i], nodeCoords[3 * i + 1], nodeCoords[3 * i + 2],
        gradient[3 * i], gradient[3 * i + 1], gradient[3 * i + 2]}});
  result.sortedNodesGradients = gather_sorted(owned);
  return result;
}

} // namespace

// The same value, and the same gradient at every node up to a float rounding step of the largest gradient:
// where a node's contributions span more than the double mantissa (as near-zero contributions of nearly
// ideal elements with larger ones), its sum depends on their order at the level of the smallest of them.
TEST(DistributionIndependentMeshObjective, gradedDistortedBox_sameValueAndGradientForReversedIdsAndOtherDecomposition)
{
  const MeshObjectiveValueAndGradient generated = mean_ratio_norm_objective_of_graded_distorted_box(false);
  const MeshObjectiveValueAndGradient varied = mean_ratio_norm_objective_of_graded_distorted_box(true);

  EXPECT_EQ(generated.value, varied.value);
  ASSERT_EQ(generated.sortedNodesGradients.size(), varied.sortedNodesGradients.size());
  double maxGradientComponent = 0.;
  for (auto & nodeGradient : generated.sortedNodesGradients)
    for (int k = 3; k < 6; ++k)
      maxGradientComponent = std::max(maxGradientComponent, std::abs(nodeGradient[k]));
  const double floatRoundingStep = maxGradientComponent * std::numeric_limits<float>::epsilon();
  size_t numDifferent = 0;
  for (size_t i = 0; i < generated.sortedNodesGradients.size(); ++i)
  {
    const auto & a = generated.sortedNodesGradients[i];
    const auto & b = varied.sortedNodesGradients[i];
    bool isSame = std::equal(a.begin(), a.begin() + 3, b.begin());
    for (int k = 3; isSame && k < 6; ++k)
      isSame = std::abs(a[k] - b[k]) <= floatRoundingStep;
    if (!isSame) ++numDifferent;
  }
  EXPECT_EQ(0u, numDifferent) << "of " << generated.sortedNodesGradients.size() << " nodes";
}

} // namespace krino
