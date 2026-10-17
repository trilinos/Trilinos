// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_RefinementMarkerUtils.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_RefinementManager.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <algorithm>
#include <array>
#include <limits>

namespace krino {

void initialize_marker(const stk::mesh::BulkData& /*mesh*/,
      const RefinementManager & refinement,
      const bool isDefaultCoarsen)
{
  const FieldRef elementMarkerField = refinement.get_marker_field_and_sync_to_host();
  const int initialVal = isDefaultCoarsen ? static_cast<int>(Refinement_Marker::COARSEN) : static_cast<int>(Refinement_Marker::NOTHING);
  stk::mesh::field_fill(initialVal, elementMarkerField);
}

Refinement_Marker determine_refinement_marker(const bool isElementSelectedForRefinement, const int minRefineLevel, const int elementRefineLevel, const bool isDefaultCoarsen)
{
  auto marker = isDefaultCoarsen ? Refinement_Marker::COARSEN : Refinement_Marker::NOTHING;
  if (isElementSelectedForRefinement)
  {
    if (elementRefineLevel < minRefineLevel)
      marker = Refinement_Marker::REFINE;
    else
      marker = Refinement_Marker::NOTHING;
  }

  return marker;
}

void mark_elements_to_reach_refinement_level(const stk::mesh::BulkData& /*mesh*/,
      const RefinementManager & refinement,
      const std::vector<stk::mesh::Entity> elementsToMark,
      const int minRefineLevel,
      const bool isDefaultCoarsen)
{
  const FieldRef elementMarkerField = refinement.get_marker_field_and_sync_to_host();
  constexpr bool doMarkElement = true;

  auto elementMarkerData = elementMarkerField.data<int, stk::mesh::ReadWrite>();
  for( auto&& elem : elementsToMark )
  {
    auto marker = elementMarkerData.entity_values(elem);
    const int elementRefineLevel = refinement.fully_refined_level(elem);

    marker() = static_cast<int>(determine_refinement_marker(doMarkElement, minRefineLevel, elementRefineLevel, isDefaultCoarsen));
  }
}

bool should_continue_to_coarsen(const int refinementIterCount, const int targetRefineLevel)
{
  // Stop coarsening after num_levels of refinement to avoid infinite looping
  // from the interface position moving between elements because of snapping changes
  // with refinement
  return refinementIterCount < targetRefineLevel;
}

static double compute_edge_length(const FieldRef coordsField, const stk::topology elemTopology, const stk::mesh::Entity * const elemNodes, const unsigned iEdge)
{
  std::array<stk::mesh::Entity,3> edgeNodes;

  elemTopology.edge_nodes(elemNodes, iEdge, edgeNodes.data());

  auto coordsData = coordsField.data<double>();
  const stk::math::Vector3d edge_node1_coords(coordsData.entity_values(edgeNodes[0]), elemTopology.dimension());
  const stk::math::Vector3d edge_node2_coords(coordsData.entity_values(edgeNodes[1]), elemTopology.dimension());
  return (edge_node2_coords-edge_node1_coords).length();
}

void write_refinement_level_sizes(const stk::mesh::BulkData& mesh,
    const RefinementManager & refinement,
    const FieldRef coordsField,
    const std::vector<stk::mesh::Entity> & elements,
    const int maxRefineLevel)
{
  const int maxNumRefinementLevels = maxRefineLevel+1;
  std::vector<double> min_edge_lengths(maxNumRefinementLevels, std::numeric_limits<double>::max());
  std::vector<double> max_edge_lengths(maxNumRefinementLevels, 0.0);

  for( auto&& elem : elements )
  {
    if (!refinement.is_transition(elem))
    {
      const int elementRefineLevel = refinement.fully_refined_level(elem);
      STK_ThrowRequire(elementRefineLevel < maxNumRefinementLevels);
      const stk::topology elemTopology = mesh.bucket(elem).topology();
      const stk::mesh::Entity * const elemNodes = mesh.begin(elem, stk::topology::NODE_RANK);

      for (unsigned iEdge = 0; iEdge < elemTopology.num_edges(); ++iEdge)
      {
        const double length = compute_edge_length(coordsField, elemTopology, elemNodes, iEdge);

        min_edge_lengths[elementRefineLevel] = std::min(min_edge_lengths[elementRefineLevel], length);
        max_edge_lengths[elementRefineLevel] = std::max(max_edge_lengths[elementRefineLevel], length);
      }
    }
  }

  std::vector<double> localVals = min_edge_lengths;
  stk::all_reduce_min(mesh.parallel(), localVals.data(), min_edge_lengths.data(), min_edge_lengths.size());
  localVals = max_edge_lengths;
  stk::all_reduce_max(mesh.parallel(), localVals.data(), max_edge_lengths.data(), max_edge_lengths.size());

  for (int i=0; i<maxNumRefinementLevels; ++i)
  {
    krinolog << "Min and Max sizes for refinement level " << i <<  " " << min_edge_lengths[i] << " " << max_edge_lengths[i] << stk::diag::dendl;
  }
}

} // namespace krino
