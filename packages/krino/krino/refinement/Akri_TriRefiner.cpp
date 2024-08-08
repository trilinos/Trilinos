/*
 * Akri_TriRefiner.cpp
 *
 *  Created on: Oct 19, 2022
 *      Author: drnoble
 */
#include "Akri_TriRefiner.hpp"

#include <array>

#include <stk_topology/topology_decl.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_RefinerUtils.hpp>

namespace krino {
namespace TriRefiner {

unsigned determine_permutation_tri3(const unsigned caseId)
{
  static constexpr std::array<unsigned,8> permutations{0, 0, 2, 0, 1, 1, 2, 0};
  return permutations[caseId];
}

unsigned determine_permuted_case_id_tri3(const unsigned caseId)
{
  static constexpr std::array<unsigned,8> permutedCaseIds{0, 1, 1, 3, 1, 3, 3, 7};
  return permutedCaseIds[caseId];
}

unsigned num_new_child_elements_tri3(const int caseId)
{
  const unsigned permutedCaseId = determine_permuted_case_id_tri3(caseId);
  switch(permutedCaseId)
  {
    case 0:
        return 0;
    case 1:
        return 2;
    case 3:
        return 3;
    case 7:
        return 4;
    default:
    {
      std::ostringstream errorMsg;
      errorMsg << "Case " << caseId << " not supported in determine_permuted_case_id_tri3.";
      throw std::runtime_error(errorMsg.str());
    }
  }
}

static std::array<stk::math::Vector3d,6> calculate_refined_tri_coordinates(const std::array<stk::math::Vector3d,3> & parentElementNodeCoords, const std::array<unsigned,6> & permutedParentNodeOrdinals)
{
  std::array<stk::math::Vector3d,6> refinedTriNodeCoords;
  refinedTriNodeCoords[0] = parentElementNodeCoords[permutedParentNodeOrdinals[0]];
  refinedTriNodeCoords[1] = parentElementNodeCoords[permutedParentNodeOrdinals[1]];
  refinedTriNodeCoords[2] = parentElementNodeCoords[permutedParentNodeOrdinals[2]];
  refinedTriNodeCoords[3] = 0.5*(refinedTriNodeCoords[0]+refinedTriNodeCoords[1]);
  refinedTriNodeCoords[4] = 0.5*(refinedTriNodeCoords[1]+refinedTriNodeCoords[2]);
  refinedTriNodeCoords[5] = 0.5*(refinedTriNodeCoords[2]+refinedTriNodeCoords[0]);
  return refinedTriNodeCoords;
}

template<size_t NUMPARENTNODES, size_t NUMCHILDNODES>
static double compute_quality_of_child_elem(const QualityMetric & qualityMetric, const std::array<stk::math::Vector3d,NUMPARENTNODES> & parentNodeCoords, const std::array<int,NUMCHILDNODES> & childElemNodeIndices, std::vector<stk::math::Vector3d> & prealllocatedScratchElemNodeCoords)
{
  for (size_t i=0; i<NUMCHILDNODES; ++i)
    prealllocatedScratchElemNodeCoords[i] = parentNodeCoords[childElemNodeIndices[i]];
  return qualityMetric.get_element_quality_metric(prealllocatedScratchElemNodeCoords);
}

template<size_t NUMPARENTNODES, size_t NUMCHILDNODES, size_t NUMCHILDRENINCONFIGURATION>
double compute_quality_of_child_element_configuration(const QualityMetric & qualityMetric,
    const std::array<stk::math::Vector3d,NUMPARENTNODES> & parentNodeCoords,
    const std::array<std::array<int,NUMCHILDNODES>,NUMCHILDRENINCONFIGURATION> & elementConfigurations,
    std::vector<stk::math::Vector3d> & prealllocatedScratchElemNodeCoords)
{
  double worstQuality = qualityMetric.get_best_value_for_metric();
  for (auto && elementConfiguration : elementConfigurations)
  {
    const double elementQuality = compute_quality_of_child_elem(qualityMetric, parentNodeCoords, elementConfiguration, prealllocatedScratchElemNodeCoords);
    if (qualityMetric.is_first_quality_metric_better_than_second(worstQuality, elementQuality))
      worstQuality = elementQuality;
  }
  return worstQuality;
}

template<size_t NUMPARENTNODES, size_t NUMCHILDNODES, size_t NUMCHILDRENINCONFIGURATION>
double compute_quality_of_child_element_configuration_terminating_early_if_below_threshold(const QualityMetric & qualityMetric,
    const std::array<stk::math::Vector3d,NUMPARENTNODES> & parentNodeCoords,
    const std::array<std::array<int,NUMCHILDNODES>,NUMCHILDRENINCONFIGURATION> & elementConfigurations,
    std::vector<stk::math::Vector3d> & prealllocatedScratchElemNodeCoords,
    const double qualityThreshold)
{
  double worstQuality = qualityMetric.get_best_value_for_metric();
  for (auto && elementConfiguration : elementConfigurations)
  {
    const double elementQuality = compute_quality_of_child_elem(qualityMetric, parentNodeCoords, elementConfiguration, prealllocatedScratchElemNodeCoords);
    if (qualityMetric.is_first_quality_metric_better_than_second(worstQuality, elementQuality))
    {
      worstQuality = elementQuality;
      if (qualityMetric.is_first_quality_metric_better_than_second(qualityThreshold, worstQuality))
        break;
    }
  }
  return worstQuality;
}

static int which_configuration_of_child_tri_elements_is_best(const std::array<stk::math::Vector3d,3> & parentElementNodeCoords,
    const std::array<int,3> & elementNodeScore,
    const std::array<unsigned,6> & permutedParentNodeOrdinals,
    const std::array<std::array<std::array<int,3>,2>,2> & configs)
{
  const auto refinedTriNodeCoords = calculate_refined_tri_coordinates(parentElementNodeCoords, permutedParentNodeOrdinals);

  const ScaledJacobianQualityMetric qualityMetric;
  std::vector<stk::math::Vector3d> prealllocatedScratchTriNodeCoords(3);
  const double qual0 = compute_quality_of_child_element_configuration(qualityMetric, refinedTriNodeCoords, configs[0], prealllocatedScratchTriNodeCoords);
  const double qual1 = compute_quality_of_child_element_configuration_terminating_early_if_below_threshold(qualityMetric, refinedTriNodeCoords, configs[1], prealllocatedScratchTriNodeCoords, qual0);

  if (qualityMetric.is_first_quality_metric_better_than_second(qual0, qual1))
    return 0;
  else if (qualityMetric.is_first_quality_metric_better_than_second(qual1, qual0))
    return 1;
  else
    return elementNodeScore[permutedParentNodeOrdinals[0]] < elementNodeScore[permutedParentNodeOrdinals[1]];
}

static std::array<unsigned,6> permutation_node_ordinals_tri3(const unsigned caseId)
{
  stk::topology topo = stk::topology::TRIANGLE_6_2D;
  std::array<unsigned,6> permutation;
  topo.permutation_node_ordinals(determine_permutation_tri3(caseId), permutation.begin());
  return permutation;
}

static std::array<unsigned,3> permutation_side_ordinals_tri3(const unsigned caseId)
{
  // nodes and sides permute the same way
  stk::topology topo = stk::topology::TRIANGLE_3_2D;
  std::array<unsigned,3> permutation;
  topo.permutation_node_ordinals(determine_permutation_tri3(caseId), permutation.begin());
  return permutation;
}

std::vector<TriDescription> refinement_child_nodes_and_sides_tri3(const unsigned caseId, const std::array<stk::math::Vector3d,3> & elementNodeCoords, const std::array<int,3> & elementNodeScore)
{
  std::vector<TriDescription> childElemNodes;

  const unsigned numChild = num_new_child_elements_tri3(caseId);
  childElemNodes.reserve(numChild);

  const unsigned permutedCaseId = determine_permuted_case_id_tri3(caseId);
  const auto permutedParentNodeOrdinals = permutation_node_ordinals_tri3(caseId);
  const auto permutedParentSideOrdinals = permutation_side_ordinals_tri3(caseId);

  switch(permutedCaseId)
  {
    case 1:
    {
      static constexpr std::array<std::array<int,3>,2> childElemNodesSplit{{ {{0,3,2}}, {{1,2,3}} }};
      static constexpr std::array<std::array<int,3>,2> childElemSidesSplit{{ {{0,-1,2}}, {{1,-1,0}} }};
      append_child_elements(permutedParentNodeOrdinals, permutedParentSideOrdinals, childElemNodesSplit, childElemSidesSplit, childElemNodes);
      break;
    }
    case 3:
    {
      static constexpr std::array<std::array<int,3>,1> triChildNodes{{ {{1,4,3}} }};
      static constexpr std::array<std::array<int,3>,1> triChildSides{{ {{1,-1,0}} }};
      static constexpr std::array<std::array<std::array<int,3>,2>,2> quadChildNodesConfigs
        {{ {{ {{0,3,2}}, {{2,3,4}} }},
           {{ {{0,3,4}}, {{2,0,4}} }} }};
      static constexpr std::array<std::array<std::array<int,3>,2>,2> quadChildSideConfigs
        {{ {{ {{0,-1,2}}, {{-1,-1,1}} }},
           {{ {{0,-1,-1}}, {{2,-1,1}} }} }};
      const int bestConfig = which_configuration_of_child_tri_elements_is_best(elementNodeCoords, elementNodeScore, permutedParentNodeOrdinals, quadChildNodesConfigs);
      append_child_elements(permutedParentNodeOrdinals, permutedParentSideOrdinals, triChildNodes, triChildSides, childElemNodes);
      append_child_elements(permutedParentNodeOrdinals, permutedParentSideOrdinals, quadChildNodesConfigs[bestConfig], quadChildSideConfigs[bestConfig], childElemNodes);
      break;
    }
    case 7:
    {
      static constexpr std::array<std::array<int,3>,4> childElemNodesFullyRefined{{ {{0,3,5}}, {{1,4,3}}, {{2,5,4}}, {{3,4,5}} }};
      static constexpr std::array<std::array<int,3>,4> childElemSidesFullyRefined{{ {{0,-1,2}}, {{1,-1,0}}, {{2,-1,1}}, {{-1,-1,-1}} }};
      append_child_elements(permutedParentNodeOrdinals, permutedParentSideOrdinals, childElemNodesFullyRefined, childElemSidesFullyRefined, childElemNodes);
      break;
    }
    default:
    {
      std::ostringstream errorMsg;
      errorMsg << "Case " << caseId << " not supported in refine_tri_3.";
      throw std::runtime_error(errorMsg.str());
    }
  }

  STK_ThrowRequireMsg(numChild == childElemNodes.size(), "Mismatch of size " << numChild << "  " << childElemNodes.size() << " for case " << caseId << " " << permutedCaseId);

  return childElemNodes;
}

} // namespace TriRefiner
} // namespace krino
