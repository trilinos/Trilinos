#include <Akri_QuadRefiner.hpp>
#include <Akri_RefinerUtils.hpp>
#include <iostream>
#include <stdexcept>
#include <stk_topology/topology.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace krino {
namespace QuadRefiner {

unsigned determine_permutation_quad4(const unsigned caseId)
{
  STK_ThrowRequireMsg(caseId == 0 || caseId == 15, "Unfinished capability");
  static constexpr std::array<unsigned,16> permutations{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  return permutations[caseId];
}

unsigned determine_permuted_case_id_quad4(const unsigned caseId)
{
  STK_ThrowRequireMsg(caseId == 0 || caseId == 15, "Unfinished capability");
  static constexpr std::array<unsigned,16> permutedCaseIds{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15};
  return permutedCaseIds[caseId];
}

unsigned num_new_child_elements_quad4(const int caseId)
{
  switch(caseId)
  {
    case 0:
        return 0;
    case 15:
        return 4;
    default:
    {
      std::ostringstream errorMsg;
      errorMsg << "Case " << caseId << " non-uniform refinement not yet support for quads in num_new_child_elements_quad4.";
      throw std::runtime_error(errorMsg.str());
    }
  }
}

static std::array<unsigned,9> permutation_node_ordinals_quad4(const unsigned permutation)
{
  stk::topology topo = stk::topology::QUADRILATERAL_9_2D;
  std::array<unsigned,9> permutedNodes;
  topo.permutation_node_ordinals(permutation, permutedNodes.begin());
  return permutedNodes;
}

static std::array<unsigned,4> permutation_side_ordinals_quad4(const unsigned permutation)
{
  // nodes and sides permute the same way
  stk::topology topo = stk::topology::QUADRILATERAL_4_2D;
  std::array<unsigned,4> permutedSides;
  topo.permutation_node_ordinals(permutation, permutedSides.begin());
  return permutedSides;
}

std::vector<QuadDescription> refinement_child_nodes_and_sides_quad4(const unsigned caseId)
{
  std::vector<QuadDescription> childElems;

  const unsigned numChild = num_new_child_elements_quad4(caseId);
  childElems.reserve(numChild);

  const unsigned permutedCaseId = determine_permuted_case_id_quad4(caseId);
  const unsigned permutation = determine_permutation_quad4(caseId);
  const auto permutedParentNodeOrdinals = permutation_node_ordinals_quad4(permutation);
  const auto permutedParentSideOrdinals = permutation_side_ordinals_quad4(permutation);

  switch(permutedCaseId)
  {
    case 15:
    {
      static constexpr std::array<std::array<int,4>,4> childElemNodesFullyRefined{{ {{0,4,8,7}}, {{1,5,8,4}}, {{2,6,8,5}}, {{3,7,8,6}} }};
      static constexpr std::array<std::array<int,4>,4> childElemSidesFullyRefined{{ {{0,-1,-1,3}}, {{1,-1,-1,0}}, {{2,-1,-1,1}}, {{3,-1,-1,2}} }};
      append_child_elements(permutedParentNodeOrdinals, permutedParentSideOrdinals, childElemNodesFullyRefined, childElemSidesFullyRefined, childElems);
      break;
    }
    default:
    {
      std::ostringstream errorMsg;
      errorMsg << "Case " << permutedCaseId << " not supported in refinement_child_nodes_and_sides_quad4.";
      throw std::runtime_error(errorMsg.str());
    }
  }

  STK_ThrowRequireMsg(numChild == childElems.size(), "Mismatch of size " << numChild << "  " << childElems.size() << " for case " << caseId << " " << permutedCaseId);

  return childElems;
}

}}

