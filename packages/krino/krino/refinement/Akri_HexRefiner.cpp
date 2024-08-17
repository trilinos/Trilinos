#include <Akri_HexRefiner.hpp>
#include <Akri_RefinerUtils.hpp>
#include <iostream>
#include <stdexcept>
#include <stk_topology/topology.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace krino {
namespace HexRefiner {

unsigned determine_permutation_hex8(const unsigned caseId)
{
  STK_ThrowRequireMsg(caseId == 0 || caseId == 4095, "Unfinished capability");
  return 0;
}

unsigned determine_permuted_case_id_hex8(const unsigned caseId)
{
  STK_ThrowRequireMsg(caseId == 0 || caseId == 4095, "Unfinished capability");
  return caseId;
}

unsigned num_new_child_elements_hex8(const int caseId)
{
  switch(caseId)
  {
    case 0:
        return 0;
    case 4095:
        return 8;
    default:
    {
      std::ostringstream errorMsg;
      errorMsg << "Case " << caseId << " non-uniform refinement not yet support for quads in num_new_child_elements_hex8.";
      throw std::runtime_error(errorMsg.str());
    }
  }
}

static std::array<unsigned,27> permutation_node_ordinals_hex8(const unsigned permutation)
{
  stk::topology topo = stk::topology::HEXAHEDRON_27;
  std::array<unsigned,27> permutedNodes;
  topo.permutation_node_ordinals(permutation, permutedNodes.begin());
  return permutedNodes;
}

static std::array<unsigned,6> permutation_side_ordinals_hex8(const unsigned permutation)
{
  STK_ThrowRequireMsg(permutation == 0, "Unfinished capability");
  std::array<unsigned,6> permutedSides{0,1,2,3,4,5};
  return permutedSides;
}

std::vector<HexDescription> refinement_child_nodes_and_sides_hex8(const unsigned caseId)
{
  std::vector<HexDescription> childElems;

  const unsigned numChild = num_new_child_elements_hex8(caseId);
  childElems.reserve(numChild);

  const unsigned permutedCaseId = determine_permuted_case_id_hex8(caseId);
  const unsigned permutation = determine_permutation_hex8(caseId);
  const auto permutedParentNodeOrdinals = permutation_node_ordinals_hex8(permutation);
  const auto permutedParentSideOrdinals = permutation_side_ordinals_hex8(permutation);

  if (caseId == permutedCaseId)
  {
    static constexpr std::array<std::array<int,8>,8> childElemNodesFullyRefined{{
      {{0,8,21,11,12,25,20,23}}, {{8,1,9,21,25,13,24,20}}, {{21,9,2,10,20,24,14,26}}, {{11,21,10,3,23,20,26,15}},
      {{12,25,20,23,4,16,22,19}}, {{25,13,24,20,16,5,17,22}}, {{20,24,14,26,22,17,6,18}}, {{23,20,26,15,19,22,18,7}}
    }};
    static constexpr std::array<std::array<int,6>,8> childElemSidesFullyRefined{{
      {{0,-1,-1,3,4,-1}}, {{0,1,-1,-1,4,-1}}, {{-1,1,2,-1,4,-1}}, {{-1,-1,2,3,4,-1}},
      {{0,-1,-1,3,-1,5}}, {{0,1,-1,-1,-1,5}}, {{-1,1,2,-1,-1,5}}, {{-1,-1,2,3,-1,5}}
    }};
    append_child_elements(permutedParentNodeOrdinals, permutedParentSideOrdinals, childElemNodesFullyRefined, childElemSidesFullyRefined, childElems);
  }
  else
  {
    std::ostringstream errorMsg;
    errorMsg << "Case " << caseId << " not supported in refine_hex_8.";
    throw std::runtime_error(errorMsg.str());
  }

  STK_ThrowRequireMsg(numChild == childElems.size(), "Mismatch of size " << numChild << "  " << childElems.size() << " for case " << caseId << " " << permutedCaseId);

  return childElems;
}

}}
