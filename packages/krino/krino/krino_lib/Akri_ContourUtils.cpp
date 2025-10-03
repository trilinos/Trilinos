/*
 * Akri_ContourUtils.cpp
 *
 *  Created on: Oct 1, 2025
 *      Author: drnoble
 */
#include <Akri_ContourUtils.hpp>

#include <stk_topology/topology.hpp>

namespace krino {

unsigned ContourTri::compute_case_id(const std::array<int,3> & nodeSigns)
{
  return (nodeSigns[0]+1) +
      (nodeSigns[1]+1)*3 +
      (nodeSigns[2]+1)*9;
}

unsigned ContourTri::get_permutation_for_case(const unsigned caseId)
{
  static constexpr std::array<unsigned,27> casePermutations =
      { 0, 0, 0, 2, 0, 0, 2, 1, 1, 1, //  0-9
        1, 2, 2, 0, 2, 2, 1, 1, 1, 1, //  10-19
        2, 0, 0, 2, 0, 0, 0 };        //  20-26
  return casePermutations[caseId];
}

std::array<unsigned, 6> ContourTri::get_permuted_node_ordinals(const unsigned caseId)
{
  stk::topology topo = stk::topology::TRIANGLE_6_2D;
  std::array<unsigned,6> permutedNodeOrdinals;
  topo.permutation_node_ordinals(get_permutation_for_case(caseId), permutedNodeOrdinals.begin());
  return permutedNodeOrdinals;
}

unsigned ContourTri::get_permuted_case_id(const unsigned caseId)
{
  static constexpr std::array<unsigned,27> permutedCaseIds =
      { 0, 1, 2, 1, 4, 5, 2,21,24, 1, //  0-9
        4,21, 4,13,22, 5,22,25, 2, 5, //  10-19
       24,21,22,25,24,25,26 };        //  20-26

  return permutedCaseIds[caseId];
}

unsigned ContourTet::compute_case_id(const std::array<int,4> & nodeSigns)
{
  return (nodeSigns[0]+1) +
      (nodeSigns[1]+1)*3 +
      (nodeSigns[2]+1)*9 +
      (nodeSigns[3]+1)*27;
}

unsigned ContourTet::get_permutation_for_case(const unsigned caseId)
{
  static constexpr std::array<unsigned,81> casePermutations =
      { 0, 0, 0, 1, 0, 0, 1, 5, 0, 2, //  0-9
        2, 6, 1, 0, 0, 1, 1,10, 2, 2, //  10-19
        2,11, 2, 4, 1, 8, 4, 4, 3, 3, //  20-29
        4, 3, 3, 9, 5, 7, 7, 6, 6, 9, //  30-39
        0, 9, 9, 6, 7, 7, 7, 9,11, 3, //  40-49
        4, 3, 3, 4, 4, 8, 3, 4, 4,11, //  50-59
        4, 2, 2,10, 8, 1,10, 0, 1, 6, //  60-69
        2, 2, 7, 5, 1, 0, 0, 1, 0, 0, //  70-79
        0 };                          //  80
  return casePermutations[caseId];
}

std::array<unsigned, 10> ContourTet::get_permuted_node_ordinals_for_permutation(const unsigned permutation)
{
  stk::topology topo = stk::topology::TETRAHEDRON_10;
  std::array<unsigned,10> permutedNodeOrdinals;
  topo.permutation_node_ordinals(permutation, permutedNodeOrdinals.begin());
  return permutedNodeOrdinals;
}

std::array<unsigned, 10> ContourTet::get_permuted_node_ordinals(const unsigned caseId)
{
  return get_permuted_node_ordinals_for_permutation(get_permutation_for_case(caseId));
}

unsigned ContourTet::get_permuted_case_id(const unsigned caseId)
{
  static constexpr std::array<unsigned,81> permutedCaseIds =
      { 0, 1, 2, 1, 4, 5, 2, 5, 8, 1, //  0-9
        4, 5, 4,13,14, 5,14,75, 2, 5, //  10-19
        8, 5,14,75, 8,75,78, 1, 4, 5, //  20-29
        4,13,14, 5,14,75, 4,13,14,13, //  30-39
       40,67,14,67,76, 5,14,75,14,67, //  40-49
       76,75,76,79, 2, 5, 8, 5,14,75, //  50-59
        8,75,78, 5,14,75,14,67,76,75, //  60-69
       76,79, 8,75,78,75,76,79,78,79, //  70-79
       80 };                          //  80

  return permutedCaseIds[caseId];
}

const std::array<unsigned, 4> & ContourTet::get_permuted_side_ordinals_for_permutation(const unsigned permutation)
{
  static constexpr std::array<std::array<unsigned,4>,12> permutationSideOrdinals =
      {{ {{ 0, 1, 2, 3 }},
         {{ 1, 2, 0, 3 }},
         {{ 2, 0, 1, 3 }},
         {{ 2, 1, 3, 0 }},
         {{ 1, 3, 2, 0 }},
         {{ 3, 2, 1, 0 }},
         {{ 3, 1, 0, 2 }},
         {{ 1, 0, 3, 2 }},
         {{ 0, 3, 1, 2 }},
         {{ 0, 2, 3, 1 }},
         {{ 2, 3, 0, 1 }},
         {{ 3, 0, 2, 1 }}
      }};
  return permutationSideOrdinals[permutation];
}

const std::array<unsigned, 4> & ContourTet::get_permuted_side_ordinals(const unsigned caseId)
{
  return get_permuted_side_ordinals_for_permutation(get_permutation_for_case(caseId));
}

}


