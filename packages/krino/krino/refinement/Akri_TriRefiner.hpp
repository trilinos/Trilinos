/*
 * Akri_TriRefiner.hpp
 *
 *  Created on: Oct 19, 2022
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_TRIREFINER_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_TRIREFINER_HPP_
#include <array>
#include <vector>
#include <stk_math/StkVector.hpp>

namespace krino {
namespace TriRefiner {

struct TriDescription
{
  std::array<int, 3> nodeIds;
  std::array<int, 3> sideIds;
};

std::vector<TriDescription> refinement_child_nodes_and_sides_tri3(const unsigned encodedEdgesToRefine, const std::array<stk::math::Vector3d,6> & elementNodeCoords, const std::array<int,3> & elementNodeScore);
unsigned determine_permutation_tri3(const unsigned caseId);
unsigned determine_permuted_case_id_tri3(const unsigned caseId);
unsigned num_new_child_elements_tri3(const int caseId);

}
}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_TRIREFINER_HPP_ */
