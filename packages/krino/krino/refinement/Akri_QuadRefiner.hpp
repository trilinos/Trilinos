#ifndef KRINO_KRINO_REFINEMENT_AKRI_QUADREFINER_HPP_
#define KRINO_KRINO_REFINEMENT_AKRI_QUADREFINER_HPP_

#include <array>
#include <vector>

namespace krino {
namespace QuadRefiner {

struct QuadDescription
{
  std::array<int, 4> nodeIds;
  std::array<int, 4> sideIds;
};

unsigned determine_permutation_quad4(const unsigned caseId);
unsigned determine_permuted_case_id_quad4(const unsigned caseId);
unsigned num_new_child_elements_quad4(const int caseId);
std::vector<QuadDescription> refinement_child_nodes_and_sides_quad4(const unsigned caseId);

}
}



#endif /* KRINO_KRINO_REFINEMENT_AKRI_QUADREFINER_HPP_ */
