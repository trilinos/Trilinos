#ifndef KRINO_KRINO_REFINEMENT_AKRI_HEXREFINER_HPP_
#define KRINO_KRINO_REFINEMENT_AKRI_HEXREFINER_HPP_

#include <array>
#include <vector>

namespace krino {
namespace HexRefiner {

struct HexDescription
{
  std::array<int, 8> nodeIds;
  std::array<int, 6> sideIds;
};

unsigned determine_permutation_hex8(const unsigned caseId);
unsigned determine_permuted_case_id_hex8(const unsigned caseId);
unsigned num_new_child_elements_hex8(const int caseId);
std::vector<HexDescription> refinement_child_nodes_and_sides_hex8(const unsigned caseId);

}
}



#endif /* KRINO_KRINO_REFINEMENT_AKRI_HEXREFINER_HPP_ */
