#ifndef KRINO_KRINO_KRINO_LIB_AKRI_CONTOURUTILS_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_CONTOURUTILS_HPP_
#include <array>

namespace krino {

class ContourTri
{
public:
  static unsigned get_permutation_for_case(const unsigned caseId);
  static std::array<unsigned, 6> get_permuted_node_ordinals(const unsigned caseId);
  static unsigned get_permuted_case_id(const unsigned caseId);
  static unsigned compute_case_id(const std::array<int,3> & nodeSigns);
};

class ContourTet
{
public:
  static unsigned get_permutation_for_case(const unsigned caseId);
  static std::array<unsigned, 10> get_permuted_node_ordinals(const unsigned caseId);
  static std::array<unsigned, 10> get_permuted_node_ordinals_for_permutation(const unsigned permutation);
  static const std::array<unsigned, 4> & get_permuted_side_ordinals(const unsigned caseId);
  static const std::array<unsigned, 4> & get_permuted_side_ordinals_for_permutation(const unsigned permutation);
  static unsigned get_permuted_case_id(const unsigned caseId);
  static unsigned compute_case_id(const std::array<int,4> & nodeSigns);
};

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_CONTOURUTILS_HPP_ */
