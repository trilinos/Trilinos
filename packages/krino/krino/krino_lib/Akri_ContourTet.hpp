// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_CONTOURTET_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_CONTOURTET_HPP_
#include <array>

#include <stk_math/StkVector.hpp>

namespace krino {
class FacetedSurfaceBase;

class ContourTet
{
public:
  static unsigned get_permutation_for_case(const unsigned caseId);
  static std::array<unsigned, 4> get_permuted_tet4_node_ordinals(const unsigned caseId);
  static std::array<unsigned, 4> get_permuted_tet4_node_ordinals_for_permutation(const unsigned permutation);
  static std::array<unsigned, 10> get_permuted_tet10_node_ordinals(const unsigned caseId);
  static std::array<unsigned, 10> get_permuted_tet10_node_ordinals_for_permutation(const unsigned permutation);
  static const std::array<unsigned, 4> & get_permuted_side_ordinals(const unsigned caseId);
  static const std::array<unsigned, 4> & get_permuted_side_ordinals_for_permutation(const unsigned permutation);
  static unsigned get_permuted_case_id(const unsigned caseId);
  static unsigned compute_case_id(const std::array<int,4> & nodeSigns);

  static void append_facets_for_tet4(const std::array<stk::math::Vector3d,4> & coords,
    const std::array<double,4> & dist,
    FacetedSurfaceBase & facets);
  static void append_facets_for_tet4_with_tet10_distance(const std::array<stk::math::Vector3d,4> & coords,
    const std::array<double,10> & tet10Dist,
    FacetedSurfaceBase & facets);
};

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_CONTOURTET_HPP_ */
