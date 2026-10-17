// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_CONTOURTRI_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_CONTOURTRI_HPP_

#include <array>

#include <stk_math/StkVector.hpp>

namespace krino {
class FacetedSurfaceBase;

class ContourTri
{
public:
  static unsigned get_permutation_for_case(const unsigned caseId);
  static std::array<unsigned, 3> get_permuted_tri3_node_ordinals(const unsigned caseId);
  static std::array<unsigned, 6> get_permuted_tri6_node_ordinals(const unsigned caseId);
  static unsigned get_permuted_case_id(const unsigned caseId);
  static unsigned compute_case_id(const std::array<int,3> & nodeSigns);

  static void append_facets_for_tri3(const std::array<stk::math::Vector3d,3> & coords,
    const std::array<double,3> & dist,
    FacetedSurfaceBase & facets);
  static void append_facets_for_tri3_with_tri6_distance(const std::array<stk::math::Vector3d,3> & coords,
    const std::array<double,6> & tri6Dist,
    FacetedSurfaceBase & facets);
};

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_CONTOURTRI_HPP_ */
