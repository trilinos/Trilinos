// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_UNIT_TEST_UTILS_AKRI_UNIT_DECOMPOSITIONSIGNATURE_HPP_
#define KRINO_KRINO_UNIT_TEST_UTILS_AKRI_UNIT_DECOMPOSITIONSIGNATURE_HPP_

#include <array>
#include <string>
#include <vector>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace krino {

// An element as it appears in a decomposition: the names of its element parts (its phase block,
// and whether it is active, a parent, a child...) and its vertex coordinates in single precision,
// both sorted.  It does not depend on the ids of the element and its nodes, or on which processor
// owns them.
struct ElementSignature
{
  std::vector<std::string> partNames;
  std::vector<std::array<float, 3>> vertexCoordinates;
  bool operator<(const ElementSignature & other) const;
  bool operator==(const ElementSignature & other) const;
};

// The selected elements owned by each processor, as their signatures, sorted and gathered on
// processor 0 (empty on the others).  Two decompositions of the same mesh are the same exactly when
// their signatures are equal, whatever the ids and the parallel decomposition.
std::vector<ElementSignature> gather_decomposition_signature(
    const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector);

// The number of elements in either signature but not the other, on every processor (the signatures
// are compared on processor 0).
size_t global_num_elements_differing_between_signatures(stk::ParallelMachine comm,
    const std::vector<ElementSignature> & signature0,
    const std::vector<ElementSignature> & signature1);

// On processor 0, a description of up to maxNum elements in one signature but not the other (each
// with its parts, its vertices, and which signature it is in), for a test's failure message; empty
// on the other processors.
std::string describe_elements_differing_between_signatures(const std::vector<ElementSignature> & signature0,
    const std::vector<ElementSignature> & signature1,
    const size_t maxNum);

} // namespace krino

#endif /* KRINO_KRINO_UNIT_TEST_UTILS_AKRI_UNIT_DECOMPOSITIONSIGNATURE_HPP_ */
