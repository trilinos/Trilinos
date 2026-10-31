// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_REFINEMENT_DRIVERS_AKRI_REFINEMENTTREEELEMENTPARTS_HPP_
#define KRINO_KRINO_REFINEMENT_DRIVERS_AKRI_REFINEMENTTREEELEMENTPARTS_HPP_

#include <stk_mesh/base/Types.hpp>

#include <functional>
#include <set>
#include <utility>
#include <vector>

namespace krino {

class RefinementManager;

// The element part of each locally owned refinement parent element, determined from the element parts of its
// children by determineParentElementPartOrdinal (given the set of the children's element part ordinals), from
// the leaf children up the refinement tree.  Children owned by other processors are communicated to the
// processor that owns their parent.
std::vector<std::pair<stk::mesh::Entity,unsigned>> get_owned_refinement_parents_and_their_element_part(const stk::mesh::BulkData& mesh,
    const RefinementManager & refinement,
    const std::function<unsigned(const std::set<unsigned> &)> & determineParentElementPartOrdinal);

} // namespace krino

#endif /* KRINO_KRINO_REFINEMENT_DRIVERS_AKRI_REFINEMENTTREEELEMENTPARTS_HPP_ */
