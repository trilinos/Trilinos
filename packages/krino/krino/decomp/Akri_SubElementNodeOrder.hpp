// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_DECOMP_AKRI_SUBELEMENTNODEORDER_HPP_
#define KRINO_KRINO_DECOMP_AKRI_SUBELEMENTNODEORDER_HPP_

#include <Akri_EntityOrder.hpp>

namespace krino {

class SubElementNode;
class SubElementNodeAncestry;

// The orders of Akri_EntityOrder.hpp for sub-element nodes: a node by its location (coordinates) and
// its entity id, and an ancestry (SubElementNodeAncestry::compare) with mesh nodes by that order and
// child nodes by their parents' ancestries, with the same order at every level.

bool less_by_coordinates_because_it_impacts_structure(const SubElementNode & a, const SubElementNode & b);
bool less_conditionally_by_coordinates_because_it_impacts_structure(
    const bool globalIdsAreParallelConsistent, const SubElementNode & a, const SubElementNode & b);
bool less_by_global_id_because_it_doesnt_impact_structure(const SubElementNode & a, const SubElementNode & b);
bool less_by_global_id_because_it_is_requested(const SubElementNode & a, const SubElementNode & b);

bool less_by_coordinates_because_it_impacts_structure(
    const SubElementNodeAncestry & a, const SubElementNodeAncestry & b);
bool less_conditionally_by_coordinates_because_it_impacts_structure(const bool globalIdsAreParallelConsistent,
    const SubElementNodeAncestry & a,
    const SubElementNodeAncestry & b);
bool less_by_global_id_because_it_doesnt_impact_structure(
    const SubElementNodeAncestry & a, const SubElementNodeAncestry & b);
bool less_by_global_id_because_it_is_requested(const SubElementNodeAncestry & a, const SubElementNodeAncestry & b);

} // namespace krino

#endif /* KRINO_KRINO_DECOMP_AKRI_SUBELEMENTNODEORDER_HPP_ */
