// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_SubElementNodeOrder.hpp>

#include <Akri_SubElement.hpp>
#include <Akri_SubElementNodeAncestry.hpp>

namespace krino {

bool less_by_coordinates_because_it_impacts_structure(const SubElementNode & a, const SubElementNode & b)
{
  return less_by_coordinates_because_it_impacts_structure(a.coordinates(), a.entityId(), b.coordinates(), b.entityId());
}

bool less_conditionally_by_coordinates_because_it_impacts_structure(
    const bool globalIdsAreParallelConsistent, const SubElementNode & a, const SubElementNode & b)
{
  return less_conditionally_by_coordinates_because_it_impacts_structure(
      globalIdsAreParallelConsistent, a.coordinates(), a.entityId(), b.coordinates(), b.entityId());
}

bool less_by_global_id_because_it_doesnt_impact_structure(const SubElementNode & a, const SubElementNode & b)
{
  return less_by_global_id_because_it_doesnt_impact_structure(a.entityId(), b.entityId());
}

bool less_by_global_id_because_it_is_requested(const SubElementNode & a, const SubElementNode & b)
{
  return less_by_global_id_because_it_is_requested(a.entityId(), b.entityId());
}

bool less_by_coordinates_because_it_impacts_structure(
    const SubElementNodeAncestry & a, const SubElementNodeAncestry & b)
{
  return SubElementNodeAncestry::compare(a,
      b,
      [](const SubElementNode & nodeA, const SubElementNode & nodeB)
      { return less_by_coordinates_because_it_impacts_structure(nodeA, nodeB); });
}

bool less_conditionally_by_coordinates_because_it_impacts_structure(const bool globalIdsAreParallelConsistent,
    const SubElementNodeAncestry & a,
    const SubElementNodeAncestry & b)
{
  return SubElementNodeAncestry::compare(a,
      b,
      [globalIdsAreParallelConsistent](const SubElementNode & nodeA, const SubElementNode & nodeB)
      { return less_conditionally_by_coordinates_because_it_impacts_structure(globalIdsAreParallelConsistent, nodeA, nodeB); });
}

bool less_by_global_id_because_it_doesnt_impact_structure(
    const SubElementNodeAncestry & a, const SubElementNodeAncestry & b)
{
  return SubElementNodeAncestry::compare(a,
      b,
      [](const SubElementNode & nodeA, const SubElementNode & nodeB)
      { return less_by_global_id_because_it_doesnt_impact_structure(nodeA, nodeB); });
}

bool less_by_global_id_because_it_is_requested(const SubElementNodeAncestry & a, const SubElementNodeAncestry & b)
{
  return SubElementNodeAncestry::compare(a,
      b,
      [](const SubElementNode & nodeA, const SubElementNode & nodeB)
      { return less_by_global_id_because_it_is_requested(nodeA, nodeB); });
}

} // namespace krino
