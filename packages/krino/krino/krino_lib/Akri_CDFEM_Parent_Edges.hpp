// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_CDFEM_PARENT_EDGES_H_
#define KRINO_INCLUDE_AKRI_CDFEM_PARENT_EDGES_H_
#include <Akri_CDFEM_Parent_Edge.hpp>
#include <Akri_NodeToCapturedDomains.hpp>
#include <Akri_OrderedIdPair.hpp>

namespace krino
{
class CDFEM_Support;
class Phase_Support;
class InterfaceID;
class ParentsToChildMapper;
struct LS_Field;

typedef OrderedIdPair ParentEdgeKey;
typedef std::map<ParentEdgeKey,CDFEM_Parent_Edge> ParentEdgeMap;

std::set<int> get_phases_present_on_edges(const std::vector<const CDFEM_Parent_Edge *> & parentEdges);
std::set<int> get_phases_present_on_edges_and_interior(const std::vector<const CDFEM_Parent_Edge *> & elementParentEdges);
std::set<InterfaceID> get_interfaces_present(const std::vector<const CDFEM_Parent_Edge *> & elementParentEdges);

const CDFEM_Parent_Edge * find_parent_edge(const stk::mesh::BulkData & mesh, const ParentEdgeMap & parentEdges, stk::mesh::Entity node0, stk::mesh::Entity node1);

ParentEdgeMap
build_parent_edges_using_elements(const stk::mesh::BulkData & mesh,
    const ParentsToChildMapper & parentsToChildMapper,
    const bool shouldLinearizeEdges,
    const std::vector<stk::mesh::Entity> & elements,
    const stk::mesh::Selector & parentElementSelector,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields);

ParentEdgeMap
build_parent_edges(const stk::mesh::BulkData & mesh,
    const ParentsToChildMapper & parentsToChildMapper,
    const bool shouldLinearizeEdges,
    const stk::mesh::Selector & parentElementSelector,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields);

void fill_element_parent_edges(const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity elem,
    const ParentEdgeMap & parentEdges,
    std::vector<const CDFEM_Parent_Edge *> & elementParentEdges,
    std::vector<bool> & areParentEdgesAreOrientedSameAsElementEdges);

void fill_face_nodes_and_parent_edges(const stk::topology & elementTopology,
    const int iFace,
    const std::vector<stk::mesh::Entity> & elementNodes,
    const std::vector<const CDFEM_Parent_Edge *> & elementParentEdges,
    const std::vector<bool> & areParentEdgesOrientedSameAsElementEdges,
    std::vector<stk::mesh::Entity> & faceNodes,
    std::vector<const CDFEM_Parent_Edge *> & faceParentEdges,
    std::vector<bool> & areParentEdgesOrientedSameAsFaceEdges);

stk::mesh::Selector get_decomposed_cdfem_parent_element_selector(const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport);

stk::mesh::Selector get_potential_cdfem_parent_element_selector(const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport);

std::vector<stk::mesh::Entity> get_owned_parent_elements(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & parentElementSelector);

}




#endif /* KRINO_INCLUDE_AKRI_CDFEM_PARENT_EDGES_H_ */
