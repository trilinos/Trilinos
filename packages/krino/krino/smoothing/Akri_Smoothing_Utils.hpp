// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SMOOTHING_AKRI_SMOOTHING_UTILS_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_SMOOTHING_UTILS_HPP_

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/HashEntityAndEntityKey.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_DistributedVector.hpp>
#include <stk_math/StkVector.hpp>
#include <functional>
#include <unordered_map>
#include <vector>

namespace krino {

struct StkMeshEntities;

using DistributedEntityVector = DistributedQuantities<stk::mesh::Entity>;
using GradientFilter = std::function<void(const stk::mesh::Entity, const stk::math::Vector3d &, stk::math::Vector3d &)>;
using SolutionProjector = std::function<void(const stk::mesh::Entity, stk::math::Vector3d &)>;
using SolutionRejector = std::function<bool()>;
using NodeLocations = std::unordered_map<stk::mesh::Entity, stk::math::Vector3d>;

struct ElemNodeIndices
{
    ElemNodeIndices(const std::vector<unsigned> & elemNodeIndices, const std::vector<unsigned> & elemNodeStartingIndices, const unsigned elemIndex)
    : mBegin{&elemNodeIndices[elemNodeStartingIndices[elemIndex]]}, mEnd{&elemNodeIndices[elemNodeStartingIndices[elemIndex+1]]} {}
    typedef unsigned value_type;
    const value_type *mBegin;
    const value_type *mEnd;
    const value_type * begin() const { return mBegin; }
    const value_type * end() const { return mEnd; }
    size_t size() const { return mEnd - mBegin; }
    bool empty() const { return mEnd == mBegin; }
    value_type operator[](int i) const { return *(mBegin + i); }
};

class MeshElemNodeIndices
{
public:
    void build_for_mesh(const stk::mesh::BulkData & mesh,
        const std::vector<stk::mesh::Entity> & elems,
        const DistributedEntityVector & nodes);
    ElemNodeIndices element_node_indices(const unsigned elemIndex) const { return ElemNodeIndices(myElemNodeIndices, myElemNodeStartingIndices, elemIndex); }
    size_t num_elements() const { const size_t numEntry = myElemNodeStartingIndices.size(); return (numEntry == 0) ? 0 : (numEntry-1); }
private:
    std::vector<unsigned> myElemNodeIndices;
    std::vector<unsigned> myElemNodeStartingIndices;
};

DistributedEntityVector get_sorted_owned_nodes_and_unowned_shared_nodes(
    const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elemSelector);

// As get_sorted_owned_nodes_and_unowned_shared_nodes, with only the vertex nodes of the elements:
// the mid-edge nodes of quadratic elements (tet10, tri6) are left out, so that the elements are
// smoothed as their linear (tet4, tri3) vertices.
DistributedEntityVector get_sorted_owned_and_unowned_shared_vertex_nodes(
    const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elemSelector);

// Whether the node is a vertex of the elements it is in, rather than a mid-edge node.
bool is_vertex_node_of_its_elements(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node);

// The vertex nodes of the element: all its nodes for a linear element, the first 4 of a tet10 or
// the first 3 of a tri6.
StkMeshEntities
element_vertex_nodes(const stk::mesh::BulkData & mesh, const stk::mesh::Entity elem);

// The locations of the vertex nodes of the quadratic elements among elems, owned or not, to give
// set_mid_edge_nodes_of_moved_edges_to_edge_midpoints once the vertices are smoothed.  Empty for
// linear elements.
NodeLocations get_vertex_node_locations_of_quadratic_elements(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const std::vector<stk::mesh::Entity> & elems);

// Collective.  Moves each mid-edge node of the elements to the midpoint of its edge's vertex nodes
// when either vertex node has moved from its location in vertexNodeLocationsBeforeSmoothing, and
// syncs the coordinates.  The mid-edge node of an edge whose vertex nodes did not move is left
// where it is, so that an unchanged curved edge stays curved.  Nothing to do for linear elements.
// Each mid-edge node is moved by its owner, which need not own any element with that edge, so elems
// should be all the elements on this proc, owned or not.
void set_mid_edge_nodes_of_moved_edges_to_edge_midpoints(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const std::vector<stk::mesh::Entity> & elems,
    const NodeLocations & vertexNodeLocationsBeforeSmoothing);

// The index of node in nodes from get_sorted_owned_nodes_and_unowned_shared_nodes: in the owned
// (local) range for an owned node, otherwise local_size() plus its index in the shared (remote)
// range.  The node is only checked to be in nodes in debug builds; for a node not in nodes the
// index is not meaningful.
unsigned get_index_of_node_in_sorted_owned_or_shared_nodes(const stk::mesh::BulkData & mesh,
    const DistributedEntityVector & nodes,
    const stk::mesh::Entity node);

DistributedVector gather_node_coordinates(
    const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const DistributedEntityVector & nodes);

void communicate_gradient(
    const stk::mesh::BulkData & mesh,
    const DistributedEntityVector & nodes,
    DistributedVector & g);

void fill_elem_node_locations(const unsigned spatialDim,
    const DistributedVector& nodeLocs,
    const ElemNodeIndices & elemNodes,
    std::vector<const double*> & elemNodeLocs);
void fill_elem_node_locations(const unsigned spatialDim,
    const DistributedVector& nodeLocs,
    const MeshElemNodeIndices & meshElemNodeIndices,
    const unsigned elemIndex,
    std::vector<const double*> & elemNodeLocs);

void apply_gradient_filter(const unsigned dim,
    const GradientFilter & gradientFilter,
    const DistributedEntityVector & nodes,
    const DistributedVector & x,
    DistributedVector & g);

void apply_projection(const unsigned dim,
    const SolutionProjector & solutionProjector,
    const DistributedEntityVector & nodes,
    DistributedVector & x);

void set_node_coordinates(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const DistributedEntityVector & nodes,
    const DistributedVector & nodesCoords);

}

#endif /* KRINO_KRINO_SMOOTHING_AKRI_SMOOTHING_UTILS_HPP_ */
