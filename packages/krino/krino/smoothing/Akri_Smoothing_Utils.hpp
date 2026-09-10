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
#include <Akri_FieldRef.hpp>
#include <Akri_DistributedVector.hpp>
#include <stk_math/StkVector.hpp>
#include <functional>

namespace krino {

using DistributedEntityVector = DistributedQuantities<stk::mesh::Entity>;
using GradientFilter = std::function<void(const stk::mesh::Entity, const stk::math::Vector3d &, stk::math::Vector3d &)>;
using SolutionProjector = std::function<void(const stk::mesh::Entity, stk::math::Vector3d &)>;

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
