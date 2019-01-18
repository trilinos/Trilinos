// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
 // 
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
 // 
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef ZOLTANRCBADAPTER_HPP_
#define ZOLTANRCBADAPTER_HPP_

#include <Zoltan2_MeshAdapter.hpp>      // for MeshEntityType, etc
#include <stk_topology/topology.hpp>
#include <vector>
#include <stk_balance/internal/GeometricVertices.hpp>
#include <stk_util/util/ReportHandler.hpp>

typedef Zoltan2::BasicUserTypes<double, BalanceLocalNumber, BalanceGlobalNumber> learning_data_t;

class ZoltanGeometricAdapter : public Zoltan2::MeshAdapter<learning_data_t>
{
public:
    typedef Zoltan2::MeshAdapter<learning_data_t> base_adapter_t;


    ZoltanGeometricAdapter(const stk::balance::internal::GeometricVertices& geometricVertices) : mGeometricVertices(geometricVertices)
    {
    }

    virtual ~ZoltanGeometricAdapter() { }

    virtual size_t getLocalNumOf(Zoltan2::MeshEntityType etype) const
    {
        ThrowRequireMsg(get_mapped_topology(etype) == m_primary_rank, "getLocalNumOf couldn't return valid answer.");
        return mGeometricVertices.num_vertices();
    }

    virtual void getIDsViewOf(Zoltan2::MeshEntityType etype, BalanceGlobalNumber const *&Ids) const
    {
        ThrowRequireMsg(get_mapped_topology(etype) == m_primary_rank, "getIDsViewOf couldn't return valid answer.");
        Ids = mGeometricVertices.get_vertex_ids().data();
    }

    virtual int getDimension() const
    {
        return mGeometricVertices.get_spatial_dim();
    }

    virtual void getCoordinatesViewOf(Zoltan2::MeshEntityType etype, const scalar_t *&coords, int &stride, int coordDim) const
    {
        ThrowRequireMsg(get_mapped_topology(etype)== m_primary_rank, "Error!");
        coords = NULL;
        stride = getDimension();

        const std::vector<double> &vert_coords = mGeometricVertices.get_vertex_coords();
        if (!vert_coords.empty())
            coords = &vert_coords[coordDim];
    }

    virtual int getNumWeightsPerOf(Zoltan2::MeshEntityType etype) const
    {
        ThrowRequireMsg(get_mapped_topology(etype)== m_primary_rank, "Error!");
        return static_cast<int>(mGeometricVertices.getNumWeightsPerVertex());
    }

    virtual void getWeightsViewOf(Zoltan2::MeshEntityType etype, const scalar_t *&weights, int &stride, int idx = 0) const
    {
        ThrowRequireMsg(get_mapped_topology(etype)== m_primary_rank, "Error!");
        weights = NULL;
        stride = mGeometricVertices.getNumWeightsPerVertex();

        const std::vector<double> &vert_weights = mGeometricVertices.get_vertex_weights();
        if (!vert_weights.empty())
            weights = &vert_weights[idx];
    }

private:

    stk::mesh::EntityRank get_mapped_topology(Zoltan2::MeshEntityType etype) const
    {
        stk::mesh::EntityRank entityRank = stk::topology::INVALID_RANK;
        if(etype==Zoltan2::MESH_VERTEX)
        {
            entityRank = stk::topology::NODE_RANK;
        }
        else if(etype==Zoltan2::MESH_EDGE)
        {
            entityRank = stk::topology::EDGE_RANK;
        }
        else if(etype==Zoltan2::MESH_FACE)
        {
            entityRank = stk::topology::FACE_RANK;
        }
        else if(etype==Zoltan2::MESH_REGION)
        {
            entityRank = stk::topology::ELEM_RANK;
        }
        return entityRank;
    }

    stk::mesh::EntityRank m_primary_rank = stk::topology::ELEM_RANK;
    const stk::balance::internal::GeometricVertices &mGeometricVertices;
};



#endif /* ZOLTANRCBADAPTER_HPP_ */
