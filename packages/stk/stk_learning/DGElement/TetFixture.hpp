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

#ifndef STK_LEARNING_DGELEMENT_TETFIXTURE_HPP
#define STK_LEARNING_DGELEMENT_TETFIXTURE_HPP

#include <gtest/gtest.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

#include <stk_io/IossBridge.hpp>

class DGTetFixture : public ::testing::Test
{
protected:
    DGTetFixture()
    : communicator(MPI_COMM_WORLD), metaData(3), bulkData(nullptr)
    {
        tetPart = &metaData.declare_part_with_topology("TET", stk::topology::TETRAHEDRON_4);
        stk::io::put_io_part_attribute(*tetPart);
        skinPart = &metaData.declare_part("skin", stk::topology::FACE_RANK);
        nodePart = &metaData.declare_part_with_topology("node_part", stk::topology::NODE);
        coords = &metaData.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian>>(stk::topology::NODE_RANK, "coordinates");
        stk::mesh::put_field_on_mesh(*coords, metaData.universal_part(), metaData.spatial_dimension(), nullptr);
    }

    virtual ~DGTetFixture()
    {
        delete bulkData;
    }

    void setup_mesh(const std::vector<stk::mesh::EntityIdVector>& tet_conn, const std::vector<std::vector<double>> &node_coords, stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::NO_AUTO_AURA)
    {
        allocate_bulk(auraOption);
        stk::mesh::EntityVector nodes = setup_nodes(node_coords.size());
        setup_elements(tet_conn);
        initialize_coordinates(nodes, node_coords);
    }

    MPI_Comm get_comm()
    {
        return communicator;
    }

    stk::mesh::MetaData& get_meta()
    {
        return metaData;
    }

    stk::mesh::BulkData& get_bulk()
    {
        ThrowRequireMsg(bulkData!=nullptr, "Unit test error. Trying to get bulk data before it has been initialized.");
        return *bulkData;
    }

    void allocate_bulk(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        bulkData = new stk::mesh::BulkData(metaData, communicator, auraOption);
    }

    stk::mesh::Part* get_skin_part()
    {
        return skinPart;
    }

    stk::mesh::Field<double, stk::mesh::Cartesian>* get_coord_field()
    {
        return coords;
    }

private:

    stk::mesh::EntityVector setup_nodes(const size_t num_nodes)
    {
        stk::mesh::EntityVector nodes(num_nodes);
        get_bulk().modification_begin();
        for(unsigned int i=0;i<num_nodes;++i)
            nodes[i] = get_bulk().declare_node(i+1, stk::mesh::ConstPartVector{nodePart});
        get_bulk().modification_end();
        return nodes;
    }

    void setup_elements(const std::vector<stk::mesh::EntityIdVector>& tet_conn)
    {
        size_t num_tets = tet_conn.size();
        get_bulk().modification_begin();
        for(unsigned int i=0;i<num_tets;i++)
            stk::mesh::declare_element(get_bulk(), *tetPart, i+1, tet_conn[i]);
        get_bulk().modification_end();
    }

    void initialize_coordinates(const stk::mesh::EntityVector& nodes, const std::vector<std::vector<double>> &node_coords)
    {
        for(const stk::mesh::Entity node : nodes )
        {
            double *nodeCoord = stk::mesh::field_data(*coords, node);
            stk::mesh::EntityId id = get_bulk().identifier(node);
            nodeCoord[0] = node_coords[id-1][0];
            nodeCoord[1] = node_coords[id-1][1];
            nodeCoord[2] = node_coords[id-1][2];
        }
    }

    MPI_Comm communicator;
    stk::mesh::MetaData metaData;
    stk::mesh::BulkData *bulkData;
    stk::mesh::Part* tetPart = nullptr;
    stk::mesh::Part* nodePart = nullptr;
    stk::mesh::Part* skinPart = nullptr;
    stk::mesh::Field<double, stk::mesh::Cartesian>* coords = nullptr;
};


#endif
