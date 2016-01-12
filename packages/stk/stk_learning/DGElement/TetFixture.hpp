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
        stk::mesh::put_field(*coords, metaData.universal_part(), metaData.spatial_dimension());
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
            nodes[i] = get_bulk().declare_entity(stk::topology::NODE_RANK, i+1, *nodePart);
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
