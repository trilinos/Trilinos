/*
 *
 *
 * Stuff I didn't want to delete, but is utterly useless
 * Mostly just legacy stuff from before GoL class
 *
 */

#include <gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <stdlib.h>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include "stk_unit_test_utils/unittestMeshUtils.hpp"
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include "GameofLife.hpp"
#include "LodePNG.hpp"

#include <time.h>

typedef stk::mesh::Field<int> ScalarIntField;

std::string getGeneratedMeshString(const int xdim, const int ydim, const int zdim)
{
    std::ostringstream oss;
    oss << "generated: " << xdim << "x" << ydim << "x" << zdim;
    return oss.str();
}

// declarations
int get_num_active_neighbors(stk::mesh::BulkData& bulkData, stk::mesh::Entity& element, stk::mesh::Field<int>& lifeField);
void update_tri_life_val(stk::mesh::Entity& element, stk::mesh::Field<int>& lifeField, stk::mesh::Field<int>& activeNeighborField);
void update_2D_hex_life_val(stk::mesh::Entity& element, stk::mesh::Field<int>& lifeField, stk::mesh::Field<int>& activeNeighborField);
void activate_elements_with_these_ids(stk::mesh::BulkData& bulkData, stk::mesh::EntityIdVector& elemIdsToActivate, stk::mesh::Field<int>& field);
size_t write_output_mesh_and_return_file_handler(stk::io::StkMeshIoBroker& stkIo, std::string meshName, stk::mesh::Field<int>& lifeField);
void run_single_timestep(stk::mesh::BulkData& bulkData, stk::io::StkMeshIoBroker& stkIo, std::string& meshName, stk::mesh::EntityVector& elements,
                         stk::mesh::Field<int>& lifeField, stk::mesh::Field<int>& activeNeighborField, size_t fh, double timeStep);
void run_game_of_life(stk::mesh::BulkData& bulkData, stk::io::StkMeshIoBroker& stkIo, std::string& meshName, stk::mesh::EntityVector& elements,
                      stk::mesh::EntityIdVector& elemIdsToActivate, stk::mesh::Field<int>& lifeField, stk::mesh::Field<int>& activeNeighborField,
                      double totalTime, double timeStep);
void write_out_timestep(stk::io::StkMeshIoBroker& stkIo, size_t fh, double time);
void update_neighbor_field_data_of_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity elem, stk::mesh::Field<int>& activeNeighborField, stk::mesh::Field<int>& lifeField);
void update_life_field_data_of_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity elem, stk::mesh::Field<int>& activeNeighborField, stk::mesh::Field<int>& lifeField, stk::topology elemTopology);
void convert_normal_char_vector_into_2D_vector(const unsigned width, const unsigned height,
                                      const std::vector<unsigned char>& charVector,
                                      std::vector<std::vector<int>>& intVector);
// implementations
int get_num_active_neighbors(stk::mesh::BulkData& bulkData, stk::mesh::Entity& element, stk::mesh::Field<int>& lifeField)
{
    std::set<stk::mesh::Entity> activeNeighbors;
    int numNodes = bulkData.num_nodes(element);
    const stk::mesh::Entity* elemNode = bulkData.begin_nodes(element);
    for (int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
    {
        stk::mesh::Entity node = elemNode[nodeIndex];
        int numElems = bulkData.num_elements(node);
        const stk::mesh::Entity* nodeElem = bulkData.begin_elements(node);
        for (int index = 0; index < numElems; index++)
        {
            if (nodeElem[index] != element)
                if (*stk::mesh::field_data(lifeField, nodeElem[index]))
                    activeNeighbors.insert(nodeElem[index]);
        }
    }
    return activeNeighbors.size();
}
void update_tri_life_val(stk::mesh::Entity& element, stk::mesh::Field<int>& lifeField, stk::mesh::Field<int>& activeNeighborField)
{
    int* activeNeighborFieldVal = stk::mesh::field_data(activeNeighborField, element);
    int* lifeFieldVal = stk::mesh::field_data(lifeField, element);
    switch (*activeNeighborFieldVal)
    {
        case 2:
        case 7:
            break;
        case 3:
            *lifeFieldVal = 1;
            break;
        default:
            *lifeFieldVal = 0;
    }
}
void update_2D_hex_life_val(stk::mesh::Entity& element, stk::mesh::Field<int>& lifeField, stk::mesh::Field<int>& activeNeighborField)
{
    int* activeNeighborFieldVal = stk::mesh::field_data(activeNeighborField, element);
    int* lifeFieldVal = stk::mesh::field_data(lifeField, element);
    switch (*activeNeighborFieldVal)
    {
        case 4:
            break;
        case 5:
            *lifeFieldVal = 1;
            break;
        default:
            *lifeFieldVal = 0;
    }
}
void run_game_of_life(stk::mesh::BulkData& bulkData, stk::io::StkMeshIoBroker& stkIo, std::string& meshName, stk::mesh::EntityVector& elements,
                      stk::mesh::EntityIdVector& elemIdsToActivate, stk::mesh::Field<int>& lifeField, stk::mesh::Field<int>& activeNeighborField,
                      double totalTime, double timeStep)
{
    activate_elements_with_these_ids(bulkData, elemIdsToActivate, lifeField);
    size_t fh = write_output_mesh_and_return_file_handler(stkIo, meshName, lifeField);
    for (double time = 0; time < totalTime; time += timeStep)
        run_single_timestep(bulkData, stkIo, meshName, elements,lifeField, activeNeighborField, fh, time);
}
void activate_elements_with_these_ids(stk::mesh::BulkData& bulkData, stk::mesh::EntityIdVector& elemIdsToActivate, stk::mesh::Field<int>& field)
{
    size_t numElemsToActivate = elemIdsToActivate.size();
    for (size_t index = 0; index < numElemsToActivate; index++)
    {
        stk::mesh::Entity elem = bulkData.get_entity(stk::topology::ELEM_RANK, elemIdsToActivate[index]);
        if (bulkData.is_valid(elem))
            *stk::mesh::field_data(field, elem) = 1;
    }
}
size_t write_output_mesh_and_return_file_handler(stk::io::StkMeshIoBroker& stkIo, std::string meshName, stk::mesh::Field<int>& field)
{
    size_t fh = stkIo.create_output_mesh(meshName, stk::io::WRITE_RESULTS);
    stkIo.add_field(fh, field);
    stkIo.write_output_mesh(fh);
    return fh;
}
void run_single_timestep(stk::mesh::BulkData& bulkData, stk::io::StkMeshIoBroker& stkIo, std::string& meshName, stk::mesh::EntityVector& elements,
                         stk::mesh::Field<int>& lifeField, stk::mesh::Field<int>& activeNeighborField, size_t fh, double time)
{
    write_out_timestep(stkIo, fh, time);
    stk::mesh::communicate_field_data(bulkData, {&lifeField}); // screw you too eclipse

    for (size_t index = 0, elemSize = elements.size(); index < elemSize; index++)
        update_neighbor_field_data_of_element(bulkData, elements[index], activeNeighborField, lifeField);

    stk::topology elemTopology = bulkData.bucket(elements[0]).topology();
    for (size_t index = 0, elemSize = elements.size(); index < elemSize; index++)
        update_life_field_data_of_element(bulkData, elements[index], activeNeighborField, lifeField, elemTopology);
}
void write_out_timestep(stk::io::StkMeshIoBroker& stkIo, size_t fh, double timeStep)
{
    stkIo.begin_output_step(fh, timeStep);
    stkIo.write_defined_output_fields(fh);
    stkIo.end_output_step(fh);
}
void update_neighbor_field_data_of_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity elem, stk::mesh::Field<int>& activeNeighborField, stk::mesh::Field<int>& lifeField)
{
    if (bulkData.is_valid(elem))
    {
        int* elem_value = stk::mesh::field_data(activeNeighborField, elem);
        *elem_value = get_num_active_neighbors(bulkData, elem, lifeField);
    }
}
void update_life_field_data_of_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity elem, stk::mesh::Field<int>& activeNeighborField, stk::mesh::Field<int>& lifeField, stk::topology elemTopology)
{
    if (bulkData.is_valid(elem))
    {
        if (stk::topology::TRIANGLE_3 == elemTopology)
            update_tri_life_val(elem, lifeField, activeNeighborField);
        else
            update_2D_hex_life_val(elem, lifeField, activeNeighborField);
    }
}
void compute_triangle_node_id_array(stk::mesh::EntityIdVector nodeIds[], int pRank, int nodesPerRow, int rowsPerProc)
{
    size_t index;
    size_t elemsInFirstRow = 2*(nodesPerRow-1);
    size_t elemsAboveFirst = 2*(nodesPerRow-1)*rowsPerProc*pRank;

    for (index = 0; index < elemsInFirstRow; index+=2)
    {
        size_t trueIndex = index + elemsAboveFirst;
        int offset = nodesPerRow*pRank*rowsPerProc-index/2;
        nodeIds[trueIndex].resize(3);
        nodeIds[trueIndex+1].resize(3);
        nodeIds[trueIndex] = {index + 1 + offset, index + 2 + offset, nodesPerRow + 1 + index + offset};
        nodeIds[trueIndex+1] = {nodesPerRow + index + 2 + offset, nodesPerRow + index + 1 + offset, index + 2 + offset};
    }
    size_t numElemsInProc = elemsInFirstRow*rowsPerProc;
    for (; index < numElemsInProc; index++)
    {
        size_t trueIndex = index + elemsAboveFirst;
        nodeIds[trueIndex].resize(3);
        for (size_t nodeIndex = 0; nodeIndex < 3; nodeIndex++)
        {
            nodeIds[trueIndex][nodeIndex] = nodeIds[index%elemsInFirstRow+elemsAboveFirst][nodeIndex] + index/elemsInFirstRow*nodesPerRow;
        }
    }
}
TEST(GameOfLife, DISABLED_SetupMesh)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int pSize = stk::parallel_machine_size(comm);

    if (pSize == 1)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulkData(meta, comm);

        std::string meshSpec = "generated:4x4x1";
        stk::io::StkMeshIoBroker exodusFileReader(comm);

        exodusFileReader.set_bulk_data(bulkData);
        exodusFileReader.add_mesh_database(meshSpec, stk::io::READ_MESH);
        exodusFileReader.create_input_mesh();

        ScalarIntField& lifeField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Status");
        int initialVal = 0;
        stk::mesh::put_field(lifeField, meta.universal_part(), &initialVal);

        exodusFileReader.populate_bulk_data();
        stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEMENT_RANK, 1);

        EXPECT_TRUE(bulkData.is_valid(elem1));

    }
}
TEST(GameofLife, DISABLED_singleTimeStepUnitTestWithExodus)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int pSize = stk::parallel_machine_size(comm);

    if (pSize == 1)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulkData(meta, comm);

        std::string meshSpec = "generated:2x2x1";
        stk::io::StkMeshIoBroker stkIo(comm);

        stkIo.set_bulk_data(bulkData);
        stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
        stkIo.create_input_mesh();

        ScalarIntField& lifeField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Status");
        int initialVal = 0;
        stk::mesh::put_field(lifeField, meta.universal_part(), &initialVal);

        stkIo.populate_bulk_data();

        std::vector<stk::mesh::Entity> elements;
        stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::ELEM_RANK, elements);

        for (size_t index = 0; index < elements.size(); index++)
        {
            int* elem_value = stk::mesh::field_data(lifeField, elements[index]);
            *elem_value = stkIo.bulk_data().identifier(elements[index]);
        }

        EXPECT_EQ(1, *stk::mesh::field_data(lifeField, elements[0]));
        EXPECT_EQ(2, *stk::mesh::field_data(lifeField, elements[1]));
        EXPECT_EQ(3, *stk::mesh::field_data(lifeField, elements[2]));
        EXPECT_EQ(4, *stk::mesh::field_data(lifeField, elements[3]));

        std::string mesh_name = "grid1.e";

        double time = 0.1;

        size_t fh = stkIo.create_output_mesh(mesh_name, stk::io::WRITE_RESULTS);
        stkIo.add_field(fh, lifeField);
        stkIo.write_output_mesh(fh);

        stkIo.begin_output_step(fh, time);
        stkIo.write_defined_output_fields(fh);
        stkIo.end_output_step(fh);
    }
}
TEST(GameofLife, DISABLED_multiTimeStepUnitTestWithExodus)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int pSize = stk::parallel_machine_size(comm);

    if (pSize == 1)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulkData(meta, comm);

        std::string meshSpec = "generated:2x2x1";
        stk::io::StkMeshIoBroker stkIo(comm);

        stkIo.set_bulk_data(bulkData);
        stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
        stkIo.create_input_mesh();

        ScalarIntField& lifeField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Status");
        int initialVal = 0;
        stk::mesh::put_field(lifeField, meta.universal_part(), &initialVal);

        stkIo.populate_bulk_data();

        std::vector<stk::mesh::Entity> elements;
        stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::ELEM_RANK, elements);

        std::string mesh_name = "grid2.e";

        size_t fh = stkIo.create_output_mesh(mesh_name, stk::io::WRITE_RESULTS);
        stkIo.add_field(fh, lifeField);
        stkIo.write_output_mesh(fh);

        for (double timeStep = 1; timeStep <= 5; timeStep++)
        {
            for (size_t index = 0; index < elements.size(); index++)
            {
                int* elem_value = stk::mesh::field_data(lifeField, elements[index]);
                *elem_value = stkIo.bulk_data().identifier(elements[index])*timeStep;
            }
            stkIo.begin_output_step(fh, timeStep);
            stkIo.write_defined_output_fields(fh);
            stkIo.end_output_step(fh);
        }
    }
}
TEST(GameofLife, DISABLED_movingElementsFromPartToPart)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int pSize = stk::parallel_machine_size(comm);

    if (pSize == 1)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulkData(meta, comm);

        stk::mesh::Part& activePart = meta.declare_part_with_topology("Active", stk::topology::HEX_8);
        std::string meshSpec = "generated:2x2x1";
        stk::io::StkMeshIoBroker stkIo(comm);

        stkIo.set_bulk_data(bulkData);
        stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
        stkIo.create_input_mesh();

        ScalarIntField& lifeField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Status");
        int initialVal = 0;
        stk::mesh::put_field(lifeField, meta.universal_part(), &initialVal);

        stkIo.populate_bulk_data();
        stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEMENT_RANK, 1);

        EXPECT_TRUE(bulkData.is_valid(elem1));

        stk::mesh::EntityVector elements;
        stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::ELEM_RANK, elements);
        stk::unit_test_util::put_mesh_into_part(bulkData, activePart);

        EXPECT_TRUE(bulkData.bucket(elements[0]).member(activePart));
        EXPECT_TRUE(bulkData.bucket(elements[1]).member(activePart));
        EXPECT_TRUE(bulkData.bucket(elements[2]).member(activePart));
        EXPECT_TRUE(bulkData.bucket(elements[3]).member(activePart));

        for (size_t index = 0; index < elements.size(); index++)
        {
            int* elem_value = stk::mesh::field_data(lifeField, elements[index]);
            *elem_value = stkIo.bulk_data().identifier(elements[index]);
        }

        stk::mesh::PartVector activeVector = {&activePart};
        stk::mesh::PartVector emptyVector;

        bulkData.modification_begin();
        for (size_t index = 0; index < elements.size(); index++)
        {
            int* elemId = stk::mesh::field_data(lifeField, elements[index]);
            if (*elemId > 2)
            {
                bulkData.change_entity_parts(elements[index], emptyVector, activeVector);
            }
        }

        bulkData.modification_end();

        EXPECT_TRUE(bulkData.bucket(elements[0]).member(activePart));
        EXPECT_TRUE(bulkData.bucket(elements[1]).member(activePart));
        EXPECT_FALSE(bulkData.bucket(elements[2]).member(activePart));
        EXPECT_FALSE(bulkData.bucket(elements[3]).member(activePart));
    }
}
TEST(GameofLife, DISABLED_loopingOverNodes)
{

    MPI_Comm comm = MPI_COMM_WORLD;
    int pSize = stk::parallel_machine_size(comm);

    if (pSize == 1)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulkData(meta, comm);

        stk::mesh::Part& activePart = meta.declare_part_with_topology("Active", stk::topology::HEX_8);
        std::string meshSpec = "generated:2x2x1";
        stk::io::StkMeshIoBroker stkIo(comm);

        stkIo.set_bulk_data(bulkData);
        stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
        stkIo.create_input_mesh();

        ScalarIntField& lifeField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Status");
        int initVal = 0;
        stk::mesh::put_field(lifeField, meta.universal_part(), &initVal);


        stkIo.populate_bulk_data();

        stk::mesh::EntityVector elements;
        stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::ELEM_RANK, elements);

        stk::mesh::EntityVector nodes;
        stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK, nodes);

        stk::unit_test_util::put_mesh_into_part(bulkData, activePart);

        stk::mesh::PartVector activeVector = {&activePart};
        stk::mesh::PartVector emptyVector;

        bulkData.modification_begin();
        bulkData.change_entity_parts(elements[1], emptyVector, activeVector);
        bulkData.modification_end();

        // will be function
        std::set<stk::mesh::Entity> activeConnectedEntities;
        for (size_t elemIndex = 0; elemIndex < elements.size(); elemIndex++)
        {
            activeConnectedEntities.clear();
            int numNodes = bulkData.num_nodes(elements[elemIndex]);
            const stk::mesh::Entity* elemNode = bulkData.begin_nodes(elements[elemIndex]);
            for (int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
            {
                stk::mesh::Entity node = elemNode[nodeIndex];
                int numElems = bulkData.num_elements(node);
                const stk::mesh::Entity* nodeElem = bulkData.begin_elements(node);
                for (int index = 0; index < numElems; index++)
                {
                    if (nodeElem[index] != elements[elemIndex])
                    {
                        if (bulkData.bucket(nodeElem[index]).member(activePart))
                        {
                            activeConnectedEntities.insert(nodeElem[index]);
                        }
                    }
                }
            }
            int* elem_value = stk::mesh::field_data(lifeField, elements[elemIndex]);
            *elem_value = activeConnectedEntities.size();
        }

        EXPECT_EQ(2, *stk::mesh::field_data(lifeField, elements[0]));
        EXPECT_EQ(3, *stk::mesh::field_data(lifeField, elements[1]));
        EXPECT_EQ(2, *stk::mesh::field_data(lifeField, elements[2]));
        EXPECT_EQ(2, *stk::mesh::field_data(lifeField, elements[3]));
    }
}
TEST(GameofLife, DISABLED_deactivatingParts)
{

    MPI_Comm comm = MPI_COMM_WORLD;
    int pSize = stk::parallel_machine_size(comm);

    if (pSize == 1)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulkData(meta, comm);

        stk::mesh::Part& activePart = meta.declare_part_with_topology("Active", stk::topology::HEX_8);
        std::string meshSpec = "generated:3x3x1";
        stk::io::StkMeshIoBroker stkIo(comm);

        stkIo.set_bulk_data(bulkData);
        stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
        stkIo.create_input_mesh();

        ScalarIntField& activeNeighborField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Neighbor Field");


        int initVal = 0;
        stk::mesh::put_field(activeNeighborField, meta.universal_part(), &initVal);


        stkIo.populate_bulk_data();

        stk::mesh::EntityVector elements;
        stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::ELEM_RANK, elements);

        stk::mesh::EntityVector nodes;
        stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK, nodes);

        stk::unit_test_util::put_mesh_into_part(bulkData, activePart);


        stk::mesh::PartVector activeVector = {&activePart};
        stk::mesh::PartVector emptyVector;
        bulkData.modification_begin();
        bulkData.change_entity_parts(elements[0], emptyVector, activeVector);
        bulkData.change_entity_parts(elements[1], emptyVector, activeVector);
        bulkData.change_entity_parts(elements[2], emptyVector, activeVector);
        bulkData.change_entity_parts(elements[6], emptyVector, activeVector);
        bulkData.change_entity_parts(elements[7], emptyVector, activeVector);
        bulkData.change_entity_parts(elements[8], emptyVector, activeVector);
        bulkData.modification_end();

        // will be function
        std::set<stk::mesh::Entity> activeConnectedEntities;
        for (size_t elemIndex = 0; elemIndex < elements.size(); elemIndex++)
        {
            activeConnectedEntities.clear();
            int numNodes = bulkData.num_nodes(elements[elemIndex]);
            const stk::mesh::Entity* elemNode = bulkData.begin_nodes(elements[elemIndex]);
            for (int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
            {
                stk::mesh::Entity node = elemNode[nodeIndex];
                int numElems = bulkData.num_elements(node);
                const stk::mesh::Entity* nodeElem = bulkData.begin_elements(node);
                for (int index = 0; index < numElems; index++)
                {
                    if (nodeElem[index] != elements[elemIndex])
                    {
                        if (bulkData.bucket(nodeElem[index]).member(activePart))
                        {
                            activeConnectedEntities.insert(nodeElem[index]);
                        }
                    }
                }
            }
            int* elem_value = stk::mesh::field_data(activeNeighborField, elements[elemIndex]);
            *elem_value = activeConnectedEntities.size();
        }

        bulkData.modification_begin();
        for (size_t elemIndex = 0; elemIndex < elements.size(); elemIndex++)
        {
            int fieldVal = *stk::mesh::field_data(activeNeighborField, elements[elemIndex]);
            switch (fieldVal)
            {
                case 2:
                    break;
                case 3:
                    bulkData.change_entity_parts(elements[elemIndex], activeVector, emptyVector);
                    break;
                default:
                    bulkData.change_entity_parts(elements[elemIndex], emptyVector, activeVector);
            }
        }
        bulkData.modification_end();

        EXPECT_TRUE(bulkData.bucket(elements[1]).member(activePart));
        EXPECT_TRUE(bulkData.bucket(elements[4]).member(activePart));
        EXPECT_TRUE(bulkData.bucket(elements[7]).member(activePart));
        EXPECT_FALSE(bulkData.bucket(elements[0]).member(activePart));
        EXPECT_FALSE(bulkData.bucket(elements[2]).member(activePart));
        EXPECT_FALSE(bulkData.bucket(elements[3]).member(activePart));
        EXPECT_FALSE(bulkData.bucket(elements[5]).member(activePart));
        EXPECT_FALSE(bulkData.bucket(elements[6]).member(activePart));
        EXPECT_FALSE(bulkData.bucket(elements[8]).member(activePart));
    }
}
TEST(GameofLife, DISABLED_deactivatingEntitiesUsingFieldAndFunctions)
{

    MPI_Comm comm = MPI_COMM_WORLD;
    int pSize = stk::parallel_machine_size(comm);

    if (pSize == 1)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulkData(meta, comm);

        std::string meshSpec = "generated:3x3x1";
        stk::io::StkMeshIoBroker stkIo(comm);

        stkIo.set_bulk_data(bulkData);
        stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
        stkIo.create_input_mesh();

        ScalarIntField& activeNeighborField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Active Neighbors");
        ScalarIntField& lifeField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Field");

        int initNeighborVal = 0;
        int initLifeVal = 1;

        stk::mesh::put_field(activeNeighborField, meta.universal_part(), &initNeighborVal);
        stk::mesh::put_field(lifeField, meta.universal_part(), &initLifeVal);
        stkIo.populate_bulk_data();

        stk::mesh::EntityVector elements;
        stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::ELEM_RANK, elements);

        stk::mesh::EntityVector nodes;
        stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK, nodes);

        // start with some dead elements
        *stk::mesh::field_data(lifeField, elements[0]) = 0;
        *stk::mesh::field_data(lifeField, elements[1]) = 0;
        *stk::mesh::field_data(lifeField, elements[2]) = 0;
        *stk::mesh::field_data(lifeField, elements[6]) = 0;
        *stk::mesh::field_data(lifeField, elements[7]) = 0;
        *stk::mesh::field_data(lifeField, elements[8]) = 0;

        int elemSize = elements.size();
        for (int index = 0; index < elemSize; index++)
        {
            int* elem_value = stk::mesh::field_data(activeNeighborField, elements[index]);
            *elem_value = get_num_active_neighbors(bulkData, elements[index], lifeField);
        }

        for (int index = 0; index < elemSize; index++)
        {
            update_2D_hex_life_val(elements[index], lifeField, activeNeighborField);
        }

        EXPECT_EQ(1, *stk::mesh::field_data(lifeField, elements[1]));
        EXPECT_EQ(1, *stk::mesh::field_data(lifeField, elements[4]));
        EXPECT_EQ(1, *stk::mesh::field_data(lifeField, elements[7]));
        EXPECT_EQ(0, *stk::mesh::field_data(lifeField, elements[0]));
        EXPECT_EQ(0, *stk::mesh::field_data(lifeField, elements[2]));
        EXPECT_EQ(0, *stk::mesh::field_data(lifeField, elements[3]));
        EXPECT_EQ(0, *stk::mesh::field_data(lifeField, elements[5]));
        EXPECT_EQ(0, *stk::mesh::field_data(lifeField, elements[6]));
        EXPECT_EQ(0, *stk::mesh::field_data(lifeField, elements[8]));
    }
}
TEST(GameofLife, DISABLED_multipleTimeStepsWithExodus)
{

    MPI_Comm comm = MPI_COMM_WORLD;
    int pSize = stk::parallel_machine_size(comm);

    if (pSize == 1)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulkData(meta, comm);

        std::string meshSpec = "generated:10x10x1";
        stk::io::StkMeshIoBroker stkIo(comm);

        ScalarIntField& activeNeighborField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Active Neighbors");
        ScalarIntField& lifeField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Field");
        int initNeighborVal = 0;
        int initLifeVal = 0;
        stk::mesh::put_field(activeNeighborField, meta.universal_part(), &initNeighborVal);
        stk::mesh::put_field(lifeField, meta.universal_part(), &initLifeVal);

        stkIo.set_bulk_data(bulkData);
        stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
        stkIo.create_input_mesh();
        stkIo.populate_bulk_data();

        stk::mesh::EntityVector elements;
        stk::mesh::get_entities(bulkData, stk::topology::ELEM_RANK, elements);

        // start with some live elements
        stk::mesh::EntityIdVector elemIdsToActivate1 = {45, 35, 55, 44, 56}; //r pentomino
        stk::mesh::EntityIdVector elemIdsToActivate2 = {45, 35, 55}; // straight oscillator
        stk::mesh::EntityIdVector elemIdsToActivate3 = {46, 47, 35, 34, 37, 25, 23, 27, 53, 65, 67, 64, 63}; //weird thing

        std::string meshName = "10x10GameofLife.e";
        double totalTime = 10;
        double timeStep = .1;
        run_game_of_life(bulkData, stkIo, meshName, elements, elemIdsToActivate3, lifeField, activeNeighborField, totalTime, timeStep);
    }
}
TEST(GameofLife, DISABLED_20x20GameofLife)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int pSize = stk::parallel_machine_size(comm);

    if (pSize == 1)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulkData(meta, comm);

        std::string meshSpec = "generated:20x20x1";
        stk::io::StkMeshIoBroker stkIo(comm);

        ScalarIntField& activeNeighborField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Active Neighbors");
        ScalarIntField& lifeField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Field");
        int initNeighborVal = 0;
        int initLifeVal = 0;
        stk::mesh::put_field(activeNeighborField, meta.universal_part(), &initNeighborVal);
        stk::mesh::put_field(lifeField, meta.universal_part(), &initLifeVal);

        stkIo.set_bulk_data(bulkData);
        stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
        stkIo.create_input_mesh();
        stkIo.populate_bulk_data();

        stk::mesh::EntityVector elements;
        stk::mesh::get_entities(bulkData, stk::topology::ELEM_RANK, elements);

        // start with some live elements
        stk::mesh::EntityIdVector elemIdsToActivate1 = {211, 212, 190, 189, 192, 170, 168, 172, 228, 248, 249, 250, 252}; //weird thing
        stk::mesh::EntityIdVector elemIdsToActivate2 = // dinner table
        {210, 211, 209, 190, 191, 189, //Middle section
         167, 166, 147, 145, 125, 105, 104, //bottom left
         152, 153, 132, 113, 114, 115, 95,  // bottom right
         253, 254, 273, 275, 295, 315, 316, // top right
         268, 267, 288, 307, 306, 305, 325 // top left
        };

        std::string meshName1 = "20x20GameofLifeInfinite.e";
        std::string meshName2 = "20x20GameofLifeDinnerTable.e";
        double totalTime = 10;
        double timeStep = .1;
        run_game_of_life(bulkData, stkIo, meshName1, elements, elemIdsToActivate1, lifeField, activeNeighborField, totalTime, timeStep);
    }
}
TEST(GameofLife, DISABLED_ScalableTriangleMesh)
{
    const unsigned spacialDim = 2;
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int pSize = stk::parallel_machine_size(pm);
    if (pSize > 1)
    {
        return;
    }

    stk::mesh::MetaData meta(spacialDim);
    stk::mesh::BulkData bulkData(meta, pm);
    stk::mesh::Part* triPart = &meta.declare_part_with_topology("Triangle_Part", stk::topology::TRIANGLE_3);
    stk::io::put_io_part_attribute(*triPart);

    stk::mesh::Field<double,stk::mesh::Cartesian2d>& nodeCoord =
            meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian2d>>(stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field(nodeCoord, meta.universal_part(), 2);

    ScalarIntField& activeNeighborField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Active Neighbors");
    ScalarIntField& lifeField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Field");
    int initNeighborVal = 0;
    int initLifeVal = 0;
    stk::mesh::put_field(activeNeighborField, meta.universal_part(), &initNeighborVal);
    stk::mesh::put_field(lifeField, meta.universal_part(), &initLifeVal);

    meta.commit();

    const size_t nodesPerRow = 25;
    const size_t nodesPerCol = 25;

    const size_t numElems = 2*(nodesPerRow-1)*(nodesPerCol-1);
    stk::mesh::EntityId elemIds[numElems];
    for (size_t index = 0; index < numElems; index++)
    {
        elemIds[index] = index+1;
    }

    stk::mesh::EntityIdVector nodeIds[numElems];
    size_t index;
    size_t elemsInFirstRow = 2*(nodesPerRow-1);

    for (index = 0; index < elemsInFirstRow; index+=2) // create the first row
            {
        int offset = -index/2;
        nodeIds[index].resize(3);
        nodeIds[index+1].resize(3);
        nodeIds[index] = {index + 1 + offset, index + 2 + offset, nodesPerRow + 1 + index + offset};
        nodeIds[index+1] = {nodesPerRow + index + 2 + offset, nodesPerRow + index + 1 + offset, index + 2 + offset};
            }

    for (; index < numElems; index++)
    {
        nodeIds[index].resize(3);
        for (size_t nodeIndex = 0; nodeIndex < 3; nodeIndex++)
        {
            nodeIds[index][nodeIndex] = nodeIds[index%elemsInFirstRow][nodeIndex] + index/elemsInFirstRow*nodesPerRow;
        }
    }

    bulkData.modification_begin();
    for (size_t index = 0; index < numElems; index++)
    {
        stk::mesh::declare_element(bulkData, *triPart, elemIds[index], nodeIds[index]);
    }

    stk::mesh::EntityVector nodes;
    stk::mesh::get_entities(bulkData, stk::topology::NODE_RANK, nodes);
    size_t numNodes = nodes.size();

    for (size_t index = 0; index < numNodes; index++)
    {
        double* const coord = stk::mesh::field_data(nodeCoord, nodes[index]);
        coord[0]= index%nodesPerRow;
        coord[1]= index/nodesPerRow;
    }
    bulkData.modification_end();

    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(bulkData, stk::topology::ELEM_RANK, elements);

    stk::io::StkMeshIoBroker stkIo(pm);
    stkIo.set_bulk_data(bulkData);

    std::string meshName ="ScalableTriangleMesh.e";
    stk::mesh::EntityIdVector elemIdsToActivate;

    for (size_t index = 1; index < elements.size(); index++)
    {
        if ((!(index%5)) || (!(index%13)) || (!index&11))
            elemIdsToActivate.push_back(index);
    }

    const double totalTime = 10;
    const double timeStep = .1;
    run_game_of_life(bulkData, stkIo, meshName, elements, elemIdsToActivate, lifeField, activeNeighborField, totalTime, timeStep);
}
TEST(GameofLife, DISABLED_ParallelTriangleMesh)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int procNum = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    if (8 == procNum)
    {
        stk::mesh::MetaData meta(2);
        stk::mesh::BulkData bulkData(meta, comm);
        stk::mesh::Part* triPart = &meta.declare_part_with_topology("Triangle_Part", stk::topology::TRIANGLE_3);
        stk::io::put_io_part_attribute(*triPart);

        stk::mesh::Field<double,stk::mesh::Cartesian2d>& nodeCoord =
                meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian2d>>(stk::topology::NODE_RANK, "coordinates");
        stk::mesh::put_field(nodeCoord, meta.universal_part(), 2);

        ScalarIntField& activeNeighborField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Active Neighbors");
        ScalarIntField& lifeField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Field");
        int initNeighborVal = 0;
        int initLifeVal = 0;
        stk::mesh::put_field(activeNeighborField, meta.universal_part(), &initNeighborVal);
        stk::mesh::put_field(lifeField, meta.universal_part(), &initLifeVal);

        meta.commit();

        const unsigned elemRowsPerProc = 12;
        const size_t nodesPerRow = 101, nodesPerCol = (elemRowsPerProc*procNum)+1;
        const size_t numElems = 2*(nodesPerRow-1)*(nodesPerCol-1);

        stk::mesh::EntityIdVector elemIds;
        elemIds.resize(numElems);
        for (size_t index = 0; index < numElems; index++)
        {
            elemIds[index] = index+1;
        }

        stk::mesh::EntityIdVector nodeIds[numElems];
        compute_triangle_node_id_array(nodeIds, procRank, nodesPerRow, elemRowsPerProc);

        bulkData.modification_begin();

        unsigned numElemsPerProc = elemRowsPerProc*2*(nodesPerRow-1);
        unsigned elemOffset = numElemsPerProc*procRank;

        for (size_t index = 0; index < numElemsPerProc; index++)
            declare_element(bulkData, *triPart, elemIds[index+elemOffset], nodeIds[index+elemOffset]);

        if (procRank != 0)
        {
            unsigned bottomNodeOffset = elemRowsPerProc*nodesPerRow*procRank;
            for (unsigned index = 1; index <= nodesPerRow; index++)
            {
                stk::mesh::Entity bottomNode = bulkData.get_entity(stk::topology::NODE_RANK, index + bottomNodeOffset);
                bulkData.add_node_sharing(bottomNode, procRank-1);
            }
        }
        if (procNum-1 != procRank)
        {
            unsigned topNodeOffset = elemRowsPerProc*nodesPerRow*(procRank+1);
            for (unsigned index = 1; index <= nodesPerRow; index++)
            {
                stk::mesh::Entity topNode = bulkData.get_entity(stk::topology::NODE_RANK, index + topNodeOffset);
                bulkData.add_node_sharing(topNode, procRank+1);
            }
        }

        stk::mesh::EntityVector nodes;
        stk::mesh::get_entities(bulkData, stk::topology::NODE_RANK, nodes);
        size_t numNodes = nodes.size();
        for (size_t index = 0; index < numNodes; index++)
        {
            if (bulkData.bucket(nodes[index]).owned())
            {
                double* const coord = stk::mesh::field_data(nodeCoord, nodes[index]);
                coord[0] = index%nodesPerRow;
                coord[1] = index/nodesPerRow + elemRowsPerProc*procRank;
            }
        }
        bulkData.modification_end();

        std::string meshName = "ParallelTrianguleMesh.e";
        stk::io::StkMeshIoBroker stkIo(comm);
        stkIo.set_bulk_data(bulkData);

        stk::mesh::EntityVector elements;
        stk::mesh::get_entities(bulkData, stk::topology::ELEM_RANK, elements);
        stk::mesh::EntityIdVector elemIdsToActivate;

        for (size_t index = 1; index <= numElems; index++)
        {
            if ((!(index%19)) || (!(index%13)) || (!index&11))
                elemIdsToActivate.push_back(index);
        }

        const double totalTime = 100;
        const double timeStep = .1;
        run_game_of_life(bulkData, stkIo, meshName, elements, elemIdsToActivate, lifeField, activeNeighborField, totalTime, timeStep);
    }
}
TEST(GameofLife, DISABLED_Parallel10x10W2procs)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int pSize = stk::parallel_machine_size(comm);
    if (pSize == 2)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulkData(meta, comm, stk::mesh::BulkData::AUTO_AURA);

        std::string meshSpec = "generated:1x10x10";
        stk::io::StkMeshIoBroker stkIo(comm);

        ScalarIntField& activeNeighborField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Active Neighbors");
        ScalarIntField& lifeField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Field");
        int initNeighborVal = 0;
        int initLifeVal = 0;
        stk::mesh::put_field(activeNeighborField, meta.universal_part(), &initNeighborVal);
        stk::mesh::put_field(lifeField, meta.universal_part(), &initLifeVal);

        stkIo.set_bulk_data(bulkData);
        stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
        stkIo.create_input_mesh();
        stkIo.populate_bulk_data();

        std::string meshName = "Parallel10x10W2procs.e";

        stk::mesh::EntityVector elements;
        stk::mesh::get_entities(bulkData, stk::topology::ELEM_RANK, elements);
        stk::mesh::EntityIdVector elemIdsToActivate = {46, 47, 35, 34, 37, 25, 23, 27, 53, 65, 67, 64, 63};

        double totalTime = 1;
        double timeStep = .01;
        run_game_of_life(bulkData, stkIo, meshName, elements,elemIdsToActivate, lifeField, activeNeighborField,totalTime, timeStep);
    }

}
TEST(GameofLife, DISABLED_QuadParallelGameofLife)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    const int procSize = stk::parallel_machine_size(comm);
    const int procRank = stk::parallel_machine_rank(comm);
    if (8 == procSize)
    {
        stk::mesh::MetaData meta(2);
        stk::mesh::BulkData bulkData(meta, comm, stk::mesh::BulkData::AUTO_AURA);

        stk::mesh::Part* quadPart = &meta.declare_part_with_topology("Quad_Part", stk::topology::QUAD_4);
        stk::io::put_io_part_attribute(*quadPart);

        stk::mesh::Field<double,stk::mesh::Cartesian2d>& nodeCoord =
                meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian2d>>(stk::topology::NODE_RANK, "coordinates");
        stk::mesh::put_field(nodeCoord, meta.universal_part(), 2);

        ScalarIntField& activeNeighborField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Active Neighbors");
        ScalarIntField& lifeField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Field");
        int initNeighborVal = 0;
        int initLifeVal = 0;
        stk::mesh::put_field(activeNeighborField, meta.universal_part(), &initNeighborVal);
        stk::mesh::put_field(lifeField, meta.universal_part(), &initLifeVal);

        meta.commit();

        const unsigned nodesPerRow = 101;
        const unsigned elemsPerRow = nodesPerRow-1;
        const unsigned elemRowsPerProc = 12;
        const unsigned numElemsPerProc = elemRowsPerProc*elemsPerRow;
        const unsigned numElems = numElemsPerProc*procSize;
        const unsigned procElemIndexOffset = numElemsPerProc*procRank;

        stk::mesh::EntityId elemIds[numElems];
        stk::mesh::EntityIdVector elemNodeIds [numElems];

        unsigned nodeIdOffset = nodesPerRow*procRank*elemRowsPerProc;
        unsigned index;
        for (index = 0; index < elemsPerRow; index++)
        {
            elemIds[index + procElemIndexOffset] = index + procElemIndexOffset+1;
            elemNodeIds[index + procElemIndexOffset].resize(4);
            elemNodeIds[index + procElemIndexOffset] = {nodeIdOffset+index+1, nodeIdOffset+index+2, nodeIdOffset+index+2+nodesPerRow, nodeIdOffset+index+1+nodesPerRow};
        }
        for (; index < numElemsPerProc; index++)
        {
            elemIds[index + procElemIndexOffset] = index + procElemIndexOffset+1;
            elemNodeIds[index + procElemIndexOffset].resize(4);
            for (unsigned nodeIndex = 0; nodeIndex < 4; nodeIndex++)
                elemNodeIds[index + procElemIndexOffset][nodeIndex] = elemNodeIds[index%elemsPerRow + procElemIndexOffset][nodeIndex] + index/elemsPerRow * nodesPerRow;
        }

        bulkData.modification_begin();
        for (unsigned index = 0; index < numElemsPerProc; index++)
            declare_element(bulkData, *quadPart, elemIds[index + procElemIndexOffset], elemNodeIds[index + procElemIndexOffset]);

        stk::mesh::EntityVector elements;
        stk::mesh::get_entities(bulkData, stk::topology::ELEM_RANK, elements);
        stk::mesh::EntityVector nodes;
        stk::mesh::get_entities(bulkData, stk::topology::NODE_RANK, nodes);

        stk::mesh::EntityIdVector elementsIds;
        for (auto R: elements)
            elementsIds.push_back(bulkData.identifier(R));

        stk::mesh::EntityIdVector nodesIds;
        for (auto R: nodes)
            nodesIds.push_back(bulkData.identifier(R));

        if (0 != procRank)
        {
            unsigned bottomNodeIdOffset = nodesPerRow*procRank*elemRowsPerProc;
            for (unsigned nodeId = 1; nodeId <= nodesPerRow; nodeId++)
            {
                stk::mesh::Entity node = bulkData.get_entity(stk::topology::NODE_RANK, nodeId + bottomNodeIdOffset);
                ThrowRequire(bulkData.is_valid(node));
                bulkData.add_node_sharing(node, procRank-1);
            }
        }
        if (procSize-1 != procRank)
        {
            unsigned topNodeIdOffset = nodesPerRow*elemRowsPerProc*(procRank+1);
            for (unsigned nodeId = 1; nodeId <= nodesPerRow; nodeId++)
            {
                stk::mesh::Entity node = bulkData.get_entity(stk::topology::NODE_RANK, nodeId + topNodeIdOffset);
                ThrowRequire(bulkData.is_valid(node));
                bulkData.add_node_sharing(node, procRank+1);
            }
        }

        size_t numNodes = nodes.size();
        for (size_t index = 0; index < numNodes; index++)
        {
            if (bulkData.is_valid(nodes[index]))
            {
                double* const coord = stk::mesh::field_data(nodeCoord, nodes[index]);
                coord[0] = index%nodesPerRow;
                coord[1] = index/nodesPerRow + elemRowsPerProc*procRank;
            }
        }
        bulkData.modification_end();

        stk::mesh::EntityIdVector elemIdsToActivate;
        for (size_t index = 1; index <= numElems; index++)
        {
            if ((!(index%19)) || (!(index%13)) || (!index&11))
                elemIdsToActivate.push_back(index);
        }

        std::string meshName = "QuadParallelGameofLife.e";
        stk::io::StkMeshIoBroker stkIo(comm);
        stkIo.set_bulk_data(bulkData);

        double timeStep = .1;
        double totalTime = 10;
        run_game_of_life(bulkData, stkIo, meshName, elements,
                              elemIdsToActivate, lifeField, activeNeighborField,
                             totalTime, timeStep);
    }
}
TEST(GameofLife, DISABLED_ParallelPartTriangleMesh)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);

    {
        stk::mesh::MetaData meta(2);
        stk::mesh::BulkData bulkData(meta, comm);
        stk::mesh::Part* triPart = &meta.declare_part_with_topology("Triangle_Part", stk::topology::TRIANGLE_3);
        stk::mesh::Part& activePart = meta.declare_part_with_topology("Active", stk::topology::TRIANGLE_3);

        stk::io::put_io_part_attribute(meta.universal_part());
        //stk::io::put_io_part_attribute(*triPart);
        stk::io::put_io_part_attribute(activePart);

        ScalarIntField& activeNeighborField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Status");
        int initVal = 0;
        stk::mesh::put_field(activeNeighborField, meta.universal_part(), &initVal);
        ScalarIntField& lifeField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Field");
        stk::mesh::put_field(lifeField, meta.universal_part(), &initVal);

        stk::mesh::Field<double,stk::mesh::Cartesian2d>& nodeCoord =
                meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian2d>>(stk::topology::NODE_RANK, "coordinates");
        stk::mesh::put_field(nodeCoord, meta.universal_part(), 2);

        meta.commit();

        const unsigned elemRowsPerProc = 12;
        const size_t nodesPerRow = 101, nodesPerCol = (elemRowsPerProc*numProcs)+1;
        const size_t numElems = 2*(nodesPerRow-1)*(nodesPerCol-1);
        unsigned numElemsPerProc = elemRowsPerProc*2*(nodesPerRow-1);
        unsigned elemOffset = numElemsPerProc*procRank;

        stk::mesh::EntityIdVector elemIds;
        elemIds.resize(numElems);
        for (size_t index = 0; index < numElems; index++)
        {
            elemIds[index] = index+1;
        }

        //bulkData.generate_new_ids(stk::topology::ELEM_RANK, numElemsPerProc, elemIds);
        stk::mesh::EntityIdVector nodeIds[numElems];
        compute_triangle_node_id_array(nodeIds, procRank, nodesPerRow, elemRowsPerProc);

        bulkData.modification_begin();

        for (size_t index = 0; index < numElemsPerProc; index++)
            declare_element(bulkData, *triPart, elemIds[index+elemOffset], nodeIds[index+elemOffset]);

        if (procRank != 0)
        {
            unsigned bottomNodeOffset = elemRowsPerProc*nodesPerRow*procRank;
            for (unsigned index = 1; index <= nodesPerRow; index++)
            {
                stk::mesh::Entity bottomNode = bulkData.get_entity(stk::topology::NODE_RANK, index + bottomNodeOffset);
                bulkData.add_node_sharing(bottomNode, procRank-1);
            }
            }
            if (numProcs-1 != procRank)
            {
                unsigned topNodeOffset = elemRowsPerProc*nodesPerRow*(procRank+1);
                for (unsigned index = 1; index <= nodesPerRow; index++)
                {
                    stk::mesh::Entity topNode = bulkData.get_entity(stk::topology::NODE_RANK, index + topNodeOffset);
                    bulkData.add_node_sharing(topNode, procRank+1);
                }
            }

            stk::mesh::EntityVector elems;
            stk::mesh::get_entities(bulkData, stk::topology::ELEM_RANK, elems);

            stk::mesh::EntityIdVector elemIds1;
            for (auto R: elems)
                elemIds1.push_back(bulkData.identifier(R));

            stk::mesh::EntityIdVector localElemIds1;
            for (auto R: elems)
                localElemIds1.push_back(bulkData.local_id(R));


            stk::mesh::EntityVector nodes;
            stk::mesh::get_entities(bulkData, stk::topology::NODE_RANK, nodes);
            size_t numNodes = nodes.size();
            for (size_t index = 0; index < numNodes; index++)
            {
                if (bulkData.is_valid(nodes[index]))
                {
                    double* const coord = stk::mesh::field_data(nodeCoord, nodes[index]);
                    coord[0] = index%nodesPerRow;
                    coord[1] = index/nodesPerRow + elemRowsPerProc*procRank;
                }
            }
            bulkData.modification_end();


            stk::mesh::EntityVector elements;
            stk::mesh::get_selected_entities(meta.locally_owned_part(), bulkData.buckets(stk::topology::ELEM_RANK), elements);

            stk::mesh::PartVector activeVector = {&activePart};
            stk::mesh::PartVector emptyVector;

            bulkData.modification_begin();
            for (size_t index = 1; index <= elements.size(); index++)
            {
                if (bulkData.is_valid(elements[index]))
                    if ((!(index%19)) || (!(index%13)) || (!index&11))
                        bulkData.change_entity_parts(elements[index], activeVector, emptyVector);
            }
            bulkData.modification_end();

            const double totalTime = 1000;
            const double timeStep = 1;

            std::string meshName = "ParallelPartTriangleMesh.e";
            stk::io::StkMeshIoBroker stkIo(comm);
            stkIo.set_bulk_data(bulkData);
            size_t fh = stkIo.create_output_mesh(meshName, stk::io::WRITE_RESULTS);
            stkIo.add_field(fh, lifeField);
            stkIo.write_output_mesh(fh);


            for (double time = 0; time < totalTime; time += timeStep)
            {
                stkIo.begin_output_step(fh, time);
                stkIo.write_defined_output_fields(fh);
                stkIo.end_output_step(fh);
                std::set<stk::mesh::Entity> activeConnectedEntities;
                for (size_t elemIndex = 0, numElems = elements.size(); elemIndex < numElems; elemIndex++)
                {
                    activeConnectedEntities.clear();
                    if (!bulkData.is_valid(elements[elemIndex]))
                        continue;
                    int numNodes = bulkData.num_nodes(elements[elemIndex]);
                    const stk::mesh::Entity* elemNode = bulkData.begin_nodes(elements[elemIndex]);
                    for (int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
                    {
                        stk::mesh::Entity node = elemNode[nodeIndex];
                        int numElems = bulkData.num_elements(node);
                        const stk::mesh::Entity* nodeElem = bulkData.begin_elements(node);
                        for (int index = 0; index < numElems; index++)
                            if (nodeElem[index] != elements[elemIndex])
                                if (bulkData.bucket(nodeElem[index]).member(activePart))
                                    activeConnectedEntities.insert(nodeElem[index]);
                    }
                    int* elem_value = stk::mesh::field_data(activeNeighborField, elements[elemIndex]);
                    *elem_value = activeConnectedEntities.size();
                }
                bulkData.modification_begin();
                for (size_t elemIndex = 0, numElems = elements.size(); elemIndex < numElems; elemIndex++)
                {
                    int* numNeighbors = stk::mesh::field_data(activeNeighborField, elements[elemIndex]);
                    switch (*numNeighbors)
                    {
                        case 2:
                        case 7:
                            break;
                        case 3:
                            bulkData.change_entity_parts(elements[elemIndex], activeVector, emptyVector);
                            *stk::mesh::field_data(lifeField, elements[elemIndex]) = 1;
                            break;
                        default:
                            bulkData.change_entity_parts(elements[elemIndex], emptyVector, activeVector);
                            *stk::mesh::field_data(lifeField, elements[elemIndex]) = 0;
                    }
                }
                bulkData.modification_end();
            }
        }
}
TEST(GameofLife, DISABLED_LazyHexMesh)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    //    const int numProcs = stk::parallel_machine_size(comm);
    //    const int procNum = stk::parallel_machine_rank(comm);

    stk::mesh::MetaData meta(3);
    stk::mesh::BulkData bulk(meta, comm);


    ScalarIntField& activeNeighborField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Active Neighbors");
    ScalarIntField& lifeField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Field");
    int initNeighborVal = 0;
    int initLifeVal = 0;
    stk::mesh::put_field(activeNeighborField, meta.universal_part(), &initNeighborVal);
    stk::mesh::put_field(lifeField, meta.universal_part(), &initLifeVal);

    std::string meshSpec = getGeneratedMeshString(4, 4, 8);
    stk::io::StkMeshIoBroker stkIo(comm);

    stkIo.set_bulk_data(bulk);
    stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    std::string meshName = "LazyHexMesh.e";

    size_t fh = stkIo.create_output_mesh(meshName, stk::io::WRITE_RESULTS);
    stkIo.write_output_mesh(fh);

//    stk::mesh::EntityVector elements;
//    get_entities(bulk, stk::topology::ELEM_RANK, elements);
//    stk::mesh::EntityIdVector elemIdsToActivate;
//    for (size_t index = 1; index <= elements.size(); index++)
//        for (unsigned index = 1; index <= elements.size(); index++)
//                if ((!(index%9)) || (!(index%5)) || (!index&11))
//                    elemIdsToActivate.push_back(index);
//    const double totalTime = 100;
//    const double timeStep = 1;

//    run_game_of_life(bulk, stkIo, meshName, elements,
//                          elemIdsToActivate, lifeField, activeNeighborField,
//                         totalTime, timeStep);
}
TEST(GameofLife, DISABLED_TryhardHexMesh)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    const int numProcs = stk::parallel_machine_size(comm);
    const int procNum = stk::parallel_machine_rank(comm);

    stk::mesh::MetaData meta(3);
    stk::mesh::BulkData bulk(meta, comm);
    stk::io::StkMeshIoBroker stkIo(comm);

    stk::mesh::Part* hexPart = &meta.declare_part_with_topology("Hex_Part", stk::topology::HEX_8);
    stk::io::put_io_part_attribute(meta.universal_part());
    stk::io::put_io_part_attribute(*hexPart);

    stk::mesh::Field<double,stk::mesh::Cartesian>& nodeCoord =
            meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian>>(stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field(nodeCoord, meta.universal_part(), 3);

    ScalarIntField& activeNeighborField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Active Neighbors");
    ScalarIntField& lifeField = meta.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Field");
    int initNeighborVal = 0;
    int initLifeVal = 0;
    stk::mesh::put_field(activeNeighborField, meta.universal_part(), &initNeighborVal);
    stk::mesh::put_field(lifeField, meta.universal_part(), &initLifeVal);
    meta.commit();

    const unsigned nodeWidth = 9;
    const unsigned nodeHeight = 9;
    const unsigned elemWidth = nodeWidth-1;
    const unsigned elemHeight = nodeHeight-1;
    const unsigned nodesPerSlice = nodeWidth*nodeHeight;
    const unsigned slicesPerProc = 1;
    const unsigned elemsPerSlice = elemWidth*elemHeight;
    const unsigned elemsPerProc = slicesPerProc*elemsPerSlice;
    const unsigned nodesPerProc = slicesPerProc*nodesPerSlice;
    const unsigned procElemOffset = elemsPerProc*procNum;
    const unsigned procNodeOffset = nodesPerSlice*slicesPerProc*procNum;

    stk::mesh::EntityIdVector elemIds(elemsPerProc);
    std::vector<stk::mesh::EntityIdVector> elemNodeIds(elemsPerProc);

    unsigned index;
    for (index = 0; index < elemsPerProc; index++)
    {
        elemIds[index] = index + procElemOffset + 1;
    }

    for (index = 0; index < elemsPerSlice; index++)
    {
        unsigned rowOffset = (index)/elemWidth;
        unsigned totalOffset = procNodeOffset+rowOffset+index;

        elemNodeIds[index].resize(8);
        elemNodeIds[index] = {totalOffset+1, totalOffset+2, totalOffset+nodesPerSlice+2, totalOffset+nodesPerSlice+1,
                       totalOffset+1+nodeWidth, totalOffset+2+nodeWidth, totalOffset+nodesPerSlice+2+nodeWidth,totalOffset+nodesPerSlice+1+nodeWidth};

    }
    for (; index < elemsPerProc; index++)
    {
        elemNodeIds[index].resize(8);
        for (unsigned nodeIndex = 0; nodeIndex<8; nodeIndex++)
            elemNodeIds[index][nodeIndex] = elemNodeIds[index%elemsPerSlice][nodeIndex] + index/elemsPerSlice*nodesPerSlice;
    }

    bulk.modification_begin();

    for (index = 0; index < elemsPerProc; index++)
        declare_element(bulk, *hexPart, elemIds[index], elemNodeIds[index]);

    stk::mesh::EntityVector nodes;
    stk::mesh::get_entities(bulk, stk::topology::NODE_RANK, nodes);

    if (numProcs-1 != procNum)
    {
        unsigned backNodeOffset = nodesPerProc*(procNum+1);
        for (unsigned nodeIndex = 1; nodeIndex <= nodesPerSlice; nodeIndex++)
        {
            stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK,
                                                     nodeIndex + backNodeOffset);
            bulk.add_node_sharing(node, procNum+1);
        }
    }
    if (0 != procNum)
    {
        unsigned frontNodeOffset = nodesPerProc*procNum;
        for (unsigned nodeIndex = 1; nodeIndex <= nodesPerSlice; nodeIndex++)
        {
            stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK,
                                                     nodeIndex + frontNodeOffset);
            bulk.add_node_sharing(node, procNum-1);
        }
    }
    bulk.modification_end();

    for (index = 0; index < nodes.size(); index++)
    {
        double* const coord = stk::mesh::field_data(nodeCoord, nodes[index]);
        coord[0]=index%nodesPerSlice%nodeWidth;
        coord[1]=index%nodesPerSlice/nodeWidth;
        coord[2]=index/nodesPerSlice + procNum*slicesPerProc;
    }

    std::string meshName = "TryhardHexMesh.e";

    stkIo.set_bulk_data(bulk);
    size_t fh = stkIo.create_output_mesh(meshName, stk::io::WRITE_RESULTS);
    stkIo.write_output_mesh(fh);
}

//Bitmap stuff
void convert_normal_char_vector_into_2D_vector(const unsigned width, const unsigned height,
                                      const std::vector<unsigned char>& charVector,
                                      std::vector<std::vector<unsigned>>& intVector)
{
    intVector.resize(height);
    unsigned olderIndex = 0;
    for (unsigned rowIndex = 0; rowIndex < height; rowIndex++)
    {
        intVector[rowIndex].resize(width);
        for (unsigned colIndex = 0; colIndex < width; colIndex ++)
        {
            int value = 0;
            for (int counter = 3; counter >= 0; counter--, olderIndex++)
                value |= (charVector[olderIndex]<<(counter*8));
            intVector[rowIndex][colIndex]=value;
        }
    }
}
void fill_element_id_vector_from_2D_int_vector(unsigned width, unsigned height,
                                               const std::vector<std::vector<unsigned>>& intVector,
                                               stk::mesh::EntityIdVector& elemIds)

{
    unsigned id = 1;
    for (int rowIndex = height-1; rowIndex >= 0; rowIndex--)
    {
        for (unsigned colIndex = 0; colIndex < width; colIndex++, id++)
            if (0x000000ff == intVector[rowIndex][colIndex])
                elemIds.push_back(id);
    }
}
TEST(Bitmap, DISABLED_LoadingBitmap1)
{
    const char* filename = "P94s2.png";
    std::vector<unsigned char> image;
    unsigned width, height;
    lodepng::decode(image, width, height, filename);
    unsigned index = 0;
    for (unsigned i = 0; i < height; i++)
    {
        for (unsigned j = 0; j < width; j++)
        {
            for (unsigned k = 0; k < 4; k++, index++)
                printf("%2x", image[index]);
            printf(" ");
        }
        printf("\n");
    }
}
TEST(Bitmap, DISABLED_LoadingBitmap2)
{
    const char* filename = "P94s2.png";
    std::vector<unsigned char> image;
    unsigned width, height;
    lodepng::decode(image, width, height, filename);

    std::vector<std::vector<unsigned>> intVector;
    convert_normal_char_vector_into_2D_vector(width, height, image, intVector);

    for (unsigned i = 0; i < height; i++)
    {
        for (unsigned j = 0; j < width; j++)
            printf("%8x ", intVector[i][j]);
        printf("\n");
    }
}
TEST(Bitmap, TuringMachine)
{
    const char* filename = "turing.png";
    std::vector<unsigned char> bits;
    unsigned width, height;
    lodepng::decode(bits, width, height, filename);

    std::vector<std::vector<unsigned>> intVector;
    convert_normal_char_vector_into_2D_vector(width, height, bits, intVector);

    stk::ParallelMachine comm = MPI_COMM_WORLD;
    unsigned rowsPerProc = (width/stk::parallel_machine_size(comm))+1;
    std::string meshName = "TuringMachine.e";

    stk::mesh::EntityIdVector elemIdsToActivate;
    fill_element_id_vector_from_2D_int_vector(width, height, intVector, elemIdsToActivate);

    QuadMesh* Mesh = new QuadMesh(comm, width, rowsPerProc, meshName);
    Mesh->fill_mesh();
    Mesh->write_output_mesh();
    Mesh->activate_following_element_ids_using_fields(elemIdsToActivate);
    Mesh->run_field_game_of_life(500);
}
TEST(Bitmap, DISABLED_GunsGameofLife)
{
    const char* filename = "P44guns.png";
    std::vector<unsigned char> bits;
    unsigned width, height;
    lodepng::decode(bits, width, height, filename);

    std::vector<std::vector<unsigned>> intVector;
    convert_normal_char_vector_into_2D_vector(width, height, bits, intVector);

    stk::ParallelMachine comm = MPI_COMM_WORLD;
    unsigned rowsPerProc = (height/stk::parallel_machine_size(comm))+1;
    std::string meshName = "Guns.e";

    stk::mesh::EntityIdVector elemIdsToActivate;
    fill_element_id_vector_from_2D_int_vector(width, height, intVector, elemIdsToActivate);

    QuadMesh* Mesh = new QuadMesh(comm, width, rowsPerProc, meshName);
    Mesh->fill_mesh();
    Mesh->write_output_mesh();
    Mesh->activate_following_element_ids_using_parts(elemIdsToActivate);
    // consider load balancing here, maybe?
    Mesh->run_part_game_of_life(50);
}
TEST(Bitmap, DISABLED_LoadBitmap3)
{
    const char* filename = "Boss.png";
    std::vector<unsigned char> bits;
    unsigned width, height;
    lodepng::decode(bits, width, height, filename);

    std::vector<std::vector<unsigned>> fakeIntVector;
    convert_normal_char_vector_into_2D_vector(width, height, bits, fakeIntVector);

    std::vector<std::vector<unsigned>> realIntVector;

    unsigned numVerticalLines = 0;
    unsigned numHorizontalLines = 0;

    for (unsigned i = 0; i < width; i++)
        if (0xc6c6c6ff == fakeIntVector[1][i])
            numVerticalLines++;

    for (unsigned i = 0; i < height; i++)
        if (0xc6c6c6ff == fakeIntVector[i][1])
            numHorizontalLines++;

    unsigned actualWidth = numVerticalLines-1;
    unsigned actualHeight = numHorizontalLines-1;
    unsigned squareWidth = (width-numVerticalLines)/actualWidth;
    unsigned squareHeight = (height-numHorizontalLines)/actualHeight;

    realIntVector.resize(actualHeight);
    for (unsigned verticalIndex = 0; verticalIndex < actualHeight; verticalIndex++)
    {
        unsigned verticalOffset = squareHeight*verticalIndex + verticalIndex + 1;
        realIntVector[verticalIndex].resize(actualWidth);
        for (unsigned horizontalIndex = 0; horizontalIndex < actualWidth; horizontalIndex++)
        {
            unsigned horizontalOffset = squareWidth*horizontalIndex + horizontalIndex + 1;
            realIntVector[verticalIndex][horizontalIndex] =
                    fakeIntVector[verticalOffset][horizontalOffset];
        }
    }
    for (unsigned i = 0; i < actualHeight; i++)
    {
        for (unsigned j = 0; j < actualWidth; j++)
            printf("%8d ", realIntVector[i][j]);
        printf("\n");
    }
    stk::mesh::EntityIdVector elemIds;
    fill_element_id_vector_from_2D_int_vector(actualWidth, actualHeight, realIntVector,elemIds);
    for (auto R : elemIds)
        printf("%lu ", R);

    stk::ParallelMachine comm = MPI_COMM_WORLD;
    unsigned rowsPerProc = (actualHeight/stk::parallel_machine_size(comm))+1;
    std::string meshName = "Boss.e";

    QuadMesh* Mesh = new QuadMesh(comm, actualWidth, rowsPerProc, meshName);
    Mesh->fill_mesh();
    Mesh->write_output_mesh();
    Mesh->activate_following_element_ids_using_fields(elemIds);
    Mesh->run_field_game_of_life(10);
}

TEST(Trash, Test)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int procRank = stk::parallel_machine_rank(comm);
    stk::mesh::MetaData meta(3);
    stk::mesh::BulkData bulkData(meta, comm, stk::mesh::BulkData::NO_AUTO_AURA);
    std::string meshSpec = "generated:16x16x16";

    stk::io::StkMeshIoBroker stkIo(comm);
    stkIo.set_bulk_data(bulkData);
    stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
    stkIo.create_input_mesh();

    int initialVal = 0;
    ScalarIntField& lifeField =
            meta.declare_field<ScalarIntField>(stk::topology::ELEM_RANK, "Life");
    ScalarIntField& neighborField =
            meta.declare_field<ScalarIntField>(stk::topology::ELEM_RANK, "Friend");
    stk::mesh::put_field(lifeField, meta.universal_part(), &initialVal);
    stk::mesh::put_field(neighborField, meta.universal_part(), &initialVal);

    stkIo.populate_bulk_data();

    stk::mesh::EntityVector elements;
    stk::mesh::EntityVector nodes;
    stk::mesh::get_entities(bulkData, stk::topology::ELEM_RANK, elements);
    stk::mesh::get_entities(bulkData, stk::topology::NODE_RANK, nodes);
    unsigned numElems = elements.size();

    std::unordered_map<stk::mesh::Entity, std::unordered_set<stk::mesh::Entity,
    std::hash<stk::mesh::Entity>>, std::hash<stk::mesh::Entity>>
    localElementToLocalNeighborElements;

    std::unordered_map<stk::mesh::EntityKey, std::unordered_set<stk::mesh::EntityKey,
    std::hash<stk::mesh::EntityKey>>, std::hash<stk::mesh::EntityKey>>
    remoteElementKeyToLocalNodeKeys;

    std::unordered_set<stk::mesh::EntityKey, std::hash<stk::mesh::EntityKey>>
    remoteElementKeys;

    std::unordered_map<stk::mesh::EntityKey, int, std::hash<stk::mesh::EntityKey>>
    remoteElementKeyToProcessor;

    std::unordered_map<stk::mesh::EntityKey, std::unordered_set<stk::mesh::Entity,
    std::hash<stk::mesh::Entity>>, std::hash<stk::mesh::EntityKey>>
    remoteElementKeyToLocalNeighborElements;

    //local elements
    for (unsigned index = 0; index < numElems; index++)
    {
        stk::mesh::Entity elem = elements[index];
        unsigned numNodes = bulkData.num_nodes(elem);
        const stk::mesh::Entity* elemNodes = bulkData.begin_nodes(elem);
        for (unsigned nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
        {
            stk::mesh::Entity node = elemNodes[nodeIndex];
            unsigned numElems = bulkData.num_elements(node);
            const stk::mesh::Entity* nodeElems = bulkData.begin_elements(node);
            for (unsigned elemIndex = 0; elemIndex < numElems; elemIndex++)
                if (elem != nodeElems[elemIndex])
                    localElementToLocalNeighborElements[elem].insert(nodeElems[elemIndex]);
        }
    }

    //pack local elem keys with shared node keys to remote processors
    stk::CommSparse send(bulkData.parallel());
    for (int phase = 0; phase < 2; phase++)
    {
        for (stk::mesh::Entity elem : elements)
        {
            stk::mesh::EntityKey elemKey = bulkData.entity_key(elem);
            unsigned numNodes = bulkData.num_nodes(elem);
            const stk::mesh::Entity* nodeBegin = bulkData.begin_nodes(elem);
            for (unsigned nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
            {
                //parallel elements
                std::vector<int> sharingProcs;
                stk::mesh::EntityKey nodeKey = bulkData.entity_key(nodeBegin[nodeIndex]);
                bulkData.comm_shared_procs(nodeKey, sharingProcs);
                for (int procNum : sharingProcs)
                {
                    send.send_buffer(procNum).pack<int>(procRank);
                    send.send_buffer(procNum).pack<stk::mesh::EntityKey>(elemKey);
                    send.send_buffer(procNum).pack<stk::mesh::EntityKey>(nodeKey);
                }
            }
        }
        if (0 == phase)
            send.allocate_buffers();
        else
            send.communicate();
    }

    //unpack remote elem keys with shared node keys and add to the map
    for (int procRank = 0; procRank < numProcs; procRank++)
    {

        stk::CommBuffer& sendBuf = send.recv_buffer(procRank);
        while (sendBuf.remaining())
        {
            int procNum;
            stk::mesh::EntityKey elemKey;
            stk::mesh::EntityKey nodeKey;
            sendBuf.unpack<int>(procNum);
            sendBuf.unpack<stk::mesh::EntityKey>(elemKey);
            sendBuf.unpack<stk::mesh::EntityKey>(nodeKey);
            remoteElementKeys.insert(elemKey);
            remoteElementKeyToProcessor[elemKey]=procNum;
            remoteElementKeyToLocalNodeKeys[elemKey].insert(nodeKey);
        }
    }

    // set up map of remote elem keys to local elements
    for (stk::mesh::EntityKey remoteElemKey : remoteElementKeys)
    {
        for (stk::mesh::EntityKey localNodeKey : remoteElementKeyToLocalNodeKeys[remoteElemKey])
        {
            stk::mesh::Entity node = bulkData.get_entity(localNodeKey);
            unsigned numElems = bulkData.num_elements(node);
            const stk::mesh::Entity* nodeElem = bulkData.begin_elements(node);
            for (unsigned elemIndex = 0; elemIndex < numElems; elemIndex++)
                remoteElementKeyToLocalNeighborElements[remoteElemKey].insert(nodeElem[elemIndex]);
        }
    }

    // start to run game of life?
    for (unsigned index = 0; index < numElems; index++)
        if (!(index%5) || !(index%7))
            *stk::mesh::field_data(lifeField, elements[index]) = 1;

    size_t fh = stkIo.create_output_mesh("Trash.e", stk::io::WRITE_RESULTS);
    stkIo.add_field(fh, lifeField);
    stkIo.add_field(fh, neighborField);
    stkIo.write_output_mesh(fh);

    for (int time = 0; time <= 100; time++)
    {
        stkIo.begin_output_step(fh, time);
        stkIo.write_defined_output_fields(fh);
        stkIo.end_output_step(fh);

        //update neighbor values with local elements
        for (stk::mesh::Entity localElem : elements)
        {
            int* neighborVal = stk::mesh::field_data(neighborField, localElem);
            *neighborVal = 0;
            for (stk::mesh::Entity elemElem : localElementToLocalNeighborElements[localElem])
                if (*stk::mesh::field_data(lifeField, elemElem))
                    (*neighborVal)++;
        }

        //update neighbor values with remote elements
        stk::CommSparse sombrero(bulkData.parallel());
        for (int phase = 0 ; phase < 2; phase++)
        {
            for (stk::mesh::EntityKey remoteElemKey : remoteElementKeys)
            {
                int neighborVal = 0;
                int procNum = remoteElementKeyToProcessor[remoteElemKey];

                for (stk::mesh::Entity localElem :
                        remoteElementKeyToLocalNeighborElements[remoteElemKey])
                    if (*stk::mesh::field_data(lifeField, localElem))
                        neighborVal++;

                sombrero.send_buffer(procNum).pack<stk::mesh::EntityKey>(remoteElemKey);
                sombrero.send_buffer(procNum).pack<int>(neighborVal);
            }
            if (0 == phase)
                sombrero.allocate_buffers();
            else if (1 == phase)
                sombrero.communicate();
        }

        for (int procNum = 0; procNum < numProcs; procNum++)
        {
            stk::CommBuffer& hat = sombrero.recv_buffer(procNum);
            while(hat.remaining())
            {
                stk::mesh::EntityKey localElemKey;
                int neighborVal;
                hat.unpack<stk::mesh::EntityKey>(localElemKey);
                hat.unpack<int>(neighborVal);
                *stk::mesh::field_data(neighborField, bulkData.get_entity(localElemKey))
                += neighborVal;
            }
        }
        for (stk::mesh::Entity elem : elements)
        {
            int* lifeVal = stk::mesh::field_data(lifeField, elem);
            switch (*stk::mesh::field_data(neighborField, elem))
            {
                case 4:
                    break;
                case 5:
                    *lifeVal = 1;
                    break;
                default:
                    *lifeVal = 0;
                    break;
            }

        }
    }
}

TEST(Trash, PartTestTrash)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int myProc = stk::parallel_machine_rank(comm);
    stk::mesh::MetaData meta(3);
    stk::mesh::BulkData bulkData(meta, comm, stk::mesh::BulkData::NO_AUTO_AURA);
    std::string meshSpec = "generated:16x16x16";

    stk::io::StkMeshIoBroker stkIo(comm);
    stkIo.set_bulk_data(bulkData);
    stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
    stkIo.create_input_mesh();

    int initialVal = 0;
    ScalarIntField& lifeField =
            meta.declare_field<ScalarIntField>(stk::topology::ELEM_RANK, "Life");
    ScalarIntField& neighborField =
            meta.declare_field<ScalarIntField>(stk::topology::ELEM_RANK, "Friend");
    stk::mesh::put_field(lifeField, meta.universal_part(), &initialVal);
    stk::mesh::put_field(neighborField, meta.universal_part(), &initialVal);

    stk::mesh::Part& activePart = meta.declare_part_with_topology("active", stk::topology::HEX_8);

    stkIo.populate_bulk_data();

    stk::mesh::EntityVector elements;
    stk::mesh::EntityVector nodes;
    stk::mesh::get_entities(bulkData, stk::topology::ELEM_RANK, elements);
    stk::mesh::get_entities(bulkData, stk::topology::NODE_RANK, nodes);
    unsigned numElems = elements.size();

    std::unordered_map<stk::mesh::Entity, std::unordered_set<stk::mesh::Entity,
    std::hash<stk::mesh::Entity>>, std::hash<stk::mesh::Entity>>
    localElementToLocalNeighborElements;

    std::unordered_set<stk::mesh::EntityKey, std::hash<stk::mesh::EntityKey>>
    remoteElementKeys;

    std::unordered_map<stk::mesh::EntityKey, std::unordered_set<stk::mesh::EntityKey,
    std::hash<stk::mesh::EntityKey>>, std::hash<stk::mesh::EntityKey>>
    remoteElementKeyToLocalNodeKeys;

    std::unordered_map<stk::mesh::EntityKey, int, std::hash<stk::mesh::EntityKey>>
    remoteElementKeyToOwningProcessor;

    std::unordered_map<stk::mesh::EntityKey, std::unordered_set<stk::mesh::Entity,
    std::hash<stk::mesh::Entity>>, std::hash<stk::mesh::EntityKey>>
    remoteElementKeyToLocalNeighborElements;

    std::unordered_map<stk::mesh::Entity, std::unordered_set<stk::mesh::EntityKey,
    std::hash<stk::mesh::EntityKey>>, std::hash<stk::mesh::Entity>>
    localElementToRemoteElementKeys;

    std::vector<stk::mesh::Entity>
    localActiveElements;

    std::unordered_set<stk::mesh::Entity, std::hash<stk::mesh::Entity>>
    localElementsToVisit;

    std::unordered_set<stk::mesh::EntityKey, std::hash<stk::mesh::EntityKey>>
    remoteElementKeysToVisit;

    //local elements
    for (unsigned index = 0; index < numElems; index++)
    {
        stk::mesh::Entity elem = elements[index];
        unsigned numNodes = bulkData.num_nodes(elem);
        const stk::mesh::Entity* elemNodes = bulkData.begin_nodes(elem);

        for (unsigned nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
        {
            stk::mesh::Entity node = elemNodes[nodeIndex];
            unsigned numElems = bulkData.num_elements(node);
            const stk::mesh::Entity* nodeElems = bulkData.begin_elements(node);

            for (unsigned elemIndex = 0; elemIndex < numElems; elemIndex++)
                if (elem != nodeElems[elemIndex])
                    localElementToLocalNeighborElements[elem].insert(nodeElems[elemIndex]);
        }
    }

    //pack local elem keys with shared node keys to remote processors
    stk::CommSparse send(bulkData.parallel());
    for (int phase = 0; phase < 2; phase++)
    {
        std::unordered_map<int,std::unordered_set<stk::mesh::EntityKey, std::hash<stk::mesh::EntityKey>>>
        remoteProcessorNumberToSharedNodes;

        for (stk::mesh::Entity elem : elements)
        {
            remoteProcessorNumberToSharedNodes.clear();
            stk::mesh::EntityKey elemKey = bulkData.entity_key(elem);
            unsigned numNodes = bulkData.num_nodes(elem);
            const stk::mesh::Entity* nodeBegin = bulkData.begin_nodes(elem);
            for (unsigned nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
            {
                //parallel elements
                std::vector<int> sharingProcs;
                stk::mesh::EntityKey nodeKey = bulkData.entity_key(nodeBegin[nodeIndex]);
                bulkData.comm_shared_procs(nodeKey, sharingProcs);

                for (int procNum : sharingProcs)
                    remoteProcessorNumberToSharedNodes[procNum].insert(nodeKey);
            }

            for (std::pair< const int,std::unordered_set<stk::mesh::EntityKey,
                    std::hash<stk::mesh::EntityKey>>>& pair : remoteProcessorNumberToSharedNodes)
            {
                int remoteProc = pair.first;
                send.send_buffer(remoteProc).pack<int>(myProc);
                send.send_buffer(remoteProc).pack<stk::mesh::EntityKey>(elemKey);
                send.send_buffer(remoteProc).pack<size_t>(pair.second.size());

                for (stk::mesh::EntityKey nodeKey : pair.second)
                    send.send_buffer(remoteProc).pack<stk::mesh::EntityKey>(nodeKey);
            }
        }
        if (0 == phase)
            send.allocate_buffers();
        else
            send.communicate();
    }

    //unpack remote elem keys with shared node keys and add to the map
    for (int procRank = 0; procRank < numProcs; procRank++)
    {

        stk::CommBuffer& sendBuf = send.recv_buffer(procRank);
        while (sendBuf.remaining())
        {
            int procNum;
            size_t numNodes;
            stk::mesh::EntityKey elemKey;
            stk::mesh::EntityKey nodeKey;
            sendBuf.unpack<int>(procNum);
            sendBuf.unpack<stk::mesh::EntityKey>(elemKey);
            sendBuf.unpack<size_t>(numNodes);
            remoteElementKeys.insert(elemKey);
            remoteElementKeyToOwningProcessor[elemKey]=procNum;
            for (unsigned nodeNum = 0; nodeNum < numNodes; nodeNum++)
            {
                sendBuf.unpack<stk::mesh::EntityKey>(nodeKey);
                remoteElementKeyToLocalNodeKeys[elemKey].insert(nodeKey);
            }
        }
    }

    // set up map of remote elem keys to local elements
    for (stk::mesh::EntityKey remoteElemKey : remoteElementKeys)
    {
        for (stk::mesh::EntityKey localNodeKey : remoteElementKeyToLocalNodeKeys[remoteElemKey])
        {
            stk::mesh::Entity node = bulkData.get_entity(localNodeKey);
            unsigned numElems = bulkData.num_elements(node);
            const stk::mesh::Entity* nodeElem = bulkData.begin_elements(node);
            for (unsigned elemIndex = 0; elemIndex < numElems; elemIndex++)
                remoteElementKeyToLocalNeighborElements[remoteElemKey].insert(nodeElem[elemIndex]);
        }
    }

    // reply with the local element keys of remote element keys to create local elements to remote element key map
    stk::CommSparse reply(bulkData.parallel());
    for (int phase = 0; phase < 2; phase++)
    {
        for (stk::mesh::EntityKey remoteElemKey : remoteElementKeys)
        {
            int procNum = remoteElementKeyToOwningProcessor[remoteElemKey];
            size_t numNeighbors = remoteElementKeyToLocalNeighborElements[remoteElemKey].size();
            reply.send_buffer(procNum).pack<stk::mesh::EntityKey>(remoteElemKey);
            reply.send_buffer(procNum).pack<size_t>(numNeighbors);
            for (stk::mesh::Entity localElem : remoteElementKeyToLocalNeighborElements[remoteElemKey])
            {
                stk::mesh::EntityKey localElemKey = bulkData.entity_key(localElem);
                reply.send_buffer(procNum).pack<stk::mesh::EntityKey>(localElemKey);
            }

        }
        if (0 == phase)
            reply.allocate_buffers();
        else
            reply.communicate();
    }

    //now get the local elements to store the neighboring remote elements mapping to their local elements
    for (int procRank = 0; procRank < numProcs; procRank++)
    {
        stk::CommBuffer& replyBuf = reply.recv_buffer(procRank);
        while (replyBuf.remaining())
        {
            size_t numNeighbors;
            stk::mesh::EntityKey localElemKey;
            stk::mesh::EntityKey remoteElemKey;
            replyBuf.unpack<stk::mesh::EntityKey>(localElemKey);
            replyBuf.unpack<size_t>(numNeighbors);
            stk::mesh::Entity localElem = bulkData.get_entity(localElemKey);
            for (unsigned neighborNum = 0; neighborNum < numNeighbors; neighborNum++)
            {
                replyBuf.unpack<stk::mesh::EntityKey>(remoteElemKey);
                localElementToRemoteElementKeys[localElem].insert(remoteElemKey);
            }
        }
    }

    // start to run game of life?
    bulkData.modification_begin();
    unsigned lollipop;
    for (lollipop = 0; lollipop < numElems; lollipop++)
    {
        if (!(lollipop%5) || !(lollipop%7))
        {
            *stk::mesh::field_data(lifeField, elements[lollipop]) = 1;
            localActiveElements.push_back(elements[lollipop]);
            bulkData.change_entity_parts(elements[lollipop], {&activePart}, {});
        }
    }
    bulkData.modification_end();

    size_t fh = stkIo.create_output_mesh("PartTestTrash.e", stk::io::WRITE_RESULTS);
    stkIo.add_field(fh, lifeField);
    stkIo.add_field(fh, neighborField);
    stkIo.write_output_mesh(fh);

    for (int time = 0; time <= 10; time++)
    {
        //determine which elements to check
        localElementsToVisit.clear();
        remoteElementKeysToVisit.clear();

        //get local elements
//        bulkData.get_entities(stk::topology::ELEM_RANK, activePart, localActiveElements);

        //get local and remote elements to visit
        for (stk::mesh::Entity localElem : localActiveElements)
        {
            localElementsToVisit.insert(localElem);
            for (stk::mesh::Entity localElemElem : localElementToLocalNeighborElements[localElem])
                localElementsToVisit.insert(localElemElem);

            for (stk::mesh::EntityKey remoteElemKey : localElementToRemoteElementKeys[localElem])
                remoteElementKeysToVisit.insert(remoteElemKey);
        }

        //communicate said local elements
        stk::CommSparse sneaker(bulkData.parallel());
        for (int phase = 0; phase < 2; phase++)
        {
            for (stk::mesh::EntityKey remoteElemKey : remoteElementKeysToVisit)
            {
                int procNum = remoteElementKeyToOwningProcessor[remoteElemKey];
                sneaker.send_buffer(procNum).pack<stk::mesh::EntityKey>(remoteElemKey);
            }

            if (0 == phase)
                sneaker.allocate_buffers();
            else
                sneaker.communicate();
        }
        for (int procNum = 0; procNum < numProcs; procNum++)
        {
            stk::CommBuffer& shoe = sneaker.recv_buffer(procNum);

            while (shoe.remaining())
            {
                stk::mesh::EntityKey localElemKey;
                shoe.unpack<stk::mesh::EntityKey>(localElemKey);
                localElementsToVisit.insert(bulkData.get_entity(localElemKey));
            }
        }

        //update neighbor values with local elements
        for (stk::mesh::Entity localElem : localElementsToVisit)
        {
            int* neighborVal = stk::mesh::field_data(neighborField, localElem);
            *neighborVal = 0;
            for (stk::mesh::Entity localElemElem : localElementToLocalNeighborElements[localElem])
                if (*stk::mesh::field_data(lifeField, localElemElem))
                    (*neighborVal)++;
        }

        //update neighbor values with remote elements
        stk::CommSparse sombrero(bulkData.parallel());
        for (int phase = 0 ; phase < 2; phase++)
        {
            for (stk::mesh::EntityKey remoteElemKey : remoteElementKeysToVisit)
            {
                int neighborVal = 0;
                int procNum = remoteElementKeyToOwningProcessor[remoteElemKey];

                for (stk::mesh::Entity localElem :
                        remoteElementKeyToLocalNeighborElements[remoteElemKey])
                    if (*stk::mesh::field_data(lifeField, localElem))
                        neighborVal++;

                sombrero.send_buffer(procNum).pack<stk::mesh::EntityKey>(remoteElemKey);
                sombrero.send_buffer(procNum).pack<int>(neighborVal);
            }
            if (0 == phase)
                sombrero.allocate_buffers();
            else if (1 == phase)
                sombrero.communicate();
        }

        for (int procNum = 0; procNum < numProcs; procNum++)
        {
            stk::CommBuffer& hat = sombrero.recv_buffer(procNum);
            while(hat.remaining())
            {
                stk::mesh::EntityKey localElemKey;
                int neighborVal;
                hat.unpack<stk::mesh::EntityKey>(localElemKey);
                hat.unpack<int>(neighborVal);
                *stk::mesh::field_data(neighborField, bulkData.get_entity(localElemKey))
                += neighborVal;
            }
        }

        //write output whenever
        stkIo.begin_output_step(fh, time);
        stkIo.write_defined_output_fields(fh);
        stkIo.end_output_step(fh);

        localActiveElements.clear();

        bulkData.modification_begin();
        for (stk::mesh::Entity localElem : localElementsToVisit)
        {
            int* lifeVal = stk::mesh::field_data(lifeField, localElem);
            switch (*stk::mesh::field_data(neighborField, localElem))
            {
                case 4:
                    break;
                case 5:
                    bulkData.change_entity_parts(localElem, {&activePart}, {});
                    *lifeVal = 1;
                    break;
                default:
                    bulkData.change_entity_parts(localElem, {}, {&activePart});
                    *lifeVal = 0;
                    break;
            }
            if (*lifeVal)
                localActiveElements.push_back(localElem);
        }
        bulkData.modification_end();
    }
}

TEST(Trash, FieldTestTrash)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    int myProc = stk::parallel_machine_rank(comm);
    stk::mesh::MetaData meta(3);
    stk::mesh::BulkData bulkData(meta, comm, stk::mesh::BulkData::NO_AUTO_AURA);
    std::string meshSpec = "generated:16x16x16";

    stk::io::StkMeshIoBroker stkIo(comm);
    stkIo.set_bulk_data(bulkData);
    stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
    stkIo.create_input_mesh();

    int initialVal = 0;
    ScalarIntField& lifeField =
            meta.declare_field<ScalarIntField>(stk::topology::ELEM_RANK, "Life");
    ScalarIntField& neighborField =
            meta.declare_field<ScalarIntField>(stk::topology::ELEM_RANK, "Friend");
    stk::mesh::put_field(lifeField, meta.universal_part(), &initialVal);
    stk::mesh::put_field(neighborField, meta.universal_part(), &initialVal);

    stkIo.populate_bulk_data();

    stk::mesh::EntityVector elements;
    stk::mesh::EntityVector nodes;
    stk::mesh::get_entities(bulkData, stk::topology::ELEM_RANK, elements);
    stk::mesh::get_entities(bulkData, stk::topology::NODE_RANK, nodes);
    unsigned numElems = elements.size();

    std::unordered_map<stk::mesh::Entity, std::unordered_set<stk::mesh::Entity,
    std::hash<stk::mesh::Entity>>, std::hash<stk::mesh::Entity>>
    localElementToLocalNeighborElements;

    std::unordered_set<stk::mesh::EntityKey, std::hash<stk::mesh::EntityKey>>
    remoteElementKeys;

    std::unordered_map<stk::mesh::EntityKey, std::unordered_set<stk::mesh::EntityKey,
    std::hash<stk::mesh::EntityKey>>, std::hash<stk::mesh::EntityKey>>
    remoteElementKeyToLocalNodeKeys;

    std::unordered_map<stk::mesh::EntityKey, int, std::hash<stk::mesh::EntityKey>>
    remoteElementKeyToOwningProcessor;

    std::unordered_map<stk::mesh::EntityKey, std::unordered_set<stk::mesh::Entity,
    std::hash<stk::mesh::Entity>>, std::hash<stk::mesh::EntityKey>>
    remoteElementKeyToLocalNeighborElements;

    std::unordered_map<stk::mesh::Entity, std::unordered_set<stk::mesh::EntityKey,
    std::hash<stk::mesh::EntityKey>>, std::hash<stk::mesh::Entity>>
    localElementToRemoteElementKeys;

    std::vector<stk::mesh::Entity>
    localActiveElements;

    std::unordered_set<stk::mesh::Entity, std::hash<stk::mesh::Entity>>
    localElementsToVisit;

    std::unordered_set<stk::mesh::EntityKey, std::hash<stk::mesh::EntityKey>>
    remoteElementKeysToVisit;

    //local elements
    for (unsigned index = 0; index < numElems; index++)
    {
        stk::mesh::Entity elem = elements[index];
        unsigned numNodes = bulkData.num_nodes(elem);
        const stk::mesh::Entity* elemNodes = bulkData.begin_nodes(elem);

        for (unsigned nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
        {
            stk::mesh::Entity node = elemNodes[nodeIndex];
            unsigned numElems = bulkData.num_elements(node);
            const stk::mesh::Entity* nodeElems = bulkData.begin_elements(node);

            for (unsigned elemIndex = 0; elemIndex < numElems; elemIndex++)
                if (elem != nodeElems[elemIndex])
                    localElementToLocalNeighborElements[elem].insert(nodeElems[elemIndex]);
        }
    }

    //pack local elem keys with shared node keys to remote processors
    stk::CommSparse send(bulkData.parallel());
    for (int phase = 0; phase < 2; phase++)
    {
        std::unordered_map<int,std::unordered_set<stk::mesh::EntityKey, std::hash<stk::mesh::EntityKey>>>
        remoteProcessorNumberToSharedNodes;

        for (stk::mesh::Entity elem : elements)
        {
            remoteProcessorNumberToSharedNodes.clear();
            stk::mesh::EntityKey elemKey = bulkData.entity_key(elem);
            unsigned numNodes = bulkData.num_nodes(elem);
            const stk::mesh::Entity* nodeBegin = bulkData.begin_nodes(elem);
            for (unsigned nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
            {
                //parallel elements
                std::vector<int> sharingProcs;
                stk::mesh::EntityKey nodeKey = bulkData.entity_key(nodeBegin[nodeIndex]);
                bulkData.comm_shared_procs(nodeKey, sharingProcs);

                for (int procNum : sharingProcs)
                    remoteProcessorNumberToSharedNodes[procNum].insert(nodeKey);
            }

            for (std::pair< const int,std::unordered_set<stk::mesh::EntityKey,
                    std::hash<stk::mesh::EntityKey>>>& pair : remoteProcessorNumberToSharedNodes)
            {
                int remoteProc = pair.first;
                send.send_buffer(remoteProc).pack<int>(myProc);
                send.send_buffer(remoteProc).pack<stk::mesh::EntityKey>(elemKey);
                send.send_buffer(remoteProc).pack<size_t>(pair.second.size());

                for (stk::mesh::EntityKey nodeKey : pair.second)
                    send.send_buffer(remoteProc).pack<stk::mesh::EntityKey>(nodeKey);
            }
        }
        if (0 == phase)
            send.allocate_buffers();
        else
            send.communicate();
    }

    //unpack remote elem keys with shared node keys and add to the map
    for (int procRank = 0; procRank < numProcs; procRank++)
    {

        stk::CommBuffer& sendBuf = send.recv_buffer(procRank);
        while (sendBuf.remaining())
        {
            int procNum;
            size_t numNodes;
            stk::mesh::EntityKey elemKey;
            stk::mesh::EntityKey nodeKey;
            sendBuf.unpack<int>(procNum);
            sendBuf.unpack<stk::mesh::EntityKey>(elemKey);
            sendBuf.unpack<size_t>(numNodes);
            remoteElementKeys.insert(elemKey);
            remoteElementKeyToOwningProcessor[elemKey]=procNum;
            for (unsigned nodeNum = 0; nodeNum < numNodes; nodeNum++)
            {
                sendBuf.unpack<stk::mesh::EntityKey>(nodeKey);
                remoteElementKeyToLocalNodeKeys[elemKey].insert(nodeKey);
            }
        }
    }

    // set up map of remote elem keys to local elements
    for (stk::mesh::EntityKey remoteElemKey : remoteElementKeys)
    {
        for (stk::mesh::EntityKey localNodeKey : remoteElementKeyToLocalNodeKeys[remoteElemKey])
        {
            stk::mesh::Entity node = bulkData.get_entity(localNodeKey);
            unsigned numElems = bulkData.num_elements(node);
            const stk::mesh::Entity* nodeElem = bulkData.begin_elements(node);
            for (unsigned elemIndex = 0; elemIndex < numElems; elemIndex++)
                remoteElementKeyToLocalNeighborElements[remoteElemKey].insert(nodeElem[elemIndex]);
        }
    }

    // reply with the local element keys of remote element keys to create local elements to remote element key map
    stk::CommSparse reply(bulkData.parallel());
    for (int phase = 0; phase < 2; phase++)
    {
        for (stk::mesh::EntityKey remoteElemKey : remoteElementKeys)
        {
            int procNum = remoteElementKeyToOwningProcessor[remoteElemKey];
            size_t numNeighbors = remoteElementKeyToLocalNeighborElements[remoteElemKey].size();
            reply.send_buffer(procNum).pack<stk::mesh::EntityKey>(remoteElemKey);
            reply.send_buffer(procNum).pack<size_t>(numNeighbors);
            for (stk::mesh::Entity localElem : remoteElementKeyToLocalNeighborElements[remoteElemKey])
            {
                stk::mesh::EntityKey localElemKey = bulkData.entity_key(localElem);
                reply.send_buffer(procNum).pack<stk::mesh::EntityKey>(localElemKey);
            }

        }
        if (0 == phase)
            reply.allocate_buffers();
        else
            reply.communicate();
    }

    //now get the local elements to store the neighboring remote elements mapping to their local elements
    for (int procRank = 0; procRank < numProcs; procRank++)
    {
        stk::CommBuffer& replyBuf = reply.recv_buffer(procRank);
        while (replyBuf.remaining())
        {
            size_t numNeighbors;
            stk::mesh::EntityKey localElemKey;
            stk::mesh::EntityKey remoteElemKey;
            replyBuf.unpack<stk::mesh::EntityKey>(localElemKey);
            replyBuf.unpack<size_t>(numNeighbors);
            stk::mesh::Entity localElem = bulkData.get_entity(localElemKey);
            for (unsigned neighborNum = 0; neighborNum < numNeighbors; neighborNum++)
            {
                replyBuf.unpack<stk::mesh::EntityKey>(remoteElemKey);
                localElementToRemoteElementKeys[localElem].insert(remoteElemKey);
            }
        }
    }

    // start to run game of life?
    unsigned lollipop;
    for (lollipop = 0; lollipop < numElems; lollipop++)
    {
        if (!(lollipop%5) || !(lollipop%7))
        {
            *stk::mesh::field_data(lifeField, elements[lollipop]) = 1;
            localActiveElements.push_back(elements[lollipop]);
        }
    }


    size_t fh = stkIo.create_output_mesh("FieldTestTrash.e", stk::io::WRITE_RESULTS);
    stkIo.add_field(fh, lifeField);
    stkIo.add_field(fh, neighborField);
    stkIo.write_output_mesh(fh);

    for (int time = 0; time <= 10; time++)
    {
        //determine which elements to check
        localElementsToVisit.clear();
        remoteElementKeysToVisit.clear();

        //get local and remote elements to visit
        for (stk::mesh::Entity localElem : localActiveElements)
        {
            localElementsToVisit.insert(localElem);
            for (stk::mesh::Entity localElemElem : localElementToLocalNeighborElements[localElem])
                localElementsToVisit.insert(localElemElem);

            for (stk::mesh::EntityKey remoteElemKey : localElementToRemoteElementKeys[localElem])
                remoteElementKeysToVisit.insert(remoteElemKey);
        }

        //communicate said local elements
        stk::CommSparse sneaker(bulkData.parallel());
        for (int phase = 0; phase < 2; phase++)
        {
            for (stk::mesh::EntityKey remoteElemKey : remoteElementKeysToVisit)
            {
                int procNum = remoteElementKeyToOwningProcessor[remoteElemKey];
                sneaker.send_buffer(procNum).pack<stk::mesh::EntityKey>(remoteElemKey);
            }

            if (0 == phase)
                sneaker.allocate_buffers();
            else
                sneaker.communicate();
        }
        for (int procNum = 0; procNum < numProcs; procNum++)
        {
            stk::CommBuffer& shoe = sneaker.recv_buffer(procNum);

            while (shoe.remaining())
            {
                stk::mesh::EntityKey localElemKey;
                shoe.unpack<stk::mesh::EntityKey>(localElemKey);
                localElementsToVisit.insert(bulkData.get_entity(localElemKey));
            }
        }

        //update neighbor values with local elements
        for (stk::mesh::Entity localElem : localElementsToVisit)
        {
            int* neighborVal = stk::mesh::field_data(neighborField, localElem);
            *neighborVal = 0;
            for (stk::mesh::Entity localElemElem : localElementToLocalNeighborElements[localElem])
                if (*stk::mesh::field_data(lifeField, localElemElem))
                    (*neighborVal)++;
        }

        //update neighbor values with remote elements
        stk::CommSparse sombrero(bulkData.parallel());
        for (int phase = 0 ; phase < 2; phase++)
        {
            for (stk::mesh::EntityKey remoteElemKey : remoteElementKeysToVisit)
            {
                int neighborVal = 0;
                int procNum = remoteElementKeyToOwningProcessor[remoteElemKey];

                for (stk::mesh::Entity localElem :
                        remoteElementKeyToLocalNeighborElements[remoteElemKey])
                    if (*stk::mesh::field_data(lifeField, localElem))
                        neighborVal++;

                sombrero.send_buffer(procNum).pack<stk::mesh::EntityKey>(remoteElemKey);
                sombrero.send_buffer(procNum).pack<int>(neighborVal);
            }
            if (0 == phase)
                sombrero.allocate_buffers();
            else if (1 == phase)
                sombrero.communicate();
        }

        for (int procNum = 0; procNum < numProcs; procNum++)
        {
            stk::CommBuffer& hat = sombrero.recv_buffer(procNum);
            while(hat.remaining())
            {
                stk::mesh::EntityKey localElemKey;
                int neighborVal;
                hat.unpack<stk::mesh::EntityKey>(localElemKey);
                hat.unpack<int>(neighborVal);
                *stk::mesh::field_data(neighborField, bulkData.get_entity(localElemKey))
                += neighborVal;
            }
        }

        //write output whenever
        stkIo.begin_output_step(fh, time);
        stkIo.write_defined_output_fields(fh);
        stkIo.end_output_step(fh);

        localActiveElements.clear();

        for (stk::mesh::Entity localElem : localElementsToVisit)
        {
            int* lifeVal = stk::mesh::field_data(lifeField, localElem);
            switch (*stk::mesh::field_data(neighborField, localElem))
            {
                case 4:
                    break;
                case 5:
                    *lifeVal = 1;
                    break;
                default:
                    *lifeVal = 0;
                    break;
            }
            if (*lifeVal)
                localActiveElements.push_back(localElem);
        }
    }
}





