#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphUpdater.hpp>
#include <tokenize.h>                   // for tokenize
#include <cmath>                        // for atan2, cos, sin
#include <cstdlib>                      // for strtod, nullptr, strtol, exit, etc
#include <cstring>                      // for memcpy
#include <iomanip>                      // for operator<<, setw
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <string>                       // for string, operator==, etc
#include <vector>
#include <Ioss_EntityType.h>            // for EntityType, etc

namespace
{

std::string get_file_name(const std::string &parameters)
{
    std::string filename = parameters;
    size_t colon = filename.find(':');
    if (colon != std::string::npos && colon > 0) {
      filename = filename.substr(colon+1);
    }

    return filename;
}

void get_mesh_dimensions(const std::string &parameters, int &numX, int &numY, int &numZ)
{
    std::string filename = get_file_name(parameters);

    auto params = Ioss::tokenize(filename, "/");
    auto groups = Ioss::tokenize(params[params.size()-1], "|+");

    // First 'group' is the interval specification -- IxJxK
    auto tokens = Ioss::tokenize(groups[0], "x");
    assert(tokens.size() == 3);

    numX = std::strtol(tokens[0].c_str(), nullptr, 10);
    numY = std::strtol(tokens[1].c_str(), nullptr, 10);
    numZ = std::strtol(tokens[2].c_str(), nullptr, 10);
}

struct DeletedMeshInfo
{
    std::vector<stk::mesh::EntityId> boundaryIds;
    std::vector< std::vector<stk::mesh::Entity> > elemNodes;
};

void get_ids_along_boundary(const stk::mesh::BulkData &bulkData, const std::string &meshSpec, DeletedMeshInfo &meshInfo)
{
    int nx, ny, nz;

    get_mesh_dimensions(meshSpec, nx, ny, nz);

    meshInfo.boundaryIds.resize(ny*nz);

    int count = 0;

    for(int z = 0; z < nz; ++z)
    {
        for(int y = 0; y < ny; ++y)
        {
            if((y%2 == 0) || ((ny-1) == y)) continue;

            int id = 1 + nx*(y + ny*z);

            stk::mesh::EntityId elemId = static_cast<stk::mesh::EntityId> (id);
            stk::mesh::Entity elem = bulkData.get_entity(stk::topology::ELEM_RANK, elemId);

            if(bulkData.is_valid(elem) && bulkData.bucket(elem).owned())
            {
                meshInfo.boundaryIds[count] = static_cast<stk::mesh::EntityId> (id);
                ++count;
            }
        }
    }

    meshInfo.boundaryIds.resize(count);
}


class ElementGraphPerformance : public stk::unit_test_util::PerformanceTester
{
public:
    ElementGraphPerformance(stk::mesh::BulkData &bulk, const std::string &fileSpec) :
            stk::unit_test_util::PerformanceTester(bulk.parallel()),
            bulkData(bulk),
            meshSpec(fileSpec),
            addElemTimer("add elements", CHILDMASK1, childTimer),
            deleteElemTimer("delete elements", CHILDMASK1, childTimer),
            createGraphTimer("create graph", CHILDMASK1, childTimer)
    {
        {
            stk::diag::TimeBlockSynchronized timerStartSynchronizedAcrossProcessors(createGraphTimer, communicator);
            elemGraph = new stk::mesh::ElemElemGraph(bulk, bulk.mesh_meta_data().universal_part());
        }

        elemGraphUpdater = new stk::mesh::ElemElemGraphUpdater(bulk, *elemGraph);
        bulk.register_observer(elemGraphUpdater);
    }

    ~ElementGraphPerformance()
    {
        delete elemGraph;
        delete elemGraphUpdater;
    }

protected:
    void delete_boundary_elements(DeletedMeshInfo &meshInfo)
    {
        stk::diag::TimeBlockSynchronized timerStartSynchronizedAcrossProcessors(deleteElemTimer, communicator);
        meshInfo.elemNodes.resize(meshInfo.boundaryIds.size());

        bulkData.modification_begin();
        for(size_t i = 0; i < meshInfo.boundaryIds.size(); ++i)
        {
            stk::mesh::EntityId elemId = meshInfo.boundaryIds[i];
            stk::mesh::Entity elem = bulkData.get_entity(stk::topology::ELEM_RANK, elemId);

            EXPECT_TRUE(bulkData.is_valid(elem));

            meshInfo.elemNodes[i].assign(bulkData.begin_nodes(elem), bulkData.end_nodes(elem));

            bulkData.destroy_entity(elem);
        }
        bulkData.modification_end();
    }

    void add_boundary_elements(DeletedMeshInfo &meshInfo)
    {
        stk::diag::TimeBlockSynchronized timerStartSynchronizedAcrossProcessors(addElemTimer, communicator);

        bulkData.modification_begin();
        stk::mesh::Part &hexPart = bulkData.mesh_meta_data().get_topology_root_part(stk::topology::HEX_8);
        stk::mesh::Part *block_1 = bulkData.mesh_meta_data().get_part("block_1");
        stk::mesh::PartVector parts{&hexPart, block_1};
        for(size_t i = 0; i < meshInfo.boundaryIds.size(); ++i)
        {
            stk::mesh::EntityId elemId = meshInfo.boundaryIds[i];
            stk::mesh::Entity elem = bulkData.declare_entity(stk::topology::ELEM_RANK, elemId, parts);

            for(size_t j = 0; j<meshInfo.elemNodes[i].size(); ++j)
            {
                stk::mesh::Entity node = meshInfo.elemNodes[i][j];
                bulkData.declare_relation(elem, node, j);
            }
        }
        bulkData.modification_end();
    }

    virtual void run_algorithm_to_time()
    {
        DeletedMeshInfo meshInfo;

        unsigned oldGraphEdgeCount = elemGraph->num_edges();
        get_ids_along_boundary(bulkData, meshSpec, meshInfo);

        delete_boundary_elements(meshInfo);

        unsigned newGraphEdgeCount = elemGraph->num_edges();
        ASSERT_TRUE(oldGraphEdgeCount > newGraphEdgeCount);

        add_boundary_elements(meshInfo);

        unsigned finalGraphEdgeCount = elemGraph->num_edges();
        ASSERT_EQ(oldGraphEdgeCount, finalGraphEdgeCount);
    }
    virtual size_t get_value_to_output_as_iteration_count()
    {
        return 0;
    }

    stk::mesh::BulkData &bulkData;
    stk::mesh::ElemElemGraph *elemGraph = nullptr;
    stk::mesh::ElemElemGraphUpdater *elemGraphUpdater = nullptr;
    const std::string meshSpec;

    stk::diag::Timer addElemTimer;
    stk::diag::Timer deleteElemTimer;
    stk::diag::Timer createGraphTimer;
};

class ElementGraphPerformanceTest : public stk::unit_test_util::MeshFixture
{
protected:
    void run_element_graph_perf_test()
    {
        ElementGraphPerformance perfTester(get_bulk(), get_mesh_spec());
        perfTester.run_performance_test();
    }
    std::string get_mesh_spec()
    {
        return unitTestUtils::getOption("-file", "NO_FILE_SPECIFIED");
    }
};

TEST_F(ElementGraphPerformanceTest, read_mesh_with_auto_decomp)
{
    allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
    stk::unit_test_util::read_from_serial_file_and_decompose(get_mesh_spec(), get_bulk(), "rcb");

    run_element_graph_perf_test();
}

TEST_F(ElementGraphPerformanceTest, read_mesh)
{
    setup_mesh(get_mesh_spec(), stk::mesh::BulkData::AUTO_AURA);

    run_element_graph_perf_test();
}

}

