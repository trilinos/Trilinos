#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_util/environment/memory_util.hpp>
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
    std::vector<stk::mesh::EntityId> elemIds;
    std::vector< std::vector<stk::mesh::Entity> > elemNodes;
};

struct NodeConnectivityInfo
{
    std::map<stk::mesh::Entity, int> nodeConnMap;
    std::vector<int> nodeConn;
};

void get_ids_along_boundary(const stk::mesh::BulkData &bulkData, const std::string &meshSpec, DeletedMeshInfo &meshInfo)
{
    int nx, ny, nz;

    get_mesh_dimensions(meshSpec, nx, ny, nz);

    meshInfo.elemIds.resize(ny*nz);

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
                meshInfo.elemIds[count] = static_cast<stk::mesh::EntityId> (id);
                ++count;
            }
        }
    }

    meshInfo.elemIds.resize(count);
}

bool is_element_on_processor_boundary(const stk::mesh::ElemElemGraph& eeGraph, stk::mesh::Entity localElement)
{
    size_t numConn = eeGraph.get_num_connected_elems(localElement);

    bool onProcBoundary = false;
    for(size_t i=0; i<numConn; ++i)
    {
        if(!eeGraph.is_connected_elem_locally_owned(localElement, i))
            onProcBoundary = true;
    }
    return onProcBoundary;
}

bool will_deleting_element_create_orphaned_nodes(const stk::mesh::BulkData &bulkData,
                                                 stk::mesh::Entity localElement,
                                                 const NodeConnectivityInfo &nodeConnInfo)
{
    unsigned numNodes = bulkData.num_nodes(localElement);
    const stk::mesh::Entity * nodes = bulkData.begin(localElement, stk::topology::NODE_RANK);

    for(unsigned i=0; i<numNodes; ++i)
    {
        stk::mesh::Entity node = nodes[i];
        int index = -1;
        std::map<stk::mesh::Entity, int>::const_iterator it = nodeConnInfo.nodeConnMap.find(node);

        if(nodeConnInfo.nodeConnMap.end() != it)
            index = it->second;

        if((-1 != index) && (nodeConnInfo.nodeConn[index] <= 1))
            return true;
    }

    return false;
}

void decrement_node_connectivity(const stk::mesh::BulkData &bulkData,
                                 stk::mesh::Entity localElement,
                                 NodeConnectivityInfo &nodeConnInfo)
{
    unsigned numNodes = bulkData.num_nodes(localElement);
    const stk::mesh::Entity * nodes = bulkData.begin(localElement, stk::topology::NODE_RANK);

    for(unsigned i=0; i<numNodes; ++i)
    {
        stk::mesh::Entity node = nodes[i];
        int index = nodeConnInfo.nodeConnMap[node];
        nodeConnInfo.nodeConn[index]--;
    }
}

unsigned num_local_elements(const stk::mesh::BulkData &bulkData, stk::mesh::Entity node)
{
    unsigned numElems = bulkData.num_elements(node);
    const stk::mesh::Entity * elements = bulkData.begin(node, stk::topology::ELEM_RANK);

    for(unsigned i=0; i<numElems; ++i)
    {
        stk::mesh::Entity elem = elements[i];
        if(!bulkData.bucket(elem).owned())
            numElems--;
    }

    return numElems;
}

void initialize_node_connectivity_map(const stk::mesh::BulkData &bulkData,
                                      NodeConnectivityInfo &nodeConnInfo)
{
    stk::mesh::EntityVector allNodes;
    stk::mesh::get_selected_entities(bulkData.mesh_meta_data().universal_part(), bulkData.buckets(stk::topology::NODE_RANK), allNodes);

    nodeConnInfo.nodeConn.resize(allNodes.size());

    for(unsigned i=0; i<allNodes.size(); ++i)
    {
        stk::mesh::Entity node = allNodes[i];
        nodeConnInfo.nodeConnMap[node] = i;
        nodeConnInfo.nodeConn[i] = num_local_elements( bulkData, node );
    }
}

void populate_interior_perforated_mesh_ids(const stk::mesh::BulkData &bulkData,
                                          const stk::mesh::ElemElemGraph& eeGraph,
                                          const stk::mesh::EntityVector &localElems,
                                          NodeConnectivityInfo &nodeConnInfo,
                                          DeletedMeshInfo &meshInfo)
{
    for(unsigned lastIndex = 0; lastIndex < localElems.size(); ++lastIndex)
    {
        stk::mesh::Entity localElement = localElems[lastIndex];
        stk::mesh::EntityId elemId = bulkData.identifier( localElement );

        if(!is_element_on_processor_boundary(eeGraph, localElement) &&
           !will_deleting_element_create_orphaned_nodes(bulkData, localElement, nodeConnInfo))
        {
            meshInfo.elemIds.push_back(elemId);
            decrement_node_connectivity(bulkData, localElement, nodeConnInfo);
        }
    }
}

void populate_processor_boundary_perforated_mesh_ids(const stk::mesh::BulkData &bulkData,
                                                    const stk::mesh::ElemElemGraph& eeGraph,
                                                    const stk::mesh::EntityVector &localElems,
                                                    NodeConnectivityInfo &nodeConnInfo,
                                                    DeletedMeshInfo &meshInfo)
{
    for(unsigned lastIndex = 0; lastIndex < localElems.size(); ++lastIndex)
    {
        stk::mesh::Entity localElement = localElems[lastIndex];
        stk::mesh::EntityId elemId = bulkData.identifier( localElement );

        if(is_element_on_processor_boundary(eeGraph, localElement) &&
           !will_deleting_element_create_orphaned_nodes(bulkData, localElement, nodeConnInfo))
        {
            meshInfo.elemIds.push_back(elemId);
            decrement_node_connectivity(bulkData, localElement, nodeConnInfo);
        }
    }
}

void get_perforated_mesh_ids(const stk::mesh::BulkData &bulkData, const stk::mesh::ElemElemGraph& eeGraph, DeletedMeshInfo &meshInfo)
{

    stk::mesh::EntityVector localElems;
    stk::mesh::get_selected_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData.buckets(stk::topology::ELEM_RANK), localElems);

    meshInfo.elemIds.reserve(localElems.size());

    NodeConnectivityInfo nodeConnInfo;

    initialize_node_connectivity_map(bulkData, nodeConnInfo);

    populate_interior_perforated_mesh_ids(bulkData, eeGraph, localElems, nodeConnInfo, meshInfo);

    populate_processor_boundary_perforated_mesh_ids(bulkData, eeGraph, localElems, nodeConnInfo, meshInfo);

    //std::cout << "Deleting " << meshInfo.elemIds.size() << " perforated elements" << std::endl;
}

void print_memory(stk::ParallelMachine communicator, const std::string &preamble)
{
    double maxHwmInMB = stk::get_max_hwm_across_procs(communicator) / (1024.0 * 1024.0);

    size_t curr_max, curr_min, curr_avg;
    stk::get_current_memory_usage_across_processors(communicator, curr_max, curr_min, curr_avg);
    double max = curr_max / (1024.0 * 1024.0);
    curr_min /= (1024.0 * 1024.0);
    curr_avg /= (1024.0 * 1024.0);

    if(stk::parallel_machine_rank(communicator) == 0)
    {
        std::ostringstream os;
        os << preamble << " high water mark = " << maxHwmInMB << " Mb, current max = " << max << " Mb" << std::endl;
        std::cerr << os.str();
    }
}

class ElemElemGraphUpdaterWithTiming : public stk::mesh::ElemElemGraphUpdater
{
public:
    ElemElemGraphUpdaterWithTiming(stk::mesh::BulkData &bulk, stk::mesh::ElemElemGraph &elemGraph_, stk::diag::Timer &addElemTimer_, stk::diag::Timer &deleteElemTimer_)
    : ElemElemGraphUpdater(bulk, elemGraph_), addElemTimer(addElemTimer_), deleteElemTimer(deleteElemTimer_), bulkData(bulk)

    {
    }

    virtual void finished_modification_end_notification()
    {
        stk::diag::TimeBlockSynchronized timerStartSynchronizedAcrossProcessors(addElemTimer, bulkData.parallel());
        stk::mesh::ElemElemGraphUpdater::finished_modification_end_notification();
    }

    virtual void started_modification_end_notification()
    {
        stk::diag::TimeBlockSynchronized timerStartSynchronizedAcrossProcessors(deleteElemTimer, bulkData.parallel());
        stk::mesh::ElemElemGraphUpdater::started_modification_end_notification();
    }

private:
    stk::diag::Timer &addElemTimer;
    stk::diag::Timer &deleteElemTimer;
    stk::mesh::BulkData &bulkData;
};

class ElementGraphPerformance : public stk::unit_test_util::PerformanceTester
{
public:
    ElementGraphPerformance(stk::mesh::BulkData &bulk, const std::string &fileSpec) :
            stk::unit_test_util::PerformanceTester(bulk.parallel()),
            bulkData(bulk),
            meshSpec(fileSpec),
            CHILDMASK2(2),
            addElemTimer("add elements", CHILDMASK2, rootTimer),
            deleteElemTimer("delete elements", CHILDMASK2, rootTimer),
            createGraphTimer("create graph", CHILDMASK2, rootTimer)
    {
        double startTime = stk::wall_time();

        enabledTimerSet.setEnabledTimerMask(CHILDMASK2);
        {
            print_memory(bulk.parallel(), "Before elem graph creation");
            stk::diag::TimeBlockSynchronized timerStartSynchronizedAcrossProcessors(createGraphTimer, communicator);
            elemGraph = new stk::mesh::ElemElemGraph(bulk);
            print_memory(bulk.parallel(), "After elem graph creation");
        }

        elemGraphUpdater = new ElemElemGraphUpdaterWithTiming(bulk, *elemGraph, addElemTimer, deleteElemTimer);
        bulk.register_observer(elemGraphUpdater);

        duration += stk::wall_time() - startTime;
    }

    virtual ~ElementGraphPerformance()
    {
        delete elemGraph;
        delete elemGraphUpdater;
    }

protected:
    void delete_elements(DeletedMeshInfo &meshInfo)
    {
        meshInfo.elemNodes.resize(meshInfo.elemIds.size());

        bulkData.modification_begin();
        for(size_t i = 0; i < meshInfo.elemIds.size(); ++i)
        {
            stk::mesh::EntityId elemId = meshInfo.elemIds[i];
            stk::mesh::Entity elem = bulkData.get_entity(stk::topology::ELEM_RANK, elemId);

            EXPECT_TRUE(bulkData.is_valid(elem));

            meshInfo.elemNodes[i].assign(bulkData.begin_nodes(elem), bulkData.end_nodes(elem));

            bulkData.destroy_entity(elem);
        }
        bulkData.modification_end();
    }

    void add_elements(DeletedMeshInfo &meshInfo)
    {
        bulkData.modification_begin();
        stk::mesh::Part &hexPart = bulkData.mesh_meta_data().get_topology_root_part(stk::topology::HEX_8);
        stk::mesh::Part *block_1 = bulkData.mesh_meta_data().get_part("block_1");
        stk::mesh::PartVector parts{&hexPart, block_1};
        for(size_t i = 0; i < meshInfo.elemIds.size(); ++i)
        {
            stk::mesh::EntityId elemId = meshInfo.elemIds[i];
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

        delete_elements(meshInfo);
        stk::unit_test_util::write_mesh_using_stk_io("afterDelete.g", bulkData, bulkData.parallel());
        unsigned newGraphEdgeCount = elemGraph->num_edges();
        ASSERT_TRUE(oldGraphEdgeCount > newGraphEdgeCount);

        add_elements(meshInfo);

        unsigned finalGraphEdgeCount = elemGraph->num_edges();
        ASSERT_EQ(oldGraphEdgeCount, finalGraphEdgeCount);
    }
    virtual size_t get_value_to_output_as_iteration_count()
    {
        return 1;
    }

    stk::mesh::BulkData &bulkData;
    stk::mesh::ElemElemGraph *elemGraph = nullptr;
    ElemElemGraphUpdaterWithTiming *elemGraphUpdater = nullptr;
    const std::string meshSpec;

    const int CHILDMASK2;
    stk::diag::Timer addElemTimer;
    stk::diag::Timer deleteElemTimer;
    stk::diag::Timer createGraphTimer;
};

class PerforatedElementGraphPerformance : public ElementGraphPerformance
{
public:
    PerforatedElementGraphPerformance(stk::mesh::BulkData &bulk, const std::string &fileSpec) :
        ElementGraphPerformance(bulk, fileSpec)
    {

    }

    ~PerforatedElementGraphPerformance()
    {

    }

protected:
    virtual void run_algorithm_to_time()
    {
        DeletedMeshInfo meshInfo;

        unsigned oldGraphEdgeCount = elemGraph->num_edges();
        get_perforated_mesh_ids(bulkData, *elemGraph, meshInfo);

        delete_elements(meshInfo);
        stk::unit_test_util::write_mesh_using_stk_io("afterDelete.g", bulkData, bulkData.parallel());
        unsigned newGraphEdgeCount = elemGraph->num_edges();
        ASSERT_TRUE(oldGraphEdgeCount > newGraphEdgeCount);

        add_elements(meshInfo);
        stk::unit_test_util::write_mesh_using_stk_io("afterAdd.g", bulkData, bulkData.parallel());
        unsigned finalGraphEdgeCount = elemGraph->num_edges();
        ASSERT_EQ(oldGraphEdgeCount, finalGraphEdgeCount);
    }
};

class PerforatedElementGraphPerformanceTest : public stk::unit_test_util::MeshFixture
{
protected:
    void run_element_graph_perf_test()
    {
        PerforatedElementGraphPerformance perfTester(get_bulk(), get_mesh_spec());
        perfTester.run_performance_test();
    }
    std::string get_mesh_spec()
    {
        return unitTestUtils::getOption("-file", "NO_FILE_SPECIFIED");
    }
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

TEST_F(ElementGraphPerformanceTest, read_mesh_with_auto_decomp_no_aura)
{
    allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::unit_test_util::read_from_serial_file_and_decompose(get_mesh_spec(), get_bulk(), "rcb");

    run_element_graph_perf_test();
}

TEST_F(ElementGraphPerformanceTest, read_mesh)
{
    print_memory(MPI_COMM_WORLD, "Before mesh creation");
    setup_mesh(get_mesh_spec(), stk::mesh::BulkData::AUTO_AURA);

    run_element_graph_perf_test();
}

TEST_F(ElementGraphPerformanceTest, read_mesh_no_aura)
{
    print_memory(MPI_COMM_WORLD, "Before mesh creation");
    setup_mesh(get_mesh_spec(), stk::mesh::BulkData::NO_AUTO_AURA);

    run_element_graph_perf_test();
}

TEST_F(PerforatedElementGraphPerformanceTest, read_mesh_with_auto_decomp)
{
    allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
    stk::unit_test_util::read_from_serial_file_and_decompose(get_mesh_spec(), get_bulk(), "rcb");

    run_element_graph_perf_test();
}

TEST_F(PerforatedElementGraphPerformanceTest, read_mesh_with_auto_decomp_no_aura)
{
    allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::unit_test_util::read_from_serial_file_and_decompose(get_mesh_spec(), get_bulk(), "rcb");

    run_element_graph_perf_test();
}

}

