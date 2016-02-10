/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t, nullptr
#include <string>                       // for string
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include "../FaceCreationTestUtils.hpp"
#include <mpi.h>
#include <map>
#include <string>

namespace
{

class ElemElemGraphTester : public stk::mesh::ElemElemGraph
{
public:
    ElemElemGraphTester(stk::mesh::BulkData& bulkData, const stk::mesh::Selector &selector, const stk::mesh::Selector *air = nullptr) :
        stk::mesh::ElemElemGraph(bulkData, selector, air) { }

    stk::mesh::impl::ParallelGraphInfo& get_parallel_info() { return m_parallelInfoForGraphEdges.get_parallel_graph_info(); }
    stk::mesh::Entity get_entity(stk::mesh::impl::LocalId local_id) const { return m_local_id_to_element_entity[local_id]; }
    const stk::mesh::impl::SparseGraph& get_coincident_graph() const { return m_coincidentGraph; }
};

struct split_element_info
{
    stk::mesh::EntityId localElementId;
    stk::mesh::EntityId remoteElementId;
    int neighboringProc;
};

std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > get_split_coincident_elements(stk::mesh::BulkData& bulkData)
{
    stk::mesh::Selector sel = bulkData.mesh_meta_data().locally_owned_part();
    ElemElemGraphTester graph(bulkData, sel);
    const stk::mesh::impl::SparseGraph& coingraph = graph.get_coincident_graph();

    std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > badElements;

    for(const stk::mesh::impl::SparseGraph::value_type& extractedEdgesForElem : coingraph)
    {
        //const stk::mesh::impl::LocalId possibleMCE = extractedEdgesForElem.first;
        const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem = extractedEdgesForElem.second;
        for(const stk::mesh::GraphEdge& edge : coincidentEdgesForElem)
        {
            if(edge.elem2 < 0)
            {
                stk::mesh::Entity entity = graph.get_entity(edge.elem1);
                stk::mesh::EntityId id = bulkData.identifier(entity);
                stk::mesh::impl::ParallelGraphInfo& par_info = graph.get_parallel_info();
                stk::mesh::impl::ParallelGraphInfo::iterator iter = par_info.find(edge);
                ThrowRequireMsg(iter!=par_info.end(), "Program error. Contact sierra-help@sandia.gov for support.");
                badElements[id] = std::make_pair(-edge.elem2, iter->second.m_other_proc);
            }
        }
    }
    return badElements;
}

void write_mesh_diagnostics(const std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > & splitCoincidentElements, stk::ParallelMachine comm)
{
    int my_proc_id = stk::parallel_machine_rank(comm);
    std::ofstream out("mesh_diagnostics_failures_" + std::to_string(my_proc_id) + ".txt");
    for(const auto& item : splitCoincidentElements)
        out << "[" << my_proc_id << "] Element " << item.first << " is coincident with element " << item.second.first << " on processor " << item.second.second << std::endl;
    out.close();
}

struct TestCase
{
    std::string filename;
    int maxNumProcs;
    std::vector<split_element_info> expected_split_elements;
};

typedef std::vector<TestCase> TestCaseData;

const TestCaseData badDecomps =
{
  /* filename, max#procs, #side,  sideset */
    {"Aef.e",     2,        { {3u, 2u, 1}, {2u, 3u, 0} }},
    {"ef.e",      2,        { {1u, 2u, 1}, {2u, 1u, 0} }},
    {"AefB.e",    3,        { {}, {2u, 3u, 2}, {3u, 2u, 1} }},
};

class MeshChecker : public ::testing::Test
{
public:
    MeshChecker()
    {
        comm = MPI_COMM_WORLD;
    }

    void run_all_test_cases(const TestCaseData &testCases, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        for(const TestCase& testCase : testCases)
            if(stk::parallel_machine_size(get_comm()) == testCase.maxNumProcs)
                EXPECT_THROW(test_one_case(testCase, auraOption), std::logic_error);
    }

    void test_one_case(const TestCase &testCase,
                       stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        stk::mesh::MetaData metaData;
        stk::mesh::BulkData bulkData(metaData, get_comm(), auraOption);
        SideTestUtil::read_and_decompose_mesh(testCase.filename, bulkData);

        std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > splitCoincidentElements = get_split_coincident_elements(bulkData);

        for(const auto& item : splitCoincidentElements)
        {
            stk::mesh::EntityId localElementId = item.first;
            stk::mesh::EntityId remoteElementId = item.second.first;
            int neighboringProc = item.second.second;

            EXPECT_EQ(testCase.expected_split_elements[bulkData.parallel_rank()].localElementId, localElementId);
            EXPECT_EQ(testCase.expected_split_elements[bulkData.parallel_rank()].remoteElementId, remoteElementId);
            EXPECT_EQ(testCase.expected_split_elements[bulkData.parallel_rank()].neighboringProc, neighboringProc);
        }

        bool is_all_ok_locally = splitCoincidentElements.empty();
        bool is_all_ok_globally = stk::is_true_on_all_procs(bulkData.parallel(), is_all_ok_locally);
        if(!is_all_ok_locally)
            write_mesh_diagnostics(splitCoincidentElements, bulkData.parallel());
        ThrowRequireMsg(is_all_ok_globally, "Mesh diagnostics failed.");
    }

    MPI_Comm get_comm() const
    {
        return comm;
    }

private:
    MPI_Comm comm;
};


TEST_F(MeshChecker, diagnose_bad_meshes)
{
    run_all_test_cases(badDecomps, stk::mesh::BulkData::NO_AUTO_AURA);
}


}
