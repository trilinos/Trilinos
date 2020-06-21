#include <stk_balance/internal/SubdomainFileWriter.hpp>
#include <test_utils/MeshFixtureMxNRebalance.hpp>
#include <stk_balance/internal/MtoNRebalancer.hpp>

#include "stk_unit_test_utils/getOption.h"
#include "stk_mesh/base/Comm.hpp"
#include "stk_util/parallel/ParallelVectorConcat.hpp"
#include "stk_mesh/base/DestroyElements.hpp"

namespace
{

class TestBalanceBalanceSmallToLarge : public MeshFixtureMxNRebalance
{
protected:
    TestBalanceBalanceSmallToLarge() : MeshFixtureMxNRebalance() {}

    virtual unsigned get_x() const { return 3; }
    virtual unsigned get_y() const { return 3; }
    virtual unsigned get_z() const { return 3; }

    virtual unsigned get_num_procs_initial_decomp() const { return 2; }
    virtual unsigned get_num_procs_target_decomp()  const { return 3; }
};

TEST_F(TestBalanceBalanceSmallToLarge, MxN_decompositionWithAura)
{
    if(stk::parallel_machine_size(get_comm()) == static_cast<int>(get_num_procs_initial_decomp()))
        setup_and_test_balance_of_mesh(stk::mesh::BulkData::AUTO_AURA);
}

TEST_F(TestBalanceBalanceSmallToLarge, MxN_decompositionWithoutAura)
{
    if(stk::parallel_machine_size(get_comm()) == static_cast<int>(get_num_procs_initial_decomp()))
        setup_and_test_balance_of_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
}

class TestBalanceMtoM : public MeshFixtureMxNRebalance
{
protected:
    TestBalanceMtoM() : MeshFixtureMxNRebalance() {}

    virtual unsigned get_x() const { return 3; }
    virtual unsigned get_y() const { return 3; }
    virtual unsigned get_z() const { return 3; }

    virtual unsigned get_num_procs_initial_decomp() const { return 2; }
    virtual unsigned get_num_procs_target_decomp()  const { return 2; }
};

TEST_F(TestBalanceMtoM, DISABLED_MxM_decompositionWithoutAura)
{
    if(stk::parallel_machine_size(get_comm()) == static_cast<int>(get_num_procs_initial_decomp())) {
        setup_initial_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        EXPECT_FALSE(stk::balance::internal::rebalanceMtoN(get_bulk(), *targetDecompField, get_num_procs_target_decomp(), get_output_filename()));
    }
}

class Mesh1x1x4 : public MeshFixtureMxNRebalance
{
protected:
    virtual unsigned get_num_procs_initial_decomp() const { return 2; }
    virtual unsigned get_num_procs_target_decomp()  const { return 4; }
    virtual std::string get_output_filename() const { return "junk.g"; }
    virtual std::string get_input_mesh_file_name() const { return "generated:1x1x4"; }
};

TEST_F(Mesh1x1x4, read2procsWrite4procsFilesUsingGeneratedMesh)
{
    if(stk::parallel_machine_size(get_comm()) == 2)
    {
        std::vector<stk::io::EntitySharingInfo> goldSharedNodesPerSubdomain =
        {
            {{ 5, 1}, { 6, 1}, { 7, 1}, { 8, 1}},
            {{ 5, 0}, { 6, 0}, { 7, 0}, { 8, 0}, { 9, 2}, {10, 2}, {11, 2}, {12, 2}},
            {{ 9, 1}, {10, 1}, {11, 1}, {12, 1}, {13, 3}, {14, 3}, {15, 3}, {16, 3}},
            {{13, 2}, {14, 2}, {15, 2}, {16, 2}}
        };

        std::vector<unsigned> targetProc_to_startingProc = stk::balance::internal::assign_target_subdomains_roundrobin_to_procs(stk::parallel_machine_size(get_comm()), get_num_procs_target_decomp());
        setup_initial_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

        stk::balance::BasicZoltan2Settings graphSettings;
        stk::balance::internal::MtoNRebalancer rebalancer(get_bulk(), *targetDecompField, graphSettings, get_num_procs_target_decomp());

        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(get_bulk(), counts);
        int global_num_nodes = counts[stk::topology::NODE_RANK];
        int global_num_elems = counts[stk::topology::ELEM_RANK];

        rebalancer.move_subdomains_such_that_entire_subdomain_doesnt_span_proc_boundaries(targetProc_to_startingProc);

        for(unsigned subdomain = 0; subdomain < targetProc_to_startingProc.size(); subdomain++)
        {
            if(rebalancer.does_this_proc_own_subdomain(targetProc_to_startingProc[subdomain]))
            {
                stk::io::EntitySharingInfo nodeSharingInfo = rebalancer.get_node_sharing_info(subdomain);
                EXPECT_TRUE(goldSharedNodesPerSubdomain[subdomain] == nodeSharingInfo);
                rebalancer.create_subdomain_and_write("testing.g", subdomain, global_num_nodes, global_num_elems, nodeSharingInfo);
            }
        }
    }
}

void expect_and_unlink_file(const std::string& baseName, int numProc, int procId)
{
    std::string file0Name = baseName + "." + std::to_string(numProc) + "." + std::to_string(procId);
    std::ifstream file0(file0Name);
    EXPECT_TRUE(!file0.fail());
    unlink(file0Name.c_str());
}

TEST(SomeProcessorsWithNoElements, writeSubdomains_onlyProcsWithElementsWrite)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    const std::string inputMesh = "generated:1x1x4";
    const std::string outputMesh = "reduced.g";

    if(stk::parallel_machine_size(comm) == 4)
    {
        stk::mesh::MetaData meta;
        stk::mesh::Field<double> &targetDecompField = meta.declare_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK, "junk", 1);
        stk::mesh::put_field_on_mesh(targetDecompField, meta.universal_part(), static_cast<double*>(nullptr));
        stk::mesh::BulkData bulk(meta, comm);
        stk::io::fill_mesh(inputMesh, bulk);

        stk::mesh::EntityVector elementsToDestroy;
        if(stk::parallel_machine_rank(comm)==0 || stk::parallel_machine_rank(comm)==3)
        {
            stk::mesh::get_selected_entities(meta.locally_owned_part(), bulk.buckets(stk::topology::ELEM_RANK), elementsToDestroy);
        }
        stk::mesh::destroy_elements(bulk, elementsToDestroy);

        int includeMe = 0;
        int numTarget = 0;
        std::tie(includeMe, numTarget) = stk::balance::internal::get_included_and_num_target_procs(bulk, comm);

        EXPECT_EQ(2, numTarget);

        int mySubdomain = stk::balance::internal::get_subdomain_index(includeMe, comm);

        if(stk::parallel_machine_rank(comm)==2)
        {
            EXPECT_EQ(1, includeMe);
            EXPECT_EQ(1, mySubdomain);
        }
        else if(stk::parallel_machine_rank(comm)==1)
        {
            EXPECT_EQ(1, includeMe);
            EXPECT_EQ(0, mySubdomain);
        }
        else
        {
            EXPECT_EQ(0, includeMe);
            EXPECT_EQ(-1, mySubdomain);
        }

        stk::balance::internal::write_subdomain_files(bulk, numTarget, mySubdomain, outputMesh);

        stk::parallel_machine_barrier(comm);
        if(stk::parallel_machine_rank(comm)==0)
        {
            expect_and_unlink_file(outputMesh, numTarget, 0);
            expect_and_unlink_file(outputMesh, numTarget, 1);
        }
    }
}

}
