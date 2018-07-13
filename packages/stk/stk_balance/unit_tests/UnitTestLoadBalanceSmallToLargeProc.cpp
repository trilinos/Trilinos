#include <test_utils/MeshFixtureMxNRebalance.hpp>
#include <stk_balance/internal/MtoNRebalancer.hpp>

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

class Mesh1x1x4 : public MeshFixtureMxNRebalance
{
protected:

    void verify_node_sharing_info(const stk::mesh::EntityIdVector &goldSharedNodes, const stk::mesh::EntityVector &sharedNodes)
    {
        ASSERT_EQ(goldSharedNodes.size(), sharedNodes.size());
        for(size_t nodeIndex = 0; nodeIndex < sharedNodes.size(); nodeIndex++)
        {
            stk::mesh::EntityId nodeId = get_bulk().identifier(sharedNodes[nodeIndex]);
            EXPECT_EQ(goldSharedNodes[nodeIndex], nodeId);
        }
    }

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

}
