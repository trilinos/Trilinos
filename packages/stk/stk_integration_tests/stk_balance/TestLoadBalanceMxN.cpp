#include <test_utils/OptionsForTesting.hpp>
#include <test_utils/MeshFixtureMxNRebalance.hpp>
#include <Ioss_DatabaseIO.h>
#include <Ioss_CommSet.h>
#include <Ioss_SideSet.h>
#include <Ioss_NodeBlock.h>
#include <Ioss_SideBlock.h>
#include <Ioss_Field.h>
#include <stk_io/IossBridge.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_balance/internal/MxNutils.hpp>
#include <stk_balance/internal/balanceMtoN.hpp>
#include <stk_balance/internal/StkBalanceUtils.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <stk_balance/internal/MtoNRebalancer.hpp>
namespace
{

class TestBalanceMxNRebalanceUsingInputFiles : public MeshFixtureMxNRebalance
{
protected:
    TestBalanceMxNRebalanceUsingInputFiles()
    : MeshFixtureMxNRebalance(), num_procs_initial_decomp(0),
    num_procs_target_decomp(0), mOutputFileName("subomain.exo"), mInputFileName("")
    {
    }

    virtual unsigned get_num_procs_initial_decomp() const { return num_procs_initial_decomp; }
    virtual unsigned get_num_procs_target_decomp()  const { return num_procs_target_decomp; }
    virtual std::string get_output_filename() const { return mOutputFileName; }
    virtual std::string get_input_mesh_file_name() const { return mInputFileName; }

    void set_options()
    {
        Options options = getOptionsForTest("none");
        test_options_are_set(options);
        init_member_data(options);
    }

    void init_member_data(const Options &options)
    {
        num_procs_initial_decomp = stk::parallel_machine_size(get_comm());
        num_procs_target_decomp = options.getNumTargetProcs();
        mInputFileName = options.getMeshFileName();
        mOutputFileName = options.getOutputFilename();
    }

    void setup_initial_mesh_from_last_time_step(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        create_target_decomp_field_on_entire_mesh();
        setup_mesh_from_last_time_step(get_input_mesh_file_name(), auraOption);
    }

    int numSteps = -1;
    double maxTime = 0.0;

private:
    void fill_time_data_from_last_time_step(stk::io::StkMeshIoBroker &stkIo)
    {
        numSteps = stkIo.get_num_time_steps();
        if(numSteps>0)
        {
            stkIo.read_defined_input_fields(numSteps);
            maxTime = stkIo.get_max_time();
        }
    }

    void setup_mesh_from_last_time_step(const std::string &meshSpecification, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        stk::io::StkMeshIoBroker stkIo;
        allocate_bulk(auraOption);
        stk::io::fill_mesh_preexisting(stkIo, meshSpecification, *bulkData);
        fill_time_data_from_last_time_step(stkIo);
    }

private:
    void test_options_are_set(const Options &options)
    {
        ASSERT_TRUE(options.getMeshFileName() != "none");
        ASSERT_TRUE(options.getNumTargetProcs() != 0);
    }

private:
    int num_procs_initial_decomp;
    int num_procs_target_decomp;
    std::string mOutputFileName;
    std::string mInputFileName;
};


std::vector<std::pair<stk::mesh::EntityId, int>> getSharingInfo(stk::mesh::BulkData& bulkData)
{
    stk::mesh::EntityVector sharedNodes;
    stk::mesh::get_selected_entities(bulkData.mesh_meta_data().globally_shared_part(),
                                     bulkData.buckets(stk::topology::NODE_RANK),
                                     sharedNodes);
    std::vector<std::pair<stk::mesh::EntityId, int>> nodeSharingInfo;
    nodeSharingInfo.reserve(8*sharedNodes.size());

    std::vector<int> sharingProcs;
    for(stk::mesh::Entity sharedNode : sharedNodes)
    {
        bulkData.comm_shared_procs(bulkData.entity_key(sharedNode), sharingProcs);
        for(unsigned j=0;j<sharingProcs.size();++j)
            nodeSharingInfo.push_back(std::make_pair(bulkData.identifier(sharedNode), sharingProcs[j]));
    }

    return nodeSharingInfo;
}

void verify_node_sharing_info(const std::vector<std::pair<stk::mesh::EntityId, int>> &nodeSharingInfo, const std::string& filename)
{
    stk::mesh::MetaData meta(3);
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, bulk);

    std::vector<std::pair<stk::mesh::EntityId, int>> nodeSharingInfoAfter = getSharingInfo(bulk);

    EXPECT_TRUE(nodeSharingInfo == nodeSharingInfoAfter);
}

TEST_F(TestBalanceMxNRebalanceUsingInputFiles, read4procswrite4procsFilesUsingIoss)
{
    if(stk::parallel_machine_size(get_comm())==4)
    {
        set_options();
        setup_initial_mesh_from_last_time_step(stk::mesh::BulkData::NO_AUTO_AURA);
        std::vector<std::pair<stk::mesh::EntityId, int>> nodeSharingInfo = getSharingInfo(get_bulk());
        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(get_bulk(), counts);
        int global_num_nodes = counts[stk::topology::NODE_RANK];
        int global_num_elems = counts[stk::topology::ELEM_RANK];

        stk::io::write_file_for_subdomain(get_output_filename(),
                                          get_bulk().parallel_rank(),
                                          get_bulk().parallel_size(),
                                          global_num_nodes,
                                          global_num_elems,
                                          get_bulk(),
                                          nodeSharingInfo,
                                          numSteps,
                                          maxTime);

        verify_node_sharing_info(nodeSharingInfo, get_output_filename());
    }
}

TEST_F(TestBalanceMxNRebalanceUsingInputFiles, MxN_decompositionWithoutAura)
{
    set_options();
    setup_initial_mesh_from_last_time_step(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::balance::internal::rebalanceMtoN(get_bulk(), *targetDecompField, get_num_procs_target_decomp(), get_output_filename(), numSteps, maxTime);
}

class BulkDataForBalance : public stk::mesh::BulkData
{
public:
    BulkDataForBalance(stk::mesh::MetaData& meta, MPI_Comm comm, stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::NO_AUTO_AURA) : BulkData(meta, comm, auraOption) {}
    void set_parallel(stk::Parallel input) { m_parallel = input; }
    void increment_sync_count() { m_meshModification.increment_sync_count(); }
};

class MxNRebalanceOnNProcs : public  stk::unit_test_util::MeshFixture
{
protected:
    void set_communicator(MPI_Comm comm) { communicator = comm; }

    virtual void allocate_bulk(stk::mesh::BulkData::AutomaticAuraOption auraOption,
                               unsigned bucketCapacity = stk::mesh::impl::BucketRepository::default_bucket_capacity)
    {
        if(nullptr == metaData)
            allocate_meta();

        ThrowRequireMsg(bucketCapacity == stk::mesh::impl::BucketRepository::default_bucket_capacity, "allocate_bulk: BulkDataForBalance doesn't use non-default bucket-capacity.");
        bulkData = new BulkDataForBalance(get_meta(), communicator, auraOption);
    }
};

bool thereAre16ElementsIn(stk::mesh::BulkData& bulkData)
{
    std::vector<size_t> mesh_counts;
    stk::mesh::comm_mesh_counts(bulkData, mesh_counts);
    return 16u == mesh_counts[stk::topology::ELEM_RANK];
}

void read_and_rebalance_mesh(stk::mesh::BulkData& bulk, const std::string& outputFilename)
{
    stk::balance::BasicZoltan2Settings graphSettings;
    stk::balance::balanceStkMesh(graphSettings, bulk);
    stk::io::write_mesh(outputFilename, bulk);
}

TEST_F(MxNRebalanceOnNProcs, testHexplateFrom4to8procs)
{
    std::string filename = stk::unit_test_util::get_option("-i", "hexplate.par");
    bool running_as_unit_test = filename == "hexplate.par";

    int numInput = stk::unit_test_util::get_command_line_option("-n", 4);
    ThrowRequireWithSierraHelpMsg(numInput>0);

    MPI_Comm globalComm = MPI_COMM_WORLD;

    int num_global_procs = stk::parallel_machine_size(globalComm);
    ThrowRequireWithSierraHelpMsg(num_global_procs>=numInput);

    int procId = stk::parallel_machine_rank(globalComm);

    int color = 0;
    if(procId>=numInput)
        color = 1;

    MPI_Comm localComm;
    MPI_Comm_split(globalComm, color, procId, &localComm);

    set_communicator(localComm);

    stk::mesh::BulkData* bulk_ptr = nullptr;
    stk::mesh::MetaData* meta_ptr = nullptr;

    if(color == 0)
    {
        setup_mesh(filename, stk::mesh::BulkData::NO_AUTO_AURA);

        if(running_as_unit_test)
        {
            EXPECT_TRUE(thereAre16ElementsIn(get_bulk()));
        }

        bulk_ptr = &get_bulk();
        meta_ptr = &get_meta();
    }
    else
    {
        int procFromSrc = 0;
        std::string tempFilename = stk::balance::internal::get_parallel_filename(procFromSrc, numInput, filename);
        meta_ptr = new stk::mesh::MetaData;
        bulk_ptr = new BulkDataForBalance(*meta_ptr, MPI_COMM_SELF);

        stk::io::StkMeshIoBroker stkIo;
        stkIo.set_bulk_data(*bulk_ptr);
        stkIo.add_mesh_database(tempFilename, stk::io::READ_MESH);
        stkIo.create_input_mesh();

        dynamic_cast<BulkDataForBalance*>(bulk_ptr)->initialize_face_adjacent_element_graph();
        dynamic_cast<BulkDataForBalance*>(bulk_ptr)->increment_sync_count();
    }

    std::string outputFilename = stk::unit_test_util::get_option("-o", "output.exo");

    BulkDataForBalance* bulkDataBalance = dynamic_cast<BulkDataForBalance*>(bulk_ptr);
    ThrowRequireWithSierraHelpMsg(bulkDataBalance!=nullptr);
    bulkDataBalance->set_parallel(stk::Parallel(globalComm));
    read_and_rebalance_mesh(*bulkDataBalance, outputFilename);

    if(running_as_unit_test)
    {
        EXPECT_TRUE(thereAre16ElementsIn(*bulkDataBalance));
    }

    bool doesProcHaveLocallyAllocatedMesh = color == 1;
    if(doesProcHaveLocallyAllocatedMesh)
    {
        delete bulk_ptr; bulk_ptr = nullptr;
        delete meta_ptr; meta_ptr = nullptr;
    }
}

}

