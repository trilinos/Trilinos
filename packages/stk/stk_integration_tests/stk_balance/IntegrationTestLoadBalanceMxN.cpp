// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
//
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <test_utils/OptionsForTesting.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_balance/internal/StkBalanceUtils.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/balance.hpp>

namespace
{

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

