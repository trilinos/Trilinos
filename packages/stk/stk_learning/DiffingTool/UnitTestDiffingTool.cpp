#include <gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <stdlib.h>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/CreateFaces.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/environment/EnvData.hpp>

#include <stk_util/parallel/DebugTool.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

#include "stk_unit_test_utils/unittestMeshUtils.hpp"

namespace 
{

void setup_env_data(stk::EnvData& env_data, stk::ParallelMachine global_comm, int color)
{
    env_data.m_parallelComm = splitComm(color, global_comm);
    env_data.m_worldComm = global_comm;
}

void pack_string(stk::CommBuffer& buf, const std::string& name)
{
    buf.pack<size_t>(name.size());
    buf.pack<char>(name.data(), name.size());
}

std::string unpack_string(stk::CommBuffer& buf)
{
    size_t num_chars = 0;
    buf.unpack<size_t>(num_chars);
    std::string name(num_chars, ' ');
    buf.unpack<char>(&name[0], num_chars);
    return name;
}

void pack_part_names(stk::CommBuffer& buf, const stk::mesh::PartVector& parts)
{
    for(size_t i=0; i<parts.size(); ++i)
    {
        pack_string(buf, parts[i]->name());
    }
}

int parallel_sum(stk::ParallelMachine comm, int numLocal)
{
    int numGlobal = 0;
    stk::all_reduce_sum(comm, &numLocal, &numGlobal, 1);
    return numGlobal;
}

int get_global_part_differences(stk::ParallelMachine comm, int numLocalDiffs)
{
    return parallel_sum(comm, numLocalDiffs);
}

void communicate_part_differences_to_apps(stk::ParallelMachine comm, int numLocalDiffs)
{
    get_global_part_differences(comm, numLocalDiffs);
}

int get_global_bucket_part_membership_differences(stk::ParallelMachine comm, int numLocalDiffs)
{
    return parallel_sum(comm, numLocalDiffs);
}

int get_global_part_differences_for_app(stk::ParallelMachine comm)
{
    int numLocalDiffs = 0;
    return get_global_part_differences(comm, numLocalDiffs);
}

int get_global_bucket_count_differences(stk::ParallelMachine comm, int numLocalDiffs)
{
    return parallel_sum(comm, numLocalDiffs);
}

void allocate_or_communicate(int iphase, stk::CommSparse& comm)
{
    if (iphase == 0)
    {
        comm.allocate_buffers();
    }
    else
    {
        comm.communicate();
    }
}

void wait_for_apps_to_pack_and_send_data(stk::CommSparse& comm)
{
    comm.allocate_buffers();
    comm.communicate();
}

stk::CommBuffer& get_comm_buffer_for_destination_proc(stk::CommSparse& comm)
{
    return comm.send_buffer(getDestinationProc(comm.parallel()));
}

void send_part_names_to_diffing_tool(const stk::mesh::BulkData& bulk, stk::ParallelMachine communicator)
{
    stk::CommSparse comm(communicator);
    for(int iphase = 0; iphase < 2; ++iphase)
    {
        pack_part_names(get_comm_buffer_for_destination_proc(comm), bulk.mesh_meta_data().get_parts());
        allocate_or_communicate(iphase, comm);
    }
}

bool parts_match(const stk::mesh::BulkData& bulk, stk::EnvData& env_data)
{
    send_part_names_to_diffing_tool(bulk, env_data.m_worldComm);
    return get_global_part_differences_for_app(env_data.m_worldComm) == 0;
}

std::vector<std::string> unpack_buffer_to_string_vector(stk::CommBuffer& buff)
{
    std::vector<std::string> partNames;
    while(buff.remaining())
    {
        partNames.push_back(unpack_string(buff));
    }
    return partNames;
}

std::vector<std::string> get_sorted_part_names(stk::CommBuffer& buff)
{
    std::vector<std::string> partNames = unpack_buffer_to_string_vector(buff);
    std::sort(partNames.begin(), partNames.end());
    return partNames;
}

std::string join_string_vector(const std::vector<std::string>& strings, const std::string& sep)
{
    std::ostringstream os;
    for(size_t i=0;i<strings.size();++i)
    {
        os << strings[i] << sep;
    }
    return os.str();
}

void create_output_message(std::ostringstream &os, std::vector<std::string>& part_diffs)
{
    if(!part_diffs.empty())
    {
        os << join_string_vector(part_diffs, "\t") << std::endl;
    }
}

std::vector<std::string> part_name_set_difference(const std::vector<std::string> & partNames1, const std::vector<std::string> & partNames2)
{
    std::vector<std::string> part_diffs(partNames1.size());
    auto iter = std::set_difference(partNames1.begin(), partNames1.end(), partNames2.begin(), partNames2.end(), part_diffs.begin());
    part_diffs.resize(iter - part_diffs.begin());
    return part_diffs;
}

int part_name_set_difference_message(const std::vector<std::string> & partNames1, const std::vector<std::string> & partNames2, std::ostringstream& os)
{
    std::vector<std::string> part_diffs = part_name_set_difference(partNames1, partNames2);
    create_output_message(os, part_diffs);
    return part_diffs.size();
}

int get_parts_in_first_not_in_second(const std::vector<std::string> & partNames1, const std::vector<std::string> & partNames2, std::ostringstream & globalss)
{
    std::ostringstream os("Parts in first app, not in second app: ");
    int num_diffs = part_name_set_difference_message(partNames1, partNames2, os);
    if (num_diffs > 0)
    {
        globalss << os.str();
    }
    return num_diffs;
}

int get_parts_in_second_not_in_first(const std::vector<std::string> & partNames1, const std::vector<std::string> & partNames2, std::ostringstream & globalss)
{
    std::ostringstream os("Parts in second app, not in first app: ");
    int num_diffs = part_name_set_difference_message(partNames2, partNames1, os);
    if (num_diffs > 0)
    {
        globalss << os.str();
    }
    return num_diffs;
}

int compare_part_names(const std::vector<std::string> & partNames1, const std::vector<std::string> & partNames2, std::ostringstream & globalss)
{
    return get_parts_in_first_not_in_second(partNames1, partNames2, globalss) +
           get_parts_in_second_not_in_first(partNames1, partNames2, globalss);
}

void unpack_part_names(stk::EnvData & env_data, std::vector<std::string> & partNames1, std::vector<std::string> & partNames2)
{
    stk::CommSparse comm(env_data.m_worldComm);
    wait_for_apps_to_pack_and_send_data(comm);
    std::pair<int,int> sourceProcs = getSourceProcs(env_data.m_worldComm);
    partNames1 = get_sorted_part_names(comm.recv_buffer(sourceProcs.first));
    partNames2 = get_sorted_part_names(comm.recv_buffer(sourceProcs.second));
}

void compare_parts(stk::EnvData & env_data, std::ostringstream& globalss)
{
    std::vector<std::string> partNames1, partNames2;
    unpack_part_names(env_data, partNames1, partNames2);
    int num_diffs = compare_part_names(partNames1, partNames2, globalss);
    communicate_part_differences_to_apps(env_data.m_worldComm, num_diffs);
}

std::pair<int,int> getXYdim()
{
    int xdim = 2;
    int ydim = 1;
    return std::make_pair(xdim, ydim);
}

std::string get_generated_mesh_spec(int numProcs)
{
    std::ostringstream os;
    std::pair<int,int> xydim = getXYdim();
    os << "generated:" << xydim.first << "x" << xydim.second << "x" << numProcs;
    return os.str();
}

class TestMesh
{
public:
    TestMesh(stk::ParallelMachine comm)
    :m_meta_data(), m_bulk_data(m_meta_data, comm)
    {
        stk::unit_test_util::fill_mesh_using_stk_io(get_generated_mesh_spec(m_bulk_data.parallel_size()), m_bulk_data, comm);
        stk::mesh::create_edges(m_bulk_data);
        stk::mesh::create_faces(m_bulk_data);
    }

    stk::mesh::BulkData & get_bulk_data()
    {
        return m_bulk_data;
    }

private:
    stk::mesh::MetaData m_meta_data;
    stk::mesh::BulkData m_bulk_data;
};

void expect_mpmd(stk::EnvData & env_data)
{
    int num_global_procs = stk::parallel_machine_size(env_data.m_worldComm);
    int num_local_procs = stk::parallel_machine_size(env_data.m_parallelComm);
    EXPECT_TRUE(num_global_procs%3 == 0 && num_local_procs==num_global_procs/3);
}

void expect_num_elements(const stk::mesh::BulkData & bulk)
{
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(bulk, counts);
    unsigned zdim = bulk.parallel_size();
    std::pair<int,int> xydim = getXYdim();
    size_t expected_num_elems = xydim.first*xydim.second*zdim;
    EXPECT_EQ(expected_num_elems, counts[stk::topology::ELEM_RANK]);
}

int bucket_counts_match(const stk::mesh::BulkData& bulk, stk::EnvData& env_data)
{
    stk::CommSparse comm(env_data.m_worldComm);
    int destinationProc = getDestinationProc(env_data.m_worldComm);
    stk::CommBuffer& buf = comm.send_buffer(destinationProc);
    for(int iphase = 0; iphase < 2; ++iphase)
    {
        for(size_t irank = 0; irank < bulk.mesh_meta_data().entity_rank_count(); ++irank)
        {
            stk::mesh::EntityRank rank = static_cast<stk::mesh::EntityRank>(irank);
            const stk::mesh::BucketVector& buckets = bulk.buckets(rank);
            buf.pack<size_t>(buckets.size());
        }
        allocate_or_communicate(iphase, comm);
    }
    int num_diffs = 0;
    return get_global_bucket_count_differences(env_data.m_worldComm, num_diffs);
}

std::string create_string_from_parts(const stk::mesh::PartVector& parts)
{
    std::string names;
    for(size_t i=0; i<parts.size(); ++i)
    {
        names += parts[i]->name() + " ";
    }
    return names;
}

void pack_buckets_parts(const stk::mesh::BucketVector& buckets, stk::CommBuffer &buff)
{
    for(size_t i = 0; i < buckets.size(); ++i)
    {
        const stk::mesh::PartVector& parts = buckets[i]->supersets();
        std::string part_names_for_bucket = create_string_from_parts(parts);
        pack_string(buff, part_names_for_bucket);
    }
}


bool bucket_part_memberships_match(const stk::mesh::BulkData& bulk, stk::EnvData& env_data)
{
    int numGlobalDiffs = bucket_counts_match(bulk, env_data);
    for(size_t irank = 0; irank < bulk.mesh_meta_data().entity_rank_count(); ++irank)
    {
        stk::CommSparse comm(env_data.m_worldComm);
        stk::mesh::EntityRank rank = static_cast<stk::mesh::EntityRank>(irank);
        const stk::mesh::BucketVector& buckets = bulk.buckets(rank);
        for(int iphase = 0; iphase < 2; ++iphase)
        {
            pack_buckets_parts(buckets, get_comm_buffer_for_destination_proc(comm));
            allocate_or_communicate(iphase, comm);
        }
    }
    numGlobalDiffs += get_global_bucket_part_membership_differences(env_data.m_worldComm, 0);
    return numGlobalDiffs == 0;
}

std::vector<size_t> unpack_bucket_counts(stk::CommBuffer& buf)
{
    std::vector<size_t> bucket_counts;
    while(buf.remaining())
    {
        size_t count = 0;
        buf.unpack<size_t>(count);
        bucket_counts.push_back(count);
    }
    return bucket_counts;
}

int compare_bucket_counts_message(int procId, const std::vector<size_t>& bucket_counts1, const std::vector<size_t>& bucket_counts2, std::ostringstream& globalss)
{
    std::ostringstream localss;
    int num_diffs = 0;
    if (bucket_counts1.size() != bucket_counts2.size())
    {
        num_diffs = 1;
        localss << "Number of entity ranks doesn't match: " << bucket_counts1.size() << " != " << bucket_counts2.size() << std::endl;
    }
    else
    {
        for(size_t i=0; i<bucket_counts1.size(); ++i)
        {
            if (bucket_counts1[i] != bucket_counts2[i])
            {
                num_diffs++;
                localss << "Number of buckets for entity rank "<< i << " doesn't match: " << bucket_counts1[i] << " != " << bucket_counts2[i] << std::endl;
            }
        }
    }
    if (num_diffs > 0)
    {
        globalss << localss.str();
    }
    return num_diffs;
}

bool compare_bucket_part_names_message(size_t num_entity_ranks1, size_t num_entity_ranks2, stk::EnvData& env_data, std::ostringstream& globalss)
{
    size_t num_equal_counts = std::min(num_entity_ranks1, num_entity_ranks2);
    std::pair<int,int> sourceProcs = getSourceProcs(env_data.m_worldComm);
    bool diffs_occurred = false;

    for(size_t i=0;i<num_equal_counts;++i)
    {
        std::ostringstream os;
        stk::CommSparse comm(env_data.m_worldComm);
        wait_for_apps_to_pack_and_send_data(comm);

        std::vector<std::string> part_names1_for_buckets_of_rank= get_sorted_part_names(comm.recv_buffer(sourceProcs.first));
        std::vector<std::string> part_names2_for_buckets_of_rank = get_sorted_part_names(comm.recv_buffer(sourceProcs.second));

        os << "Buckets in first app, not in second app, with parts: ";
        int num_local_diffs = part_name_set_difference_message(part_names1_for_buckets_of_rank, part_names2_for_buckets_of_rank, os);

        os << "Buckets in second app, not in first app, with parts: ";
        num_local_diffs += part_name_set_difference_message(part_names2_for_buckets_of_rank, part_names1_for_buckets_of_rank, os);

        if (num_local_diffs > 0)
        {
            globalss << os.str();
            diffs_occurred = true;
        }
    }

    return diffs_occurred;
}

int compare_bucket_counts(stk::EnvData& env_data, size_t& num_ranks1, size_t& num_ranks2, std::ostringstream& globalss)
{
    stk::CommSparse comm(env_data.m_worldComm);
    wait_for_apps_to_pack_and_send_data(comm);
    std::pair<int,int> sourceProcs = getSourceProcs(env_data.m_worldComm);
    std::vector<size_t> bucket_counts1 = unpack_bucket_counts(comm.recv_buffer(sourceProcs.first));
    std::vector<size_t> bucket_counts2 = unpack_bucket_counts(comm.recv_buffer(sourceProcs.second));
    num_ranks1 = bucket_counts1.size();
    num_ranks2 = bucket_counts2.size();
    int localProcId = stk::parallel_machine_rank(env_data.m_parallelComm);
    int num_diffs = compare_bucket_counts_message(localProcId, bucket_counts1, bucket_counts2, globalss);
    return get_global_bucket_count_differences(env_data.m_worldComm, num_diffs);
}

void build_message_for_missing_ranks(stk::EnvData& env_data, size_t num_ranks1, size_t num_ranks2, std::ostringstream& globalss)
{
    std::pair<int,int> sourceProcs = getSourceProcs(env_data.m_worldComm);
    int which_app = -1;
    size_t max = 0;
    int which_proc = -1;
    size_t num_equal_counts = std::min(num_ranks1, num_ranks2);
    if(num_ranks1 != num_equal_counts)
    {
        which_app = 0;
        max = num_ranks1;
        which_proc = sourceProcs.first;
    }
    else if (num_ranks2 != num_equal_counts)
    {
        which_app = 1;
        max = num_ranks2;
        which_proc = sourceProcs.second;
    }

    if(which_app != -1)
    {
        for(size_t i=num_equal_counts;i<max;++i)
        {
            stk::CommSparse comm(env_data.m_worldComm);
            wait_for_apps_to_pack_and_send_data(comm);

            stk::CommBuffer & buf = comm.recv_buffer(which_proc);
            std::vector<std::string> part_names_for_buckets_of_rank = get_sorted_part_names(buf);

            globalss << "Proc " << which_proc << " of app # " << which_app << " has extra entity ranks\n";
        }
    }
}

int bool_to_int(bool flag)
{
    return flag ? 1 : 0;
}

void compare_bucket_part_memberships(stk::EnvData& env_data, std::ostringstream& globalss)
{
    size_t num_ranks1 = 0, num_ranks2 = 0;
    compare_bucket_counts(env_data, num_ranks1, num_ranks2, globalss);
    bool diffs_exist = compare_bucket_part_names_message(num_ranks1, num_ranks2, env_data, globalss);
    if (diffs_exist)
    {
        build_message_for_missing_ranks(env_data, num_ranks1, num_ranks2, globalss);
    }
    get_global_bucket_part_membership_differences(env_data.m_worldComm, bool_to_int(diffs_exist));
}

void make_part_membership_different_on_app1(stk::mesh::BulkData& bulk)
{
    stk::mesh::Part& newPart = bulk.mesh_meta_data().declare_part("new_part");
    stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
    bulk.modification_begin();
    if (bulk.is_valid(elem1) && bulk.bucket(elem1).owned())
    {
        bulk.change_entity_parts(elem1, {&newPart}, {});
    }
    bulk.modification_end();
}

TEST(ParallelDebugTool, mock_app_1)
{
    int color = 0;
    stk::EnvData & env_data = stk::EnvData::instance();
    setup_env_data(env_data, MPI_COMM_WORLD, color);
    expect_mpmd(env_data);

    TestMesh mesh(env_data.m_parallelComm);
    stk::mesh::BulkData & bulk = mesh.get_bulk_data();
    expect_num_elements(bulk);

    EXPECT_TRUE(parts_match(bulk, env_data));

    make_part_membership_different_on_app1(bulk);
    EXPECT_TRUE(!bucket_part_memberships_match(bulk, env_data));
}

TEST(ParallelDebugTool, mock_app_2)
{
    int color = 1;
    stk::EnvData & env_data = stk::EnvData::instance();
    setup_env_data(env_data, MPI_COMM_WORLD, color);
    expect_mpmd(env_data);

    TestMesh mesh(env_data.m_parallelComm);
    stk::mesh::BulkData & bulk = mesh.get_bulk_data();
    expect_num_elements(bulk);

    EXPECT_TRUE(parts_match(bulk, env_data));

    EXPECT_TRUE(!bucket_part_memberships_match(bulk, env_data));
}

std::vector<std::string> get_delimited_strings_from_stream(std::stringstream& ss, std::string& scratch_space, char delim)
{
    std::vector<std::string> elems;
    while (std::getline(ss, scratch_space, delim))
    {
        elems.push_back(scratch_space);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim)
{
    std::stringstream ss(s);
    std::string scratch_space;
    return get_delimited_strings_from_stream(ss, scratch_space, delim);
}

void output_message(const std::string &message, int localProcId)
{
    if (!message.empty())
    {
        std::ostringstream massaged_output;
        std::vector<std::string> tokens = split(message,'\n');
        for(std::string& s: tokens)
        {
            massaged_output << "[" << localProcId << "] " << s << std::endl;
        }
        std::cerr << massaged_output.str();
    }
}

void run_diffing_tool(stk::EnvData &env_data)
{
    std::ostringstream globalss;
    compare_parts(env_data, globalss);
    compare_bucket_part_memberships(env_data, globalss);
    output_message(globalss.str(), stk::parallel_machine_rank(env_data.m_parallelComm));
}

TEST(ParallelDebugTool, mockDiffingToolUnit)
{
    int color = 2;
    stk::EnvData &env_data = stk::EnvData::instance();
    setup_env_data(env_data, MPI_COMM_WORLD, color);
    expect_mpmd(env_data);
    run_diffing_tool(env_data);
}

TEST(ParallelDebugTool, mockDiffingToolApp)
{
    int color = 2;
    stk::EnvData &env_data = stk::EnvData::instance();
    setup_env_data(env_data, MPI_COMM_WORLD, color);
    expect_mpmd(env_data);
    //run_diffing_tool(env_data);
}

}
