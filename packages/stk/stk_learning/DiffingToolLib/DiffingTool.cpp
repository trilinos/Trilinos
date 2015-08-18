#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <DiffingToolLib/DiffingTool.hpp>
#include <stk_util/diag/StringUtil.hpp>

namespace stk {
namespace diff {

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

void pack_part_names_except(stk::CommBuffer& buf, const stk::mesh::PartVector& parts, stk::mesh::Part* skipPart)
{
    for(size_t i=0; i<parts.size(); ++i)
    {
        if (skipPart==nullptr || parts[i]->mesh_meta_data_ordinal() != skipPart->mesh_meta_data_ordinal())
        {
            pack_string(buf, parts[i]->name());
        }
    }
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

void send_part_names_to_diffing_tool_except(const stk::mesh::BulkData& bulk, stk::ParallelMachine communicator, stk::mesh::Part* skipPart)
{
    stk::CommSparse comm(communicator);
    for(int iphase = 0; iphase < 2; ++iphase)
    {
        pack_part_names_except(get_comm_buffer_for_destination_proc(comm), bulk.mesh_meta_data().get_parts(), skipPart);
        allocate_or_communicate(iphase, comm);
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

int get_global_part_differences_for_app(stk::ParallelMachine comm)
{
    int numLocalDiffs = 0;
    return get_global_part_differences(comm, numLocalDiffs);
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

bool parts_match(const stk::mesh::BulkData& bulk, stk::EnvData& env_data)
{
    send_part_names_to_diffing_tool(bulk, env_data.m_worldComm);
    return get_global_part_differences_for_app(env_data.m_worldComm) == 0;
}

bool parts_match_except(const stk::mesh::BulkData& bulk, stk::EnvData& env_data, stk::mesh::Part* skipPart)
{
    send_part_names_to_diffing_tool_except(bulk, env_data.m_worldComm, skipPart);
    return get_global_part_differences_for_app(env_data.m_worldComm) == 0;
}

int get_global_bucket_count_differences(stk::ParallelMachine comm, int numLocalDiffs)
{
    return parallel_sum(comm, numLocalDiffs);
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
        stk::diff::pack_string(buff, part_names_for_bucket);
    }
}

int get_global_bucket_part_membership_differences(stk::ParallelMachine comm, int numLocalDiffs)
{
    return parallel_sum(comm, numLocalDiffs);
}

bool bucket_part_memberships_match(const stk::mesh::BulkData& bulk, stk::EnvData& env_data)
{
    int numGlobalDiffs = bucket_counts_match(bulk, env_data);

    if(numGlobalDiffs > 0)
    {
        for(size_t irank = 0; irank < bulk.mesh_meta_data().entity_rank_count(); ++irank)
        {
            stk::CommSparse comm(env_data.m_worldComm);
            stk::CommBuffer &buff = stk::diff::get_comm_buffer_for_destination_proc(comm);
            stk::mesh::EntityRank rank = static_cast<stk::mesh::EntityRank>(irank);
            stk::mesh::EntityVector entities;
            stk::mesh::get_entities(bulk, rank, entities);
            for(int iphase = 0; iphase < 2; ++iphase)
            {
                for(size_t i=0;i<entities.size();++i)
                {
                    const stk::mesh::PartVector& parts = bulk.bucket(entities[i]).supersets();
                    std::string part_names_for_entity = create_string_from_parts(parts);
                    std::string string_to_send;
                    if(irank != 1 && irank != 2)
                    {
                        string_to_send = sierra::to_string(bulk.identifier(entities[i])) + " " + part_names_for_entity;
                    }
                    else
                    {
                        string_to_send = part_names_for_entity;
                    }
                    stk::diff::pack_string(buff, string_to_send);
                }
                stk::diff::allocate_or_communicate(iphase, comm);
            }
        }
    }

    for(size_t irank = 0; irank < bulk.mesh_meta_data().entity_rank_count(); ++irank)
    {
        stk::CommSparse comm(env_data.m_worldComm);
        stk::mesh::EntityRank rank = static_cast<stk::mesh::EntityRank>(irank);
        const stk::mesh::BucketVector& buckets = bulk.buckets(rank);
        for(int iphase = 0; iphase < 2; ++iphase)
        {
            pack_buckets_parts(buckets, stk::diff::get_comm_buffer_for_destination_proc(comm));
            stk::diff::allocate_or_communicate(iphase, comm);
        }
    }
    numGlobalDiffs += get_global_bucket_part_membership_differences(env_data.m_worldComm, 0);
    return numGlobalDiffs == 0;
}

void communicate_run_state(stk::EnvData& env_data, bool continue_runs)
{
    int send_int = 0;
    if(continue_runs)
    {
        send_int = 1;
    }
    stk::CommSparse comm(env_data.m_worldComm);
    for(int iphase = 0; iphase < 2; ++iphase)
    {
        get_comm_buffer_for_destination_proc(comm).pack<int>(send_int);
        allocate_or_communicate(iphase, comm);
    }
}

}
}
