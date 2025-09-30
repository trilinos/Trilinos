#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/EntityCommListInfo.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_mesh/baseImpl/check_comm_list.hpp>

#include <iostream>
#include <vector>
#include <algorithm>

namespace stk {
namespace mesh {
namespace impl {

struct KeyProcGhostId {
    EntityKey key;
    int proc;
    unsigned ghostId;
    bool alreadyFound;

    bool operator<(const KeyProcGhostId& other) const {
        if (key != other.key) { return key < other.key; }
        if (proc != other.proc) { return proc < other.proc; }
        if (ghostId != other.ghostId) { return ghostId < other.ghostId; }
        return false;
    }

    bool operator==(const KeyProcGhostId& other) const {
        return key == other.key && proc == other.proc && ghostId == other.ghostId;
    }
};

void pack_key_and_ghost_id(EntityKey key, unsigned ghost_id, int proc, stk::CommSparse& comm)
{
    stk::CommBuffer& buf = comm.send_buffer(proc);
    buf.pack<EntityKey>(key);
    buf.pack<unsigned>(ghost_id);
}

void pack_key_and_ghost_ids(EntityKey key, bool locally_owned, PairIterEntityComm commInfo,
                            stk::CommSparse& comm)
{
    for(; !commInfo.empty(); ++commInfo) {
        const bool needToSend = locally_owned || (commInfo->ghost_id == 0);
        if (needToSend) {
            pack_key_and_ghost_id(key, commInfo->ghost_id, commInfo->proc, comm);
        }
    }
}

void pack_send_data(const stk::mesh::BulkData& mesh, int local_proc,
                    const EntityCommDatabase& commDB,
                    const EntityCommListInfoVector& comm_list, stk::CommSparse& comm)
{
    for(const EntityCommListInfo& commInfo : comm_list) {
        STK_ThrowAssert(commInfo.entity_comm != -1);
        STK_ThrowAssert(commInfo.entity_comm == commDB.entity_comm(mesh.entity_key(commInfo.entity)));
        pack_key_and_ghost_ids(commInfo.key, mesh.parallel_owner_rank(commInfo.entity)==local_proc,
                               commDB.comm(commInfo.entity_comm), comm);
    }
}

void push_back_key_proc_ghost_id(EntityKey key, int proc, unsigned ghost_id, std::vector<KeyProcGhostId>& data_vec)
{
    KeyProcGhostId data = {key, proc, ghost_id, false};
    data_vec.push_back(data);
}

void push_back_key_procs_ghost_ids(EntityKey key, bool locally_owned,
                                   PairIterEntityComm commInfo,
                                   std::vector<KeyProcGhostId>& key_proc_ghostid_vec)
{
    for(; !commInfo.empty(); ++commInfo) {
        const bool expectToRecv = (!locally_owned || (commInfo->ghost_id == 0))
                                && commInfo->ghost_id < BulkData::SYMM_INFO;
        if (expectToRecv) {
            push_back_key_proc_ghost_id(key, commInfo->proc, commInfo->ghost_id, key_proc_ghostid_vec);
        }
    }
}

void fill_expected_recv_data(const stk::mesh::BulkData& mesh,
                             const EntityCommDatabase& commDB,
                             const EntityCommListInfoVector& comm_list,
                             std::vector<KeyProcGhostId>& recv_data)
{
    for(const EntityCommListInfo& commInfo : comm_list) {
        STK_ThrowAssertMsg(commInfo.entity_comm != -1,"P"<<mesh.parallel_rank()<<" entity_comm==-1 for "<<commInfo.key);
        STK_ThrowAssertMsg(commInfo.entity_comm == commDB.entity_comm(mesh.entity_key(commInfo.entity)),
                           "P"<<mesh.parallel_rank()<<" mod-cycle="<<mesh.synchronized_count()
                           <<" bad comm-map data, entity-key="<<mesh.entity_key(commInfo.entity)<<", key="<<commInfo.key
                           <<", commInfo.entity_comm="<<commInfo.entity_comm<<", commDB.entity_comm(entity-key)="<<commDB.entity_comm(mesh.entity_key(commInfo.entity)));
        push_back_key_procs_ghost_ids(commInfo.key,
                                      mesh.parallel_owner_rank(commInfo.entity)==mesh.parallel_rank(),
                                      commDB.comm(commInfo.entity_comm), recv_data);
    }
    stk::util::sort_and_unique(recv_data);
}

void check_recvd_key_and_ghost_id_against_expected_recv_data(stk::CommBuffer& buf, int proc, std::vector<KeyProcGhostId>& expected_recv_data, std::ostream& os)
{
    KeyProcGhostId data;
    data.proc = proc;
    buf.unpack<EntityKey>(data.key);
    buf.unpack<unsigned>(data.ghostId);

    std::vector<KeyProcGhostId>::iterator iter = std::lower_bound(expected_recv_data.begin(), expected_recv_data.end(), data);

    const bool found = (iter != expected_recv_data.end() && *iter == data);
    if (found) {
        iter->alreadyFound = true;
    }
    else {
        os << "\trecvd {"<<data.key<<",ghost_id="<<data.ghostId<<"} from proc "<<proc<<" but not in comm list"<<std::endl;
    }
}

void check_for_expected_recv_data_that_failed_to_arrive(const std::vector<KeyProcGhostId>& expected_recv_data, int localProc, std::ostream& os)
{
    if (!expected_recv_data.empty()) {
        for(const KeyProcGhostId& data : expected_recv_data) {
           if (!data.alreadyFound) {
              os << "\tP"<<localProc<<" Failed to recv {"<<data.key<<",ghost_id="<<data.ghostId<<"} from proc "<<data.proc<<"\n";
           }
        }
    }
}

void pack_and_send_comm_list_data(const stk::mesh::BulkData& mesh,
                                  const EntityCommDatabase& commDB,
                                  const EntityCommListInfoVector& comm_list, stk::CommSparse& comm)
{
    for(int phase=0; phase<2; ++phase) {
        pack_send_data(mesh, mesh.parallel_rank(), commDB, comm_list, comm);

        if (phase==0) {
            comm.allocate_buffers();
        }
        else {
            comm.communicate();
        }
    }
}

void unpack_and_check_recvd_data(stk::CommSparse& comm, int local_proc, int num_procs, std::vector<KeyProcGhostId>& expected_recv_data, std::ostream& os)
{
    for(int p=0; p<num_procs; ++p) {
        if (p != local_proc) {
            stk::CommBuffer& buf = comm.recv_buffer(p);
            while(buf.remaining()) {
                check_recvd_key_and_ghost_id_against_expected_recv_data(buf, p, expected_recv_data, os);
            }
        }
    }
}

bool is_comm_list_globally_consistent(const stk::mesh::BulkData& mesh,
                                      const EntityCommDatabase& commDB,
                                      const EntityCommListInfoVector& comm_list,
                                      std::ostream& error_msg)
{
    int local_proc = mesh.parallel_rank();
    int num_procs = mesh.parallel_size();

    std::vector<KeyProcGhostId> expected_recv_data;
    fill_expected_recv_data(mesh, commDB, comm_list, expected_recv_data);

    stk::CommSparse comm(mesh.parallel());
    pack_and_send_comm_list_data(mesh, commDB, comm_list, comm);

    std::ostringstream os;
    unpack_and_check_recvd_data(comm, local_proc, num_procs, expected_recv_data, os);

    check_for_expected_recv_data_that_failed_to_arrive(expected_recv_data, local_proc, os);

    std::string str = os.str();
    if (!str.empty()) {
        error_msg<<"P"<<local_proc<<" mod-cycle="<<mesh.synchronized_count()<<" check_comm_list_global_consistency:\n"<<str;
    }

    return str.empty();
}

} //namespace impl
} //namespace mesh
} //namespace stk

