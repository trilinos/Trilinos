#include "BulkDataTester.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity, operator<<, etc
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, etc
#include <stk_util/parallel/CommSparse.hpp>  // for CommSparse


namespace stk { namespace mesh { namespace unit_test {

bool BulkDataTester::is_entity_in_ghosting_comm_map(stk::mesh::Entity entity)
{
    stk::mesh::EntityKey entityKey = this->entity_key(entity);
    bool is_entity_in_aura_comm_map = !this->internal_entity_comm_map(entityKey, this->aura_ghosting()).empty();
    return is_entity_in_aura_comm_map;
}

void BulkDataTester::check_sharing_comm_maps()
{
    stk::CommSparse comm(parallel());

    for(int phase = 0; phase < 2; ++phase)
    {
        for(stk::mesh::EntityCommListInfoVector::const_iterator i = this->my_internal_comm_list().begin(); i != this->my_internal_comm_list().end(); ++i)
        {
            for(stk::mesh::PairIterEntityComm ec = this->internal_entity_comm_map(i->key); !ec.empty(); ++ec)
            {
                int type = ec->ghost_id;
                if ( type == 0 )
                {
                    std::vector<int> sharingProcs;
                    this->comm_shared_procs(i->key, sharingProcs);
                    // pack shared info
                    int owner = -1;
                    if(bucket_ptr(i->entity) != 0)
                    {
                        owner = parallel_owner_rank(i->entity);
                    }

                    comm.send_buffer(ec->proc).pack<stk::mesh::EntityKey>(i->key).pack<int>(type).pack<int>(owner);
                    comm.send_buffer(ec->proc).pack<size_t>(sharingProcs.size());
                    for (size_t proc=0;proc<sharingProcs.size();++proc)
                    {
                        comm.send_buffer(ec->proc).pack<int>(sharingProcs[proc]);
                    }
                }
            }
        }

        if(phase == 0)
        {
            comm.allocate_buffers();
        }
        else
        {
            comm.communicate();
        }
    }

    // unpack

    std::ostringstream os;
    bool anyErrors = false;

    for(int i = 0; i < parallel_size(); ++i)
    {
        if ( i != parallel_rank() )
        {
            stk::mesh::EntityKey key;
            int type = -1;
            int from = -1;
            int owner = -1;
            while(comm.recv_buffer(i).remaining())
            {
                from = i;
                comm.recv_buffer(from).unpack<stk::mesh::EntityKey>(key).unpack<int>(type).unpack<int>(owner);

                size_t numSharingProcs = 0;
                comm.recv_buffer(from).unpack<size_t>(numSharingProcs);
                std::vector<int> sharingProcs(numSharingProcs);
                for (size_t proc=0;proc<numSharingProcs;++proc)
                {
                    comm.recv_buffer(from).unpack<int>(sharingProcs[proc]);
                }

                std::vector<int> localSharingProcs;
                this->comm_shared_procs(key, localSharingProcs);

                std::sort(localSharingProcs.begin(), localSharingProcs.end());
                std::sort(sharingProcs.begin(), sharingProcs.end());
                size_t maxNum = localSharingProcs.size() + sharingProcs.size();
                std::vector<int> unsharedProcs(maxNum);
                std::vector<int>::iterator iter = std::set_symmetric_difference( localSharingProcs.begin(), localSharingProcs.end(),
                        sharingProcs.begin(), sharingProcs.end(), unsharedProcs.begin());

                size_t numUnshared = iter - unsharedProcs.begin();
                unsharedProcs.resize(numUnshared);

                int counter = 0;
                {
                    for (size_t j=0;j<numUnshared;++j)
                    {
                        std::vector<int>::iterator localIter = std::find(sharingProcs.begin(), sharingProcs.end(), unsharedProcs[j]);
                        if ( localIter != sharingProcs.end() && *localIter != parallel_rank() )
                        {
                            if ( counter == 0 )
                            {
                                os << "Error in sharing between procs for entity " << key.id() << " with rank " << key.rank()  << "  between procs: " << this->parallel_rank() << " and " << from << std::endl;
                                counter++;
                            }
                            os << "\tEntity " << key << " is shared with proc " << unsharedProcs[j] << " from other proc: "
                                    << from << " but not from this proc: " << parallel_rank() << std::endl;
                            anyErrors = true;
                        }

                        localIter = std::find(localSharingProcs.begin(), localSharingProcs.end(), unsharedProcs[j]);
                        if ( localIter != localSharingProcs.end() && *localIter != from )
                        {
                            if ( counter == 0 )
                            {
                                os << "Error in sharing between procs for entity " << key.id() << " with rank " << key.rank()  << "  between procs: " << this->parallel_rank() << " and " << from << std::endl;
                                counter++;
                            }
                            os << "\tEntity " << key << " is shared with proc " << unsharedProcs[j] << " from this proc: "
                                    << parallel_rank() << " but not from other proc: " << from << std::endl;
                            anyErrors = true;
                        }
                    }
                }
            }
        }
    }

    ThrowRequireMsg(!anyErrors, os.str());
}

} } } // namespace stk mesh unit_test
