// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "BulkDataTester.hpp"
#include <algorithm>                                 // for find, sort, etc
#include <ostream>                                   // for operator<<, etc
#include <stk_unit_test_utils/FaceTestingUtils.hpp>
#include <stk_util/parallel/CommSparse.hpp>          // for CommSparse
#include <vector>                                    // for vector, etc
#include "stk_mesh/base/BulkData.hpp"                // for BulkData, etc
#include "stk_mesh/base/Entity.hpp"                  // for Entity
#include "stk_mesh/base/EntityKey.hpp"               // for EntityKey, etc
#include "stk_mesh/base/Types.hpp"                   // for EntityVector, etc
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_util/util/ReportHandler.hpp"    // for ThrowRequireMsg, etc
#include "stk_util/parallel/ParallelComm.hpp"        // for CommBuffer
#include "stk_util/util/PairIter.hpp"                // for PairIter
namespace stk { namespace mesh { class Part; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################


namespace stk { namespace unit_test_util {

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

///////////////////////////

void BulkDataFaceSharingTester::change_entity_key_and_nodes(const std::vector<stk::mesh::shared_entity_type> & potentially_shared_sides)
{
    for(size_t i = 0, e = potentially_shared_sides.size(); i < e; ++i)
    {
        stk::mesh::Entity entity = potentially_shared_sides[i].entity;
        if(potentially_shared_sides[i].need_update_nodes)
        {
            if(potentially_shared_sides[i].global_key != potentially_shared_sides[i].local_key)
            {
                my_internal_change_entity_key(potentially_shared_sides[i].local_key, potentially_shared_sides[i].global_key, entity);
            }
            change_connectivity_for_edge_or_face(entity, potentially_shared_sides[i].nodes);
        }
    }
}

void BulkDataFaceSharingTester::change_entity_key_and_update_sharing_info(std::vector<stk::mesh::shared_entity_type> & potentially_shared_sides)
{
   change_entity_key_and_nodes(potentially_shared_sides);
   insert_sharing_info_into_comm_map(potentially_shared_sides);
}

void BulkDataFaceSharingTester::update_shared_entity_this_proc2(stk::mesh::EntityKey global_key_other_proc, stk::mesh::shared_entity_type& shared_entity, int proc_id, const std::vector<stk::mesh::EntityKey>& nodes)
{
    stk::mesh::Entity entity = shared_entity.entity;
    shared_entity.sharing_procs.push_back(proc_id);
    if(proc_id < this->parallel_rank())
    {
        shared_entity.global_key = global_key_other_proc;
        shared_entity.nodes = nodes;
        shared_entity.need_update_nodes = true;
    }
    this->internal_mark_entity(entity, BulkData::IS_SHARED);
}

stk::mesh::EntityVector get_elements_connected_to_all_nodes(const BulkDataFaceSharingTester& bulkData, const stk::mesh::EntityVector& nodes)
{
    stk::mesh::EntityVector elements;
    stk::mesh::impl::find_locally_owned_elements_these_nodes_have_in_common(bulkData, nodes.size(), nodes.data(), elements);
    return elements;
}

void add_side_to_shared_entities(stk::mesh::Entity side, std::vector<stk::mesh::shared_entity_type>& shared_entities_this_proc, const stk::mesh::shared_entity_type &shared_entity_other_proc, int other_proc_id)
{
    stk::mesh::shared_entity_type sentity(shared_entity_other_proc.global_key, side, shared_entity_other_proc.topology);
    sentity.nodes = shared_entity_other_proc.nodes;
    sentity.sharing_procs.push_back(other_proc_id);
    shared_entities_this_proc.push_back(sentity);
}

stk::mesh::Entity create_side_and_add_to_shared_entity_list(stk::mesh::Entity element, const stk::mesh::EntityVector& nodes, const stk::mesh::shared_entity_type &shared_entity_other_proc,
        stk::mesh::BulkData& bulkData, std::vector<stk::mesh::shared_entity_type>& shared_entities_this_proc, int other_proc_id, stk::mesh::Part& root_topo_part)
{
    stk::mesh::Entity side = stk::unit_test_util::declare_element_side_with_nodes(bulkData, element, nodes, shared_entity_other_proc.global_key.id(), root_topo_part);
    ThrowRequireWithSierraHelpMsg(bulkData.is_valid(side));
    add_side_to_shared_entities(side, shared_entities_this_proc, shared_entity_other_proc, other_proc_id);
    return side;
}

void BulkDataFaceSharingTester::connect_side_from_other_proc_to_local_elements(const stk::mesh::EntityVector& elements, const stk::mesh::EntityVector& nodes, const stk::mesh::shared_entity_type &shared_entity_other_proc,
        stk::mesh::BulkData& bulkData, std::vector<stk::mesh::shared_entity_type>& shared_entities_this_proc, int other_proc_id)
{
    for(stk::mesh::Entity element : elements)
    {
        stk::mesh::Entity side = create_side_and_add_to_shared_entity_list(element, nodes, shared_entity_other_proc, *this,
                shared_entities_this_proc, other_proc_id, get_topology_root_part(shared_entity_other_proc.topology));
        this->internal_mark_entity(side, BulkData::IS_SHARED);
    }
}

void BulkDataFaceSharingTester::create_and_connect_shared_face_on_this_proc(const stk::mesh::shared_entity_type &shared_entity_other_proc, std::vector<stk::mesh::shared_entity_type>& shared_entities_this_proc, int other_proc_id)
{
    stk::mesh::EntityVector nodes = stk::mesh::impl::convert_keys_to_entities(*this, shared_entity_other_proc.nodes);
    stk::mesh::EntityVector elements = get_elements_connected_to_all_nodes(*this, nodes);
    connect_side_from_other_proc_to_local_elements(elements, nodes, shared_entity_other_proc, *this, shared_entities_this_proc, other_proc_id);
}

void BulkDataFaceSharingTester::check_if_entity_from_other_proc_exists_on_this_proc_and_update_info_if_shared(std::vector<stk::mesh::shared_entity_type>& shared_entities_this_proc,
        int proc_id, const stk::mesh::shared_entity_type &shared_entity_other_proc)
{
    int matching_index = does_entity_exist_in_list(shared_entities_this_proc, shared_entity_other_proc);
    bool entitiesAreTheSame = matching_index >= 0;
    if( entitiesAreTheSame )
        update_shared_entity_this_proc2(shared_entity_other_proc.global_key, shared_entities_this_proc[matching_index], proc_id, shared_entity_other_proc.nodes);
    else
        create_and_connect_shared_face_on_this_proc(shared_entity_other_proc, shared_entities_this_proc, proc_id);
}

} } // namespace stk unit_test_util
