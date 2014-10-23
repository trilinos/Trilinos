// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#ifndef _BulkDataTester_hpp_
#define _BulkDataTester_hpp_

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>      // for MeshIndex, EntityRank, etc

class BulkDataTester : public stk::mesh::BulkData
{
public:

    BulkDataTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm) :
            stk::mesh::BulkData(mesh_meta_data, comm)
    {
    }

    virtual ~BulkDataTester()
    {
    }

    void my_internal_change_entity_owner( const std::vector<stk::mesh::EntityProc> & arg_change, bool regenerate_aura = true, modification_optimization mod_optimization = MOD_END_SORT )
    {
        this->internal_change_entity_owner(arg_change,regenerate_aura,mod_optimization);
    }

    void internal_change_entity_owner_exp(const std::vector<stk::mesh::EntityProc> & arg_change,
            bool regenerate_aura = true,
            modification_optimization mod_optimization = MOD_END_SORT);

    void change_entity_owner_exp( const std::vector<stk::mesh::EntityProc> & arg_change,
                                                bool regenerate_aura = true,
                                                modification_optimization mod_optimization = MOD_END_SORT );

    void my_internal_resolve_shared_modify_delete()
    {
        this->internal_resolve_shared_modify_delete();
    }

    void my_resolve_ownership_of_modified_entities(const std::vector<stk::mesh::Entity> &shared_new)
    {
        this->resolve_ownership_of_modified_entities(shared_new);
    }

    void my_move_entities_to_proper_part_ownership( const std::vector<stk::mesh::Entity> &shared_modified )
    {
        this->move_entities_to_proper_part_ownership(shared_modified);
    }

//    void reset_closure_count(stk::mesh::Entity entity)
//    {
//        m_closure_count[entity.local_offset()] = 0;
//    }

    uint16_t closure_count(stk::mesh::Entity entity)
    {
        return m_closure_count[entity.local_offset()];
    }

    bool my_entity_comm_map_insert(stk::mesh::Entity entity, const stk::mesh::EntityCommInfo & val)
    {
        return BulkData::entity_comm_map_insert(entity, val);
    }

    bool my_entity_comm_map_erase(const stk::mesh::EntityKey& key, const stk::mesh::EntityCommInfo& commInfo)
    {
        return BulkData::entity_comm_map_erase(key, commInfo);
    }


    bool my_entity_comm_map_erase(const stk::mesh::EntityKey& key, const stk::mesh::Ghosting& ghost)
    {
        return BulkData::entity_comm_map_erase(key, ghost);
    }

    void my_entity_comm_map_clear(const stk::mesh::EntityKey& key)
    {
        BulkData::entity_comm_map_clear(key);
    }

    void my_entity_comm_map_clear_ghosting(const stk::mesh::EntityKey& key)
    {
        BulkData::entity_comm_map_clear_ghosting(key);
    }

    void my_update_comm_list_based_on_changes_in_comm_map()
    {
        this->update_comm_list_based_on_changes_in_comm_map();
    }

    void my_update_comm_list(const std::vector<stk::mesh::Entity>& shared_modified)
    {
        this->update_comm_list(shared_modified);
    }

    const stk::mesh::EntityCommDatabase my_entity_comm_map() const
    {
        return m_entity_comm_map;
    }

    void internal_resolve_parallel_create_exp();
    void internal_update_distributed_index_exp(std::vector<stk::mesh::Entity> & shared_new );

    bool internal_modification_end_for_change_entity_owner_exp( bool regenerate_aura, modification_optimization opt );

    stk::mesh::impl::BucketRepository& my_bucket_repository()
    {
        return this->bucket_repository();
    }
};


#endif
