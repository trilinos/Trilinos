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

#ifndef stk_mesh_ModificationSummary_hpp
#define stk_mesh_ModificationSummary_hpp

#include <fstream>   // for writing file
#include <algorithm> // for sort
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/EntityKey.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/Ghosting.hpp"

namespace stk { namespace mesh { class BulkData; } }

namespace stk
{

class EmptyModificationSummary
{
public:
    EmptyModificationSummary(stk::mesh::BulkData& bulkData)
    {
    }

    ~EmptyModificationSummary(){}

    // void track_create_ghosting();
    void track_induced_parts(stk::mesh::Entity entity, stk::mesh::Entity e_to, const stk::mesh::OrdinalVector& add_parts, const stk::mesh::OrdinalVector& emptyParts)
    {

    }

    void track_change_ghosting(const stk::mesh::Ghosting & ghosts, const std::vector<stk::mesh::EntityProc> & add_send , const std::vector<stk::mesh::EntityKey> & remove_receive )
    {
    }

    void track_add_to_ghosting(const stk::mesh::Ghosting & ghosts, const std::vector<stk::mesh::EntityProc> & add_send )
    {
    }

    void track_destroy_relation(stk::mesh::Entity e_from, stk::mesh::Entity e_to, stk::mesh::RelationIdentifier rel)
    {
    }

    void track_declare_relation(stk::mesh::Entity e_from, stk::mesh::Entity e_to, stk::mesh::RelationIdentifier rel, stk::mesh::Permutation permut)
    {
    }

    void track_declare_entity(stk::mesh::EntityRank rank, stk::mesh::EntityId newId, const stk::mesh::PartVector& addParts)
    {
    }

    void track_change_entity_owner(const std::vector<stk::mesh::EntityProc> &changes)
    {
    }

    void track_set_global_id(stk::mesh::Entity entity, int64_t newId)
    {
    }

    void track_change_entity_id(stk::mesh::EntityId newId, stk::mesh::Entity entity)
    {
    }

    void track_destroy_entity(stk::mesh::Entity entity)
    {
    }

    void track_change_entity_parts(stk::mesh::Entity entity, const stk::mesh::OrdinalVector& addParts, const stk::mesh::OrdinalVector& rmParts)
    {
    }

    void track_comm_map_insert(stk::mesh::Entity entity, const stk::mesh::EntityCommInfo & val)
    {

    }

    void track_comm_map_erase(stk::mesh::EntityKey key, const stk::mesh::EntityCommInfo & val)
    {

    }

    void track_comm_map_erase(stk::mesh::EntityKey key, const stk::mesh::Ghosting & val)
    {

    }

    void track_comm_map_clear_ghosting(stk::mesh::EntityKey key)
    {

    }

    void track_comm_map_clear(stk::mesh::EntityKey key)
    {

    }

    void track_set_parallel_owner_rank_but_not_comm_lists(stk::mesh::Entity entity, int old_owner, int new_owner)
    {

    }

    void track_change_owner_in_comm_data(stk::mesh::EntityKey key, int old_owner, int new_owner)
    {

    }

    void write_summary(int mod_cycle_count, bool sort=true)
    {
    }

    void set_proc_id(int proc_id)
    {
    }
};

class ModificationSummary
{
public:
    ModificationSummary(stk::mesh::BulkData& bulkData);

    ~ModificationSummary();

    // void track_create_ghosting();
    void track_induced_parts(stk::mesh::Entity entity, stk::mesh::Entity e_to, const stk::mesh::OrdinalVector& add_parts, const stk::mesh::OrdinalVector& emptyParts);

    void track_change_ghosting(const stk::mesh::Ghosting & ghosts, const std::vector<stk::mesh::EntityProc> & add_send , const std::vector<stk::mesh::EntityKey> & remove_receive );

    void track_add_to_ghosting(const stk::mesh::Ghosting & ghosts, const std::vector<stk::mesh::EntityProc> & add_send );

    void track_destroy_relation(stk::mesh::Entity e_from, stk::mesh::Entity e_to, stk::mesh::RelationIdentifier rel);

    void track_declare_relation(stk::mesh::Entity e_from, stk::mesh::Entity e_to, stk::mesh::RelationIdentifier rel, stk::mesh::Permutation permut);

    void track_declare_entity(stk::mesh::EntityRank rank, stk::mesh::EntityId newId, const stk::mesh::PartVector& addParts);

    void track_change_entity_owner(const std::vector<stk::mesh::EntityProc> &changes);

    void track_change_entity_id(stk::mesh::EntityId newId, stk::mesh::Entity entity);

    void track_set_global_id(stk::mesh::Entity entity, uint32_t newId);

    void track_destroy_entity(stk::mesh::Entity entity);

    void track_change_entity_parts(stk::mesh::Entity entity, const stk::mesh::OrdinalVector& addParts, const stk::mesh::OrdinalVector& rmParts);

    void track_comm_map_insert(stk::mesh::Entity entity, const stk::mesh::EntityCommInfo & val);

    void track_comm_map_erase(stk::mesh::EntityKey key, const stk::mesh::EntityCommInfo & val);

    void track_comm_map_erase(stk::mesh::EntityKey key, const stk::mesh::Ghosting & val);

    void track_comm_map_clear_ghosting(stk::mesh::EntityKey key);

    void track_comm_map_clear(stk::mesh::EntityKey key);

    void track_set_parallel_owner_rank_but_not_comm_lists(stk::mesh::Entity entity, int old_owner, int new_owner);

    void track_change_owner_in_comm_data(stk::mesh::EntityKey key, int old_owner, int new_owner);

    void write_summary(int mod_cycle_count, bool sort=true);

    void set_proc_id(int proc_id) { m_procId = proc_id; }

private:

    void clear_summary();

    bool isValid(stk::mesh::Entity entity) const;

    stk::mesh::EntityKey getEntityKey(stk::mesh::Entity entity) const;

    void addEntityKeyAndStringToTracker(stk::mesh::EntityKey key, const std::string& string);

    std::string get_filename(int mod_cycle_count) const;

    int my_proc_id() const;

    void writeParts(std::ostringstream& os, const std::string &label, const stk::mesh::OrdinalVector& parts);

    stk::mesh::BulkData &m_bulkData;
    std::vector<std::pair<stk::mesh::EntityKey, std::string> > m_stringTracker;
    int m_lastModCycle;
    int m_modCounter;
    int m_modificationSummaryNumber;
    int m_procId = -1;
    std::vector<size_t> watchedFaces;
};

} // namespace

#endif
