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
    void track_induced_parts(stk::mesh::Entity entity, stk::mesh::Entity e_to, const stk::mesh::PartVector& add_parts, const stk::mesh::PartVector& emptyParts)
    {

    }

    void track_change_ghosting(const stk::mesh::Ghosting & ghosts, const std::vector<stk::mesh::EntityProc> & add_send , const std::vector<stk::mesh::EntityKey> & remove_receive )
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

    void track_change_entity_id(stk::mesh::EntityId newId, stk::mesh::Entity entity)
    {
    }

    void track_destroy_entity(stk::mesh::Entity entity)
    {
    }

    void track_change_entity_parts(stk::mesh::Entity entity, const stk::mesh::PartVector& addParts, const stk::mesh::PartVector& rmParts)
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

    void write_summary(int mod_cycle_count, bool sort=true)
    {
    }
};

class ModificationSummary
{
public:
    ModificationSummary(stk::mesh::BulkData& bulkData);

    ~ModificationSummary();

    // void track_create_ghosting();
    void track_induced_parts(stk::mesh::Entity entity, stk::mesh::Entity e_to, const stk::mesh::PartVector& add_parts, const stk::mesh::PartVector& emptyParts);

    void track_change_ghosting(const stk::mesh::Ghosting & ghosts, const std::vector<stk::mesh::EntityProc> & add_send , const std::vector<stk::mesh::EntityKey> & remove_receive );

    void track_destroy_relation(stk::mesh::Entity e_from, stk::mesh::Entity e_to, stk::mesh::RelationIdentifier rel);

    void track_declare_relation(stk::mesh::Entity e_from, stk::mesh::Entity e_to, stk::mesh::RelationIdentifier rel, stk::mesh::Permutation permut);

    void track_declare_entity(stk::mesh::EntityRank rank, stk::mesh::EntityId newId, const stk::mesh::PartVector& addParts);

    void track_change_entity_owner(const std::vector<stk::mesh::EntityProc> &changes);

    void track_change_entity_id(stk::mesh::EntityId newId, stk::mesh::Entity entity);

    void track_destroy_entity(stk::mesh::Entity entity);

    void track_change_entity_parts(stk::mesh::Entity entity, const stk::mesh::PartVector& addParts, const stk::mesh::PartVector& rmParts);

    void track_comm_map_insert(stk::mesh::Entity entity, const stk::mesh::EntityCommInfo & val);

    void track_comm_map_erase(stk::mesh::EntityKey key, const stk::mesh::EntityCommInfo & val);

    void track_comm_map_erase(stk::mesh::EntityKey key, const stk::mesh::Ghosting & val);

    void track_comm_map_clear_ghosting(stk::mesh::EntityKey key);

    void track_comm_map_clear(stk::mesh::EntityKey key);

    void write_summary(int mod_cycle_count, bool sort=true);

private:

    void clear_summary();

    bool isValid(stk::mesh::Entity entity) const;

    stk::mesh::EntityKey getEntityKey(stk::mesh::Entity entity) const;

    void addEntityKeyAndStringToTracker(stk::mesh::EntityKey key, const std::string& string);

    std::string get_filename(int mod_cycle_count) const;

    int my_proc_id() const;

    void writeParts(std::ostringstream& os, const std::string &label, const stk::mesh::PartVector& parts);

    stk::mesh::BulkData &m_bulkData;
    std::vector<std::pair<stk::mesh::EntityKey, std::string> > m_stringTracker;
    int m_lastModCycle;
    int m_modCounter;
};

} // namespace

#endif
