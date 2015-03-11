#include "stk_mesh/base/ModificationSummary.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <iomanip>

namespace stk
{

ModificationSummary::ModificationSummary(stk::mesh::BulkData& bulkData) :
m_bulkData(bulkData), m_stringTracker(), m_lastModCycle(-1), m_modCounter(0)
{
}

ModificationSummary::~ModificationSummary()
{
}

void ModificationSummary::track_induced_parts(stk::mesh::Entity e_from, stk::mesh::Entity e_to, const stk::mesh::PartVector& add_parts, const stk::mesh::PartVector& emptyParts)
{
    if ( isValid(e_from) && isValid(e_to) )
    {
        std::ostringstream os;
        os << "Inducing parts to entity key " << getEntityKey(e_to) << " from entity key " << getEntityKey(e_from) << std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(e_from), os.str());
        addEntityKeyAndStringToTracker(getEntityKey(e_to), os.str());
    }
}
// void track_create_ghosting();
void ModificationSummary::track_change_ghosting(const stk::mesh::Ghosting & ghosts, const std::vector<stk::mesh::EntityProc> & add_send , const std::vector<stk::mesh::EntityKey> & remove_receive )
{
    std::ostringstream os;

    for(size_t i=0;i<add_send.size();++i)
    {
        os << "Sending ghost key " << getEntityKey(add_send[i].first) << " to processor " << add_send[i].second << " for ghosting " << ghosts.name() << std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(add_send[i].first), os.str());
        os.str("");
    }

    for (size_t i=0;i<remove_receive.size();++i)
    {
        os << "Deleting receive ghost " << remove_receive[i] << " for ghosting " << ghosts.name() << std::endl;
        addEntityKeyAndStringToTracker(remove_receive[i], os.str());
        os.str("");
    }
}

void ModificationSummary::track_destroy_relation(stk::mesh::Entity e_from, stk::mesh::Entity e_to, stk::mesh::RelationIdentifier rel)
{
    if(isValid(e_from) && isValid(e_to))
    {
        std::ostringstream os;
        os << "Destroying a relation from " << getEntityKey(e_from) << " to " << getEntityKey(e_to) << " with relation identifier: " << rel << std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(e_from), os.str());
        addEntityKeyAndStringToTracker(getEntityKey(e_to), os.str());
    }
}

void ModificationSummary::track_declare_relation(stk::mesh::Entity e_from, stk::mesh::Entity e_to, stk::mesh::RelationIdentifier rel, stk::mesh::Permutation permut)
{
    if(isValid(e_from) && isValid(e_to))
    {
        std::ostringstream os;
        os << "Declaring a relation from " << getEntityKey(e_from) << " to " << getEntityKey(e_to) << " with relation identifier: " << rel << " and permutation: " << permut << std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(e_from), os.str());
        addEntityKeyAndStringToTracker(getEntityKey(e_to), os.str());
    }
}

void ModificationSummary::track_declare_entity(stk::mesh::EntityRank rank, stk::mesh::EntityId newId, const stk::mesh::PartVector& addParts)
{
    std::ostringstream os;
    stk::mesh::EntityKey key(rank, newId);
    os << "Declaring new entity with entity key " << key << " on parts: " << std::endl;
    writeParts(os, "adding parts:", addParts);
    addEntityKeyAndStringToTracker(key, os.str());
}

void ModificationSummary::track_change_entity_owner(const std::vector<stk::mesh::EntityProc> &changes)
{
    std::ostringstream os;
    if(!changes.empty())
    {
        for(size_t i = 0; i < changes.size(); ++i)
        {
            stk::mesh::Entity entity = changes[i].first;
            if(isValid(entity))
            {
                os << "Changing owner of entity key " << getEntityKey(entity) << " from proc " << my_proc_id() << " to " << changes[i].second << std::endl;
                addEntityKeyAndStringToTracker(getEntityKey(entity), os.str());
                os.str("");
            }
        }
    }
}

void ModificationSummary::track_change_entity_id(stk::mesh::EntityId newId, stk::mesh::Entity entity)
{
    if(isValid(entity))
    {
        std::ostringstream os;
        os << "Changing id of entity key " << getEntityKey(entity) << std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(entity), os.str());
    }
}

void ModificationSummary::track_destroy_entity(stk::mesh::Entity entity)
{
    if(isValid(entity))
    {
        std::ostringstream os;
        os << "Destroying entity with key " << getEntityKey(entity) << std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(entity), os.str());
    }
}

void ModificationSummary::track_change_entity_parts(stk::mesh::Entity entity, const stk::mesh::PartVector& addParts, const stk::mesh::PartVector& rmParts)
{
    if(isValid(entity))
    {
        std::ostringstream os;
        os << "Part change for entity_key " << getEntityKey(entity) << ":\n";

        writeParts(os, "adding parts:", addParts);
        writeParts(os, "removing parts:", rmParts);

        addEntityKeyAndStringToTracker(getEntityKey(entity), os.str());
    }
}

void ModificationSummary::track_comm_map_insert(stk::mesh::Entity entity, const stk::mesh::EntityCommInfo & val)
{
    if(isValid(entity))
    {
        std::ostringstream os;
        os << "Adding entity with key " << getEntityKey(entity) << " to comm_map for ghosting id: " << val.ghost_id << " to proc " << val.proc << "\n";
        addEntityKeyAndStringToTracker(getEntityKey(entity), os.str());
    }
}

void ModificationSummary::track_comm_map_erase(stk::mesh::EntityKey key, const stk::mesh::EntityCommInfo & val)
{
    if(key != stk::mesh::EntityKey())
    {
        std::ostringstream os;
        os << "Erasing entity with key " << key << " from comm_map for ghosting id: " << val.ghost_id << " to proc " << val.proc << "\n";
        addEntityKeyAndStringToTracker(key, os.str());
    }
}

void ModificationSummary::track_comm_map_erase(stk::mesh::EntityKey key, const stk::mesh::Ghosting & val)
{
    if(key != stk::mesh::EntityKey())
    {
        std::ostringstream os;
        os << "Erasing entity with key " << key << " from comm_map for ghosting id: " << val.ordinal() << " for all procs\n";
        addEntityKeyAndStringToTracker(key, os.str());
    }
}

void ModificationSummary::track_comm_map_clear_ghosting(stk::mesh::EntityKey key)
{
    if(key != stk::mesh::EntityKey())
    {
        std::ostringstream os;
        os << "Erasing entity with key " << key << " from all ghosting comm_maps\n";
        addEntityKeyAndStringToTracker(key, os.str());
    }
}

void ModificationSummary::track_comm_map_clear(stk::mesh::EntityKey key)
{
    if(key != stk::mesh::EntityKey())
    {
        std::ostringstream os;
        os << "Erasing entity with key " << key << " from all ghosting and sharing comm_maps\n";
        addEntityKeyAndStringToTracker(key, os.str());
    }
}

void ModificationSummary::write_summary(int mod_cycle_count, bool sort)
{
    if(mod_cycle_count > m_lastModCycle)
    {
        std::string filename = get_filename(mod_cycle_count);
        std::ofstream out(filename.c_str());

        if (sort)
        {
            std::sort(m_stringTracker.begin(), m_stringTracker.end());
        }

        if(m_stringTracker.empty())
        {
            out << "*** Nothing happened this cycle on this processor ***\n";
        }
        else
        {
            for(size_t i = 0; i < m_stringTracker.size(); ++i)
            {
                out << m_stringTracker[i].first << "\t" << m_stringTracker[i].second;
            }
        }

        out.close();

        m_lastModCycle = mod_cycle_count;
    }
    else
    {
        std::cerr << "*** ERROR ***: Trying to write summary for invalid mod cycle " << mod_cycle_count << std::endl;
    }
    clear_summary();
}

void ModificationSummary::clear_summary()
{
    m_stringTracker.clear();
    m_modCounter = 0;
}

bool ModificationSummary::isValid(stk::mesh::Entity entity) const
{
    return m_bulkData.is_valid(entity);
}

stk::mesh::EntityKey ModificationSummary::getEntityKey(stk::mesh::Entity entity) const
        {
    return m_bulkData.entity_key(entity);
}

void ModificationSummary::addEntityKeyAndStringToTracker(stk::mesh::EntityKey key, const std::string& string)
{
    const int numDigitFill = 8;
    std::ostringstream os;
    os << "[" << std::setw(numDigitFill) << std::setfill('0') << m_modCounter << "] " << string;
    m_stringTracker.push_back(std::make_pair(key, os.str()));
    m_modCounter++;
}

std::string ModificationSummary::get_filename(int mod_cycle_count) const
        {
    std::ostringstream os;
    os << "modification_cycle_" << mod_cycle_count << "_P" << my_proc_id() << ".txt";
    return os.str();
}

int ModificationSummary::my_proc_id() const
{
    return m_bulkData.parallel_rank();
}

void ModificationSummary::writeParts(std::ostringstream& os, const std::string &label, const stk::mesh::PartVector& parts)
{
    if(!parts.empty())
    {
        std::vector<std::string> names(parts.size());
        for(size_t i = 0; i < parts.size(); ++i)
        {
            names[i] = parts[i]->name();
        }
        std::sort(names.begin(), names.end());

        os << "\t" << label << "\n";
        for(size_t i = 0; i < names.size(); ++i)
        {
            os << "\t\t" << names[i] << std::endl;
        }
    }
}

} // namespace stk
