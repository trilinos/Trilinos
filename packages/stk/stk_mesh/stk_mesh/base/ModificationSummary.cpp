#include "stk_mesh/base/ModificationSummary.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/baseImpl/MeshCommVerify.hpp>
#include <iomanip>

namespace stk
{

static int modificationSummaryNumber = 0;

ModificationSummary::ModificationSummary(stk::mesh::BulkData& bulkData)
 : m_bulkData(bulkData),
   m_stringTracker(),
   m_lastModCycle(-1),
   m_setModCycleCounter(false),
   m_modCounter(0)
{
    m_modificationSummaryNumber = modificationSummaryNumber++;
}

ModificationSummary::~ModificationSummary()
{
}

void ModificationSummary::track_induced_parts(stk::mesh::Entity e_from, stk::mesh::Entity e_to, const stk::mesh::OrdinalVector& /*add_parts*/, const stk::mesh::OrdinalVector& /*emptyParts*/)
{
    if ( isValid(e_from) && isValid(e_to) )
    {
        std::ostringstream os;
        os << "Induced parts to " << getEntityKey(e_to) << " from " << getEntityKey(e_from) << std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(e_from), os.str());
        addEntityKeyAndStringToTracker(getEntityKey(e_to), os.str());
    }
}

void ModificationSummary::track_change_ghosting(const stk::mesh::Ghosting & ghosts, const std::vector<stk::mesh::EntityProc> & add_send,
                                                const std::vector<stk::mesh::Entity> & remove_receive)
{
    std::ostringstream os;

    for(size_t i=0;i<add_send.size();++i)
    {
        os << "change_ghosting add-send " << getEntityKey(add_send[i].first) << " to P" << add_send[i].second << " for ghosting " << ghosts.name() << std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(add_send[i].first), os.str());
        os.str("");
    }

    for (size_t i=0;i<remove_receive.size();++i)
    {
        os << "change_ghosting rm-recv " << getEntityKey(remove_receive[i]) << " for ghosting " << ghosts.name() << std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(remove_receive[i]), os.str());
        os.str("");
    }
}

void ModificationSummary::track_change_ghosting(const stk::mesh::Ghosting & ghosts, const std::vector<stk::mesh::EntityProc> & add_send,
                                                const std::vector<stk::mesh::EntityProc> & remove_receive)
{
    std::ostringstream os;

    for(size_t i=0;i<add_send.size();++i)
    {
        os << "change_ghosting add-send " << getEntityKey(add_send[i].first) << " to P" << add_send[i].second << " for ghosting " << ghosts.name() << std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(add_send[i].first), os.str());
        os.str("");
    }

    for (size_t i=0;i<remove_receive.size();++i)
    {
        os << "change_ghosting rm-recv " << getEntityKey(remove_receive[i].first) << " from P"<<remove_receive[i].second<<" for ghosting " << ghosts.name()<<std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(remove_receive[i].first), os.str());
        os.str("");
    }
}

void ModificationSummary::track_add_to_ghosting(const stk::mesh::Ghosting & ghosts, const std::vector<stk::mesh::EntityProc> & add_send )
{
    std::ostringstream os;

    for(size_t i=0;i<add_send.size();++i)
    {
        os << "add ghost " << getEntityKey(add_send[i].first) << " to P" << add_send[i].second << " for ghosting " << ghosts.name() << std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(add_send[i].first), os.str());
        os.str("");
    }
}

void ModificationSummary::track_destroy_relation(stk::mesh::Entity e_from, stk::mesh::Entity e_to, stk::mesh::RelationIdentifier rel)
{
    if(isValid(e_from) && isValid(e_to))
    {
        std::ostringstream os;
        os << "destroy_relation " << getEntityKey(e_from) <<mesh::impl::get_sharing_str(m_bulkData,e_from)<< " to " << getEntityKey(e_to) <<mesh::impl::get_sharing_str(m_bulkData,e_to)<< " rel-id: " << rel << std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(e_from), os.str());
        addEntityKeyAndStringToTracker(getEntityKey(e_to), os.str());
    }
}

void ModificationSummary::track_declare_relation(stk::mesh::Entity e_from, stk::mesh::Entity e_to,
                                                 stk::mesh::RelationIdentifier rel, stk::mesh::Permutation permut)
{
    if(isValid(e_from) && isValid(e_to))
    {
        std::ostringstream os;
        os << "declare_relation " << getEntityKey(e_from) <<mesh::impl::get_sharing_str(m_bulkData,e_from)<< " to " << getEntityKey(e_to)
           <<mesh::impl::get_sharing_str(m_bulkData,e_to)<< " rel-id: " << rel << " perm: " << static_cast<unsigned>(permut) << std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(e_from), os.str());
        addEntityKeyAndStringToTracker(getEntityKey(e_to), os.str());
    }
}

void ModificationSummary::track_declare_entity(stk::mesh::EntityRank rank, stk::mesh::EntityId newId, const stk::mesh::PartVector& addParts)
{
    std::ostringstream os;
    stk::mesh::EntityKey key(rank, newId);
    os << "declare_entity " << key << " on parts: ";
    stk::mesh::OrdinalVector partOrdinals;
    partOrdinals.reserve(addParts.size());
    for(const stk::mesh::Part* part : addParts) {
        partOrdinals.push_back(part->mesh_meta_data_ordinal());
    }
    writeParts(os, "adding parts:", partOrdinals);
    os << std::endl;
    addEntityKeyAndStringToTracker(key, os.str());
}

void ModificationSummary::track_set_global_id(stk::mesh::Entity entity, uint32_t newId)
{
    if (isValid(entity)) {
        stk::mesh::EntityKey oldKey = getEntityKey(entity);
        std::ostringstream os;
        os << "Changing Fmwk global id " << oldKey << " to " << "(" << oldKey.rank() << "," << newId << ")" << std::endl;
        addEntityKeyAndStringToTracker(oldKey, os.str());
        stk::mesh::EntityKey newKey(oldKey.rank(),newId);
        addEntityKeyAndStringToTracker(newKey, os.str());
    }
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
                os << "change_entity_owner " << getEntityKey(entity) << " from P" << my_proc_id() << " to P" << changes[i].second << std::endl;
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
        os << "change_entity_id " << getEntityKey(entity) << " to " << newId << std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(entity), os.str());
        if (newId < stk::mesh::EntityKey::MAX_ID && getEntityKey(entity).id() != newId) {
            stk::mesh::EntityKey newKey(getEntityKey(entity).rank(),newId);
            addEntityKeyAndStringToTracker(newKey, os.str());
        }
    }
}

void ModificationSummary::track_destroy_entity(stk::mesh::Entity entity)
{
    if(isValid(entity))
    {
        std::ostringstream os;
        stk::mesh::EntityKey key = getEntityKey(entity);
        os << "destroy_entity " << key << std::endl;
        addEntityKeyAndStringToTracker(key, os.str());
    }
}

void ModificationSummary::track_change_entity_parts(stk::mesh::Entity entity, const stk::mesh::OrdinalVector& addParts, const stk::mesh::OrdinalVector& rmParts)
{
    if(isValid(entity))
    {
        std::ostringstream os;
        os << "change_entiy_parts " << getEntityKey(entity) << ":";
        writeParts(os, "add parts:", addParts);
        writeParts(os, "rm parts:", rmParts);
        os << std::endl;
        addEntityKeyAndStringToTracker(getEntityKey(entity), os.str());
    }
}

void ModificationSummary::track_comm_map_insert(stk::mesh::Entity entity, const stk::mesh::EntityCommInfo & val)
{
    if(isValid(entity))
    {
        std::ostringstream os;
        os << "Adding " << getEntityKey(entity) << " to comm_map for ghst-id: " << val.ghost_id << " to P" << val.proc << "\n";
        addEntityKeyAndStringToTracker(getEntityKey(entity), os.str());
    }
}

void ModificationSummary::track_comm_map_erase(stk::mesh::EntityKey key, const stk::mesh::EntityCommInfo & val)
{
    if(key != stk::mesh::EntityKey())
    {
        std::ostringstream os;
        os << "Erasing " << key << " from comm_map for ghst-id: " << val.ghost_id << " to P" << val.proc << "\n";
        addEntityKeyAndStringToTracker(key, os.str());
    }
}

void ModificationSummary::track_comm_map_erase(stk::mesh::EntityKey key, const stk::mesh::Ghosting & val)
{
    if(key != stk::mesh::EntityKey())
    {
        std::ostringstream os;
        os << "Erasing " << key << " from comm_map for ghst-id: " << val.ordinal() << " for all procs\n";
        addEntityKeyAndStringToTracker(key, os.str());
    }
}

void ModificationSummary::track_comm_map_clear_ghosting(stk::mesh::EntityKey key)
{
    if(key != stk::mesh::EntityKey())
    {
        std::ostringstream os;
        os << "Erasing " << key << " from all ghosting comm_maps\n";
        addEntityKeyAndStringToTracker(key, os.str());
    }
}

void ModificationSummary::track_comm_map_clear(stk::mesh::EntityKey key)
{
    if(key != stk::mesh::EntityKey())
    {
        std::ostringstream os;
        os << "Erasing " << key << " from all ghosting and sharing comm_maps\n";
        addEntityKeyAndStringToTracker(key, os.str());
    }
}

void ModificationSummary::track_comm_list_remove(stk::mesh::EntityKey key)
{
    if(key != stk::mesh::EntityKey())
    {
        std::ostringstream os;
        os << "Removing " << key << " from comm_list\n";
        addEntityKeyAndStringToTracker(key, os.str());
    }
}

void ModificationSummary::track_set_parallel_owner_rank_but_not_comm_lists(stk::mesh::Entity entity, int old_owner, int new_owner)
{
    stk::mesh::EntityKey key = getEntityKey(entity);
    if(key != stk::mesh::EntityKey())
    {
        std::ostringstream os;
        os << "Changing owner (in bucket data) " << key << " from P" << old_owner << " to P" << new_owner << std::endl;
        addEntityKeyAndStringToTracker(key, os.str());
    }
}

void ModificationSummary::track_change_owner_in_comm_data(stk::mesh::EntityKey key, int old_owner, int new_owner)
{
    if(key != stk::mesh::EntityKey())
    {
        std::ostringstream os;
        os << "Changing owner (in comm data) " << key << " from P" << old_owner << " to P" << new_owner << std::endl;
        addEntityKeyAndStringToTracker(key, os.str());
    }
}

void ModificationSummary::write_summary(int mod_cycle_count, bool sort)
{
    if(mod_cycle_count > m_lastModCycle) {
        std::string filename = get_filename(mod_cycle_count);
        std::ofstream out(filename.c_str());

        if (sort) {
            std::sort(m_stringTracker.begin(), m_stringTracker.end());
        }

        if(m_stringTracker.empty()) {
            out << "*** Nothing happened this cycle on this processor ***\n";
        }
        else {
            for(size_t i = 0; i < m_stringTracker.size(); ++i) {
              std::ostringstream oss;
              oss<<m_stringTracker[i].first;
                out << std::left << std::setw(20) << oss.str() << " " <<m_stringTracker[i].second;
            }
        }

        out.close();

        m_lastModCycle = mod_cycle_count;
        m_setModCycleCounter = false;
    }
    else {
        std::cerr << "*** ERROR ***: Trying to write summary for invalid mod cycle " << mod_cycle_count << std::endl;
    }
    clear_summary();
}

void ModificationSummary::set_sync_count(int syncCount)
{
  if (syncCount < m_lastModCycle || syncCount > (m_lastModCycle + 1)) {
    std::cerr << "*** ERROR ***: set_sync_count, syncCount=" << syncCount << " but last-mod-cycle is " << m_lastModCycle << std::endl;
  }
  m_lastModCycle = syncCount-1;
  m_setModCycleCounter = true;
  m_stringTracker.emplace_back(stk::mesh::EntityKey(stk::topology::NODE_RANK,0), "Reset Mod Cycle, probably for no-op mod-cycle " + std::to_string(m_lastModCycle));
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
    const int numDigitFill = 6;
    std::ostringstream os;
    os << "[" << std::setw(numDigitFill) << std::setfill('0') << m_modCounter << "] " << string;
    m_stringTracker.emplace_back(key, os.str());
    m_modCounter++;
}

int find_how_much_to_pad(int number, int width)
{
    std::ostringstream tempStream;
    tempStream << number;
    return width - tempStream.str().length();
}

std::string string_of_zeros(int number)
{
    std::ostringstream tempStream;
    for (int i=0 ; i<number ; ++i) tempStream << 0;
    return tempStream.str();
}

std::string pad_int_with_zeros(int number, int width)
{
    const int howMuchToPad = find_how_much_to_pad(number,width);
    std::string zerosString = string_of_zeros(howMuchToPad);
    std::ostringstream tempStream;
    tempStream << zerosString << number;
    return tempStream.str();
}

std::string ModificationSummary::get_filename(int mod_cycle_count) const
{
    std::ostringstream os;
    os << "mesh_mod_P" << pad_int_with_zeros(my_proc_id(),3)
            << "_B" << pad_int_with_zeros(m_modificationSummaryNumber,1)
            << "_C" << pad_int_with_zeros(mod_cycle_count,4);
    if (m_setModCycleCounter) {
      os    << ".1";
    }
    os << ".txt";
    return os.str();
}

int ModificationSummary::my_proc_id() const
{
    if (-1 == m_procId)  {
        return m_bulkData.parallel_rank();
    }
    return m_procId;
}

void ModificationSummary::writeParts(std::ostringstream& os, const std::string &label, const stk::mesh::OrdinalVector& parts)
{
    const stk::mesh::PartVector& allParts = m_bulkData.mesh_meta_data().get_parts();
    if(!parts.empty())
    {
        std::vector<std::string> names(parts.size());
        for(size_t i = 0; i < parts.size(); ++i)
        {
            names[i] = allParts[parts[i]]->name();
        }
        std::sort(names.begin(), names.end());

        os << "\t" << label << "  (";
        for(size_t i = 0; i < names.size(); ++i)
        {
            os << " " << names[i];
        }
        os << " )";
    }
}

} // namespace stk
