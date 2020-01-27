#include <stk_mesh/base/SideSetEntry.hpp>
#include "stk_mesh/base/BulkData.hpp"
#include "stk_util/util/ReportHandler.hpp"

namespace stk
{
namespace mesh
{
    SideSet::SideSet(const BulkData& bulk, bool fromInput)
    : m_bulk(bulk), m_fromInput(fromInput), m_part(nullptr)
    {

    }

    SideSet::SideSet(const BulkData& bulk, const std::vector<SideSetEntry>& data, bool fromInput)
    : m_bulk(bulk), m_fromInput(fromInput), m_data(data), m_part(nullptr)
    {

    }

    bool SideSet::is_from_input() const
    {
    	return m_fromInput;
    }

    void SideSet::add(const SideSetEntry& entry)
    {
    	ThrowRequireMsg(m_bulk.entity_rank(entry.element) == stk::topology::ELEMENT_RANK,
                       "ERROR, stk::mesh::SideSet::add only allows element-rank entities.");
    	stk::util::insert_keep_sorted_and_unique(entry, m_data);
    }

    void SideSet::add(stk::mesh::Entity element, stk::mesh::ConnectivityOrdinal side)
    {
        add(SideSetEntry{element, side});
    }

    void SideSet::add(const std::vector<SideSetEntry>& entries)
    {
        m_data.insert(m_data.end(), entries.begin(), entries.end());
        stk::util::sort_and_unique(m_data);
    }

    bool SideSet::contains(const SideSetEntry& entry) const
    {
        std::vector<SideSetEntry>::const_iterator beginIter = begin();
        std::vector<SideSetEntry>::const_iterator endIter = end();

        std::vector<SideSetEntry>::const_iterator lowerBound = std::lower_bound(beginIter, endIter, entry);
        std::vector<SideSetEntry>::const_iterator upperBound = std::upper_bound(beginIter, endIter, entry);
        return (lowerBound != upperBound && lowerBound != endIter);
    }

    bool SideSet::contains(stk::mesh::Entity elem, stk::mesh::ConnectivityOrdinal side) const
    {
    	return contains(SideSetEntry{elem, side});
    }

    SideSetEntry SideSet::operator[](unsigned index) const
    {
        return m_data[index];
    }

    SideSet& SideSet::operator=(const SideSet &rhs)
    {
        m_fromInput = rhs.m_fromInput;
        m_data = rhs.m_data;

        return *this;
    }

    std::vector<SideSetEntry>::iterator SideSet::erase(std::vector<SideSetEntry>::iterator iter)
    {
    	return m_data.erase(iter);
    }

    std::vector<SideSetEntry>::iterator SideSet::erase(std::vector<SideSetEntry>::iterator begin, std::vector<SideSetEntry>::iterator end)
    {
      return m_data.erase(begin, end);
    }

    std::vector<SideSetEntry>::iterator SideSet::begin()
    {
    	return m_data.begin();
    }

    std::vector<SideSetEntry>::iterator SideSet::end()
    {
    	return m_data.end();
    }

    std::vector<SideSetEntry>::const_iterator SideSet::begin() const
    {
    	return m_data.begin();
    }

    std::vector<SideSetEntry>::const_iterator SideSet::end() const
    {
    	return m_data.end();
    }

    void SideSet::clear()
    {
    	m_data.clear();
    }

    size_t SideSet::size() const
    {
    	return m_data.size();
    }

    void SideSet::resize(size_t n)
    {
    	m_data.resize(n);
    }

    const std::string& SideSet::get_name() const
    {
    	return m_name;
    }

    void SideSet::set_name(const std::string& name)
    {
    	m_name = name;
    }

    void SideSet::set_part(const stk::mesh::Part* part)
    {
        m_part = part;
        if(nullptr != part) {
          set_name(part->name());
        }
    }

    const stk::mesh::Part* SideSet::get_part() const
    {
      return m_part;
    }

    void remove_element_entries_from_sidesets(BulkData& mesh, const Entity entity, std::set<const stk::mesh::Part*> *touchedSidesetParts)
    {
      if (mesh.entity_rank(entity) == stk::topology::ELEMENT_RANK && mesh.num_sides(entity) > 0)
      {
        std::vector<SideSet* > sidesets = mesh.get_sidesets();
        for (stk::mesh::SideSet* sideset : sidesets)
        {
          std::vector<SideSetEntry>::iterator lowerBound = std::lower_bound(sideset->begin(), sideset->end(), SideSetEntry(entity, 0));
          std::vector<SideSetEntry>::iterator upperBound = std::upper_bound(sideset->begin(), sideset->end(), SideSetEntry(entity, INVALID_CONNECTIVITY_ORDINAL));
          sideset->erase(lowerBound, upperBound);

          if(nullptr != touchedSidesetParts) {
            const stk::mesh::Part* part = sideset->get_part();

            if(nullptr != part) {
              touchedSidesetParts->insert(part);
            }
          }
        }
      }
    }
}
}
