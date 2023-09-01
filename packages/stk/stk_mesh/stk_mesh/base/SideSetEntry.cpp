#include <stk_mesh/base/SideSetEntry.hpp>
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemGraphCoincidentElems.hpp"

namespace stk
{
namespace mesh
{
    SideSet::SideSet(const BulkData& bulk, bool fromInput)
    : m_bulk(bulk), m_fromInput(fromInput), m_part(nullptr), m_acceptAllInternalNonCoincidentEntries(true)
    {

    }

    SideSet::SideSet(const BulkData& bulk, const std::vector<SideSetEntry>& data, bool fromInput)
    : m_bulk(bulk), m_fromInput(fromInput), m_data(data), m_part(nullptr), m_acceptAllInternalNonCoincidentEntries(true)
    {

    }

    bool SideSet::is_from_input() const
    {
        return m_fromInput;
    }

    void SideSet::set_from_input(bool fromInput)
    {
      m_fromInput = fromInput;
    }

    bool SideSet::add(const SideSetEntry& entry)
    {
        STK_ThrowRequireMsg(m_bulk.entity_rank(entry.element) == stk::topology::ELEMENT_RANK,
                       "ERROR, stk::mesh::SideSet::add only allows element-rank entities.");

        bool modified = stk::util::insert_keep_sorted_and_unique(entry, m_data);
        m_isModified |= modified;
        return modified;
    }

    bool SideSet::add(Entity element, ConnectivityOrdinal side)
    {
        return add(SideSetEntry{element, side});
    }

    bool SideSet::add(const std::vector<SideSetEntry>& entries)
    {
        size_t oldSize = m_data.size();
        m_data.insert(m_data.end(), entries.begin(), entries.end());
        stk::util::sort_and_unique(m_data);
        bool modified = oldSize != m_data.size();
        m_isModified |= modified;
        return modified;
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
        m_isModified = true;
	return m_data.erase(iter);
    }

    std::vector<SideSetEntry>::iterator SideSet::erase(std::vector<SideSetEntry>::iterator begin, std::vector<SideSetEntry>::iterator end)
    {
      m_isModified = true;
      return m_data.erase(begin, end);
    }

    void SideSet::clear()
    {
      if(!m_data.empty()) {
        m_isModified = true;
        m_data.clear();
      }
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

    void SideSet::set_part(const Part* part)
    {
        m_part = part;
        if(nullptr != part) {
          set_name(part->name());
        }
    }

    const Part* SideSet::get_part() const
    {
      return m_part;
    }

    void SideSet::set_accept_all_internal_non_coincident_entries(bool flag)
    {
      m_acceptAllInternalNonCoincidentEntries = flag;
    }

    bool SideSet::get_accept_all_internal_non_coincident_entries() const
    {
      return m_acceptAllInternalNonCoincidentEntries;
    }
}
}
