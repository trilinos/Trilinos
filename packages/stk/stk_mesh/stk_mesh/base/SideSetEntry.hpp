// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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


#ifndef SIDESETENTRY_HPP_
#define SIDESETENTRY_HPP_

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_util/util/SortAndUnique.hpp>

namespace stk
{
namespace mesh
{

class BulkData;

struct SideSetEntry
{
  SideSetEntry()
    : element(Entity()), side(INVALID_CONNECTIVITY_ORDINAL)
  {  };
  SideSetEntry(Entity in_element)
    : element(in_element), side(INVALID_CONNECTIVITY_ORDINAL)
  {  }
  SideSetEntry(Entity in_element, ConnectivityOrdinal in_side)
    : element(in_element), side(in_side)
  {  }
  SideSetEntry(Entity in_element, int in_side)
    : SideSetEntry(in_element, static_cast<ConnectivityOrdinal>(in_side))
  {  }

  bool operator==(const SideSetEntry &rhs) const
  {
      return ((element == rhs.element) && (side == rhs.side));
  }

  bool operator!=(const SideSetEntry &rhs) const
  {
      return ((element != rhs.element) || (side != rhs.side));
  }

  bool operator<(const SideSetEntry &rhs) const
  {
      if(element < rhs.element)
          return true;
      else if (element == rhs.element && side < rhs.side)
          return true;
      else return false;
  }

  Entity element;
  ConnectivityOrdinal side;
};

class SideSet
{
public:
    typedef std::vector<SideSetEntry>::value_type value_type;

    SideSet(const BulkData& bulk, bool fromInput = false);
    SideSet(const BulkData& bulk, const std::vector<SideSetEntry>& data, bool fromInput = false);

    bool is_from_input() const;
    void set_from_input(bool fromInput = false);

    bool add(const SideSetEntry& entry);
    bool add(Entity element, ConnectivityOrdinal side);
    bool add(const std::vector<SideSetEntry>& entries);

    bool contains(const SideSetEntry& entry) const {return std::binary_search(begin(), end(), entry);}

    bool contains(Entity elem, ConnectivityOrdinal side) const
    {
        return contains(SideSetEntry{elem, side});
    }

    SideSetEntry operator[](unsigned index) const;
    SideSet(const SideSet&) = default;
    SideSet& operator=(const SideSet &rhs);

    std::vector<SideSetEntry>::iterator erase(std::vector<SideSetEntry>::iterator iter);
    std::vector<SideSetEntry>::iterator erase(std::vector<SideSetEntry>::iterator begin, std::vector<SideSetEntry>::iterator end);

    std::vector<SideSetEntry>::iterator begin() {return m_data.begin();}
    std::vector<SideSetEntry>::iterator end() {return m_data.end();}

    std::vector<SideSetEntry>::const_iterator begin() const {return m_data.begin();}
    std::vector<SideSetEntry>::const_iterator end() const {return m_data.end();}

    void clear();

    size_t size() const;
    void resize(size_t n);

    const std::string& get_name() const;
    void set_name(const std::string& name);
    void set_part(const Part* part);
    const Part* get_part() const;
    size_t capacity() const { return m_data.capacity(); }

    void set_accept_all_internal_non_coincident_entries(bool flag);
    bool get_accept_all_internal_non_coincident_entries() const;

    bool is_modified() const {return m_isModified;}
    void clear_modification_flag() {m_isModified = false;}

    bool is_automatically_updated() const {return m_isAutomaticallyUpdated;}
    void set_automatic_update(bool flag) {m_isAutomaticallyUpdated = flag;}

    const BulkData& get_bulk_data() const { return m_bulk; }

private:
    const BulkData& m_bulk;
    bool m_fromInput;
    std::vector<SideSetEntry> m_data;
    std::string m_name;
    const Part * m_part = nullptr;
    bool m_acceptAllInternalNonCoincidentEntries = true;
    bool m_isModified = true;
    bool m_isAutomaticallyUpdated = true;
};


class SideSetSelector
{
public:
  SideSetSelector(const Part& part, SideSet* sideset, const Selector* selector)
  : m_part(&part)
  , m_sideset(sideset)
  , m_selector(selector)
  { }

  SideSetSelector(const SideSetSelector&) = default;

  inline bool operator<(const stk::mesh::SideSetSelector& rhs) const
  {
    if(m_sideset->get_part()->mesh_meta_data_ordinal() < rhs.m_sideset->get_part()->mesh_meta_data_ordinal()) {
      return true;
    } else if(m_sideset->get_part()->mesh_meta_data_ordinal() == rhs.m_sideset->get_part()->mesh_meta_data_ordinal() &&
              m_part->mesh_meta_data_ordinal()                <  rhs.m_part->mesh_meta_data_ordinal()) {
      return true;
    } else {
      return false;
    }
  }

  inline bool operator!=( const stk::mesh::SideSetSelector& rhs ) const
  {
    return (m_sideset->get_part()->mesh_meta_data_ordinal() != rhs.m_sideset->get_part()->mesh_meta_data_ordinal()) ||
           (m_part->mesh_meta_data_ordinal()                != rhs.m_part->mesh_meta_data_ordinal());
  }

  inline bool operator==( const stk::mesh::SideSetSelector& rhs ) const
  {
    return (m_sideset->get_part()->mesh_meta_data_ordinal() == rhs.m_sideset->get_part()->mesh_meta_data_ordinal()) &&
           (m_part->mesh_meta_data_ordinal()                == rhs.m_part->mesh_meta_data_ordinal());
  }

  inline const Part&     part()     const {return *m_part;}
  inline const SideSet*  sideset()  const {return  m_sideset;}
  inline       SideSet*  sideset()        {return  m_sideset;}
  inline const Selector* selector() const {return  m_selector;}

private:
  const Part* m_part = nullptr;
  SideSet* m_sideset = nullptr;
  const Selector* m_selector = nullptr;

  SideSetSelector();
};

using SideSetVector = std::vector<SideSet*>;
using SideSetSelectorVector = std::vector<SideSetSelector>;

}
}

#endif /* SIDESETENTRY_HPP_ */
