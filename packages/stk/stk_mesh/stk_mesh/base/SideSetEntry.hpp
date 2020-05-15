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

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_util/util/SortAndUnique.hpp>

namespace stk
{
namespace mesh
{

class BulkData;

struct SideSetEntry
{
  SideSetEntry() : element(stk::mesh::Entity()), side(stk::mesh::INVALID_CONNECTIVITY_ORDINAL){};
  SideSetEntry(stk::mesh::Entity in_element, stk::mesh::ConnectivityOrdinal in_side)
    : element(in_element),
      side(in_side)
  {  }
  SideSetEntry(stk::mesh::Entity in_element, int in_side)
    : SideSetEntry(in_element, static_cast<stk::mesh::ConnectivityOrdinal>(in_side))
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

  stk::mesh::Entity element;
  stk::mesh::ConnectivityOrdinal side;
};

class SideSet
{
public:
    typedef std::vector<SideSetEntry>::value_type value_type;

    SideSet(const BulkData& bulk, bool fromInput = false);
    SideSet(const BulkData& bulk, const std::vector<SideSetEntry>& data, bool fromInput = false);

    bool is_from_input() const;
    void add(const SideSetEntry& entry);
    void add(stk::mesh::Entity element, stk::mesh::ConnectivityOrdinal side);
    void add(const std::vector<SideSetEntry>& entries);

    bool contains(const SideSetEntry& entry) const;
    bool contains(stk::mesh::Entity elem, stk::mesh::ConnectivityOrdinal side) const;

    SideSetEntry operator[](unsigned index) const;
    SideSet& operator=(const SideSet &rhs);

    std::vector<SideSetEntry>::iterator erase(std::vector<SideSetEntry>::iterator iter);
    std::vector<SideSetEntry>::iterator erase(std::vector<SideSetEntry>::iterator begin, std::vector<SideSetEntry>::iterator end);

    std::vector<SideSetEntry>::iterator begin();
    std::vector<SideSetEntry>::iterator end();

    std::vector<SideSetEntry>::const_iterator begin() const;
    std::vector<SideSetEntry>::const_iterator end() const;

    void clear();

    size_t size() const;
    void resize(size_t n);

    const std::string& get_name() const;
    void set_name(const std::string& name);
    void set_part(const stk::mesh::Part* part);
    const stk::mesh::Part* get_part() const;
    size_t capacity() const { return m_data.capacity(); }


private:
    const BulkData& m_bulk;
    bool m_fromInput;
    std::vector<SideSetEntry> m_data;
    std::string m_name;
    const stk::mesh::Part * m_part = nullptr;
};

typedef std::vector<SideSet*> SideSetVector;

void remove_element_entries_from_sidesets(BulkData& mesh, const Entity entity, std::set<const stk::mesh::Part*> *touchedSidesetParts = nullptr);

}
}

#endif /* SIDESETENTRY_HPP_ */
