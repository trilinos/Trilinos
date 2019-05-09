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


#ifndef SIDESETENTRY_HPP_
#define SIDESETENTRY_HPP_

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace stk
{
namespace mesh
{

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

    SideSet(bool fromInput = false) : m_fromInput(fromInput) {}
    SideSet(const std::vector<SideSetEntry>& data, bool fromInput = false) : m_fromInput(fromInput), m_data(data)  {}

    bool is_from_input() const {return m_fromInput;}
    void add(const SideSetEntry& entry) { m_data.push_back(entry); }
    void add(stk::mesh::Entity element, stk::mesh::ConnectivityOrdinal side)
    {
        add(SideSetEntry{element, side});
    }

    SideSetEntry operator[](unsigned index) const
    {
        return m_data[index];
    }

    SideSet& operator=(const SideSet &rhs)
    {
        m_fromInput = rhs.m_fromInput;
        m_data = rhs.m_data;

        return *this;
    }

    inline std::vector<SideSetEntry>::iterator erase(std::vector<SideSetEntry>::iterator iter) { return m_data.erase(iter); }

    inline std::vector<SideSetEntry>::iterator begin() { return m_data.begin(); }
    inline std::vector<SideSetEntry>::iterator end() { return m_data.end(); }

    inline std::vector<SideSetEntry>::const_iterator begin() const { return m_data.begin(); }
    inline std::vector<SideSetEntry>::const_iterator end() const { return m_data.end(); }

    void clear() { m_data.clear();}

    size_t size() const { return m_data.size(); }
    void resize(size_t n) { m_data.resize(n); }

    const std::string& get_name() const { return m_name; }
    void set_name(const std::string& name) { m_name = name; }

private:
    bool m_fromInput;
    std::vector<SideSetEntry> m_data;
    std::string m_name;
};

//typedef std::vector<SideSetEntry> SideSet;
typedef std::vector<SideSet*> SideSetVector;

}
}

#endif /* SIDESETENTRY_HPP_ */
