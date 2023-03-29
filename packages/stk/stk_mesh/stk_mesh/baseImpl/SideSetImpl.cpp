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
//

#include <stk_mesh/baseImpl/SideSetImpl.hpp>
#include <stk_mesh/base/Part.hpp>       // for insert
#include <stk_mesh/base/BulkData.hpp>

namespace stk {
namespace mesh {
namespace impl {

//----------------------------------------------------------------------

template<typename KEY>
bool SideSetImpl<KEY>::does_sideset_exist(KEY sideset_key) const
{
        return (m_sideSetData.find(sideset_key) != m_sideSetData.end());
}

template<typename KEY>
bool SideSetImpl<KEY>::does_sideset_exist(const stk::mesh::Part &part) const
{
    return does_sideset_exist(m_keyGenerator.generate_key(part));
}

template<typename KEY>
SideSet& SideSetImpl<KEY>::create_sideset(KEY sideset_key, bool fromInput)
{
  auto iter = m_sideSetData.find(sideset_key);

  if(iter == m_sideSetData.end()) {
    m_sideSetSyncCount = m_bulk.synchronized_count();
    m_sideSetData.emplace(sideset_key, stk::mesh::SideSet(m_bulk, fromInput));

    iter = m_sideSetData.find(sideset_key);
    STK_ThrowRequire(iter != m_sideSetData.end());
  }

  return iter->second;
}

template<typename KEY>
SideSet& SideSetImpl<KEY>::create_sideset(const stk::mesh::Part &part, bool fromInput)
{
    SideSet& ss = create_sideset(m_keyGenerator.generate_key(part), fromInput);
    ss.set_part(&part);
    return ss;
}

template<typename KEY>
const SideSet& SideSetImpl<KEY>::get_sideset(KEY sideset_key) const
{
    auto iter = m_sideSetData.find(sideset_key);
    STK_ThrowRequire(iter != m_sideSetData.end());

    return iter->second;
}

template<typename KEY>
const SideSet& SideSetImpl<KEY>::get_sideset(const stk::mesh::Part &part) const
{
    return get_sideset(m_keyGenerator.generate_key(part));
}

template<typename KEY>
SideSet& SideSetImpl<KEY>::get_sideset(KEY sideset_key)
{
    auto iter = m_sideSetData.find(sideset_key);
    STK_ThrowRequire(iter != m_sideSetData.end());

    return iter->second;
}

template<typename KEY>
SideSet& SideSetImpl<KEY>::get_sideset(const stk::mesh::Part &part)
{
    return get_sideset(m_keyGenerator.generate_key(part));
}

template<typename KEY>
std::vector<KEY> SideSetImpl<KEY>::get_sideset_keys() const
{
    std::vector<KEY> keys;
    keys.reserve(m_sideSetData.size());
    for(auto & keyAndSideSet : m_sideSetData)
        keys.push_back(keyAndSideSet.first);

    return keys;
}

template<typename KEY>
void SideSetImpl<KEY>::clear_sideset(KEY sideset_key)
{
    m_sideSetData.erase(sideset_key);
}

template<typename KEY>
void SideSetImpl<KEY>::clear_sideset(const stk::mesh::Part &part)
{
    clear_sideset(m_keyGenerator.generate_key(part));
}

template<typename KEY>
void SideSetImpl<KEY>::clear_sidesets()
{
    m_sideSetData.clear();
}

template<typename KEY>
size_t SideSetImpl<KEY>::size() const
{
    return m_sideSetData.size();
}

template<typename KEY>
std::vector<SideSet *> SideSetImpl<KEY>::get_sidesets()
{
    std::vector<SideSet *> sidesets;
    sidesets.reserve(m_sideSetData.size());
    for(auto & keyAndSideSet : m_sideSetData)
        sidesets.push_back(&keyAndSideSet.second);

    return sidesets;
}

template<typename KEY>
std::vector<const SideSet *> SideSetImpl<KEY>::get_sidesets() const
{
    std::vector<const SideSet *> sidesets;
    sidesets.reserve(m_sideSetData.size());
    for(auto & keyAndSideSet : m_sideSetData)
        sidesets.push_back(&keyAndSideSet.second);

    return sidesets;
}

template<typename KEY>
bool SideSetImpl<KEY>::was_mesh_modified_since_sideset_creation() const
{
    return m_sideSetSyncCount != m_bulk.synchronized_count();
}

template<typename KEY>
void SideSetImpl<KEY>::set_sideset_sync_count(unsigned syncCount)
{
    m_sideSetSyncCount = syncCount;
}

template class SideSetImpl<unsigned>;
template class SideSetImpl<int64_t>;
template class SideSetImpl<std::string>;

} // namespace impl
} // namespace mesh
} // namespace stk

