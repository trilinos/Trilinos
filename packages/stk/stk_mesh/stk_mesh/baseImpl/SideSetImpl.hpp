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


#ifndef stk_mesh_SideSetImpl_hpp
#define stk_mesh_SideSetImpl_hpp

//----------------------------------------------------------------------

#include <stddef.h>                     // for size_t, NULL
#include <stdint.h>                     // for uint16_t
#include <algorithm>                    // for max
#include <functional>                   // for less, equal_to
#include <list>                         // for list
#include <map>                          // for map, map<>::value_compare
#include <set>
#include <string>                       // for char_traits, string
#include <utility>                      // for pair
#include <vector>                       // for vector
#include <unordered_map>
#include <stk_mesh/base/SideSetEntry.hpp>
#include "stk_mesh/base/Part.hpp"       // for Part, remove, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc

namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class BulkData; } }

namespace stk {
namespace mesh {
namespace impl {


template<typename KEY>
struct SideSetKeyGenerator {
    KEY generate_key(const stk::mesh::Part &part) const;
    KEY generate_key(int id) const;
};

template<>
struct SideSetKeyGenerator<int64_t> {
    int64_t generate_key(const stk::mesh::Part &part) const {
        return part.id();
    }

    int64_t generate_key(int id) const {
        return id;
    }
};

template<>
struct SideSetKeyGenerator<unsigned> {
    unsigned generate_key(const stk::mesh::Part &part) const {
        return part.mesh_meta_data_ordinal();
    }

    unsigned generate_key(int id) const {
        return id;
    }
};

template<>
struct SideSetKeyGenerator<std::string> {
    std::string generate_key(const stk::mesh::Part &part) const {
        return part.name();
    }

    std::string generate_key(int id) const {
        return std::to_string(id);
    }
};


template<typename KEY>
class SideSetImpl {
public:
    SideSetImpl(const stk::mesh::BulkData &bulk) : m_bulk(bulk) {}

    bool does_sideset_exist(KEY sideset_key) const;

    bool does_sideset_exist(const stk::mesh::Part &part) const;

    SideSet& create_sideset(KEY sideset_key, bool fromInput);

    SideSet& create_sideset(const stk::mesh::Part &part, bool fromInput);

    const SideSet& get_sideset(KEY sideset_key) const;

    const SideSet& get_sideset(const stk::mesh::Part &part) const;

    SideSet& get_sideset(KEY sideset_key);

    SideSet& get_sideset(const stk::mesh::Part &part);

    std::vector<KEY> get_sideset_keys() const;

    bool was_mesh_modified_since_sideset_creation() const;

    void clear_sideset(KEY sideset_key);

    void clear_sideset(const stk::mesh::Part &part);

    void clear_sidesets();

    void set_sideset_sync_count(unsigned syncCount);

    size_t size() const;

    std::vector<SideSet *> get_sidesets();
    std::vector<const SideSet *> get_sidesets() const;

private:
    SideSetImpl();
    SideSetImpl( const SideSetImpl & );
    SideSetImpl & operator = ( const SideSetImpl & );

    const stk::mesh::BulkData &m_bulk;
    std::map<KEY, SideSet> m_sideSetData;
    size_t m_sideSetSyncCount = 0;
    bool m_hasSideSetData = false;

    SideSetKeyGenerator<KEY> m_keyGenerator;
};

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif
