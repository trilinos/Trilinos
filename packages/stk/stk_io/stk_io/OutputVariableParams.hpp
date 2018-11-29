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
//

#ifndef STK_IO_OUTPUTVARIABLEPARAMS_HPP_
#define STK_IO_OUTPUTVARIABLEPARAMS_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <string>                        // for string, operator<, etc
#include <vector>                        // for vector
#include "Ioss_GroupingEntity.h"

namespace stk {
namespace io {

class OutputVariableParams {
public:
    OutputVariableParams(const std::string &varName)
    : m_name(varName),
      m_subsetInclude(false),
      m_nodesetVariable(false),
      m_alwaysOutputNodeRank(false)
    {}

    void set_subset_info(bool isInclude, const std::vector<std::string>& entities)
    {
        m_subsetInclude = isInclude;
        m_entities = entities;
    }

    void add_subset_entity(const std::string &entity)
    {
        if(std::find(m_entities.begin(), m_entities.end(), entity) == m_entities.end()) {
            m_entities.push_back(entity);
        }
    }

    bool apply_to_entity(Ioss::GroupingEntity *entity) const
    {
        if (m_entities.empty())
            return true;

        // NOTE: Need to handle aliases here.  For now brute force the check.
        bool in_list = false;
        std::vector<std::string>::const_iterator I = m_entities.begin();
        while (I != m_entities.end() && !in_list) {
            if (entity->is_alias(*I)) {
                in_list = true;
            }
            ++I;
        }

        if (m_subsetInclude) {
            // List specifies the entities that are to be included...
            return in_list;
        }

        // List specifies the entities that are to be excluded...
        return !in_list;
    }

    bool is_nodeset_variable() const {return m_nodesetVariable;}
    void is_nodeset_variable(bool flag) {m_nodesetVariable = flag;}

    bool always_output_node_rank() const {return m_alwaysOutputNodeRank;}
    void always_output_node_rank(bool flag) {m_alwaysOutputNodeRank = flag;}

    const std::vector<std::string>& get_subset_entities() const {return m_entities;}
    bool has_subset_info() const {return (m_entities.size() > 0);}

    const std::string& name() const {return m_name;}

private:
    std::string m_name;
    bool m_subsetInclude;
    std::vector<std::string> m_entities;
    bool m_nodesetVariable;
    bool m_alwaysOutputNodeRank;
};

}//namespace io
}//namespace stk

#endif /* STK_IO_OUTPUTVARIABLEPARAMS_HPP_ */
