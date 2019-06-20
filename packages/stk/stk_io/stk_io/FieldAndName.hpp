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

#ifndef STK_IO_FIELDANDNAME_HPP_
#define STK_IO_FIELDANDNAME_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stddef.h>                      // for size_t
#include <string>                        // for string, operator<, etc
#include <vector>                        // for vector
#include "Ioss_EntityType.h"             // for EntityType
#include "Ioss_GroupingEntity.h"
#include "stk_mesh/base/FieldState.hpp"  // for FieldState
#include "stk_mesh/base/FieldBase.hpp"  // for FieldState
#include "stk_io/OutputVariableParams.hpp"
#include <stk_mesh/base/Types.hpp>       // for EntityId, EntityRank
#include <stk_topology/topology.hpp>     // for topology

namespace stk {
namespace io {
enum DataLocation {MESH = 0, UNIVERSAL_NODEBLOCK, GLOBAL};

struct FieldAndName
{
public:
  FieldAndName(stk::mesh::FieldBase *my_field, const std::string& my_db_name) :
    m_field(my_field),
    m_dbName(my_db_name),
    m_variableType(my_field != nullptr ? my_field->entity_rank() : stk::topology::INVALID_RANK),
    m_useAlias(true),
    m_outputParams(my_db_name),
    m_wasFound(false),
    m_forceNodeblockOutput(false) {}

  FieldAndName(stk::mesh::FieldBase *my_field, const std::string& my_db_name, stk::mesh::EntityRank my_var_type) :
    m_field(my_field),
    m_dbName(my_db_name),
    m_variableType(my_var_type),
    m_useAlias(true),
    m_outputParams(my_db_name),
    m_wasFound(false),
    m_forceNodeblockOutput(false) {}

  stk::mesh::FieldBase *field() const {return m_field;};
  std::string db_name() const {return m_dbName;}
  void set_db_name(const std::string &name) {m_dbName = name;}
  stk::mesh::EntityRank type() const {return m_variableType;}
  void set_use_alias(bool useAlias) { m_useAlias = useAlias; }
  bool get_use_alias() const { return m_useAlias; }

  void set_output_params(const OutputVariableParams& outputParams) {m_outputParams = outputParams;}
  bool has_subset_info() const {return m_outputParams.has_subset_info();}
  void add_subset_entity(const std::string &entity) {m_outputParams.add_subset_entity(entity);}
  bool always_output_node_rank() const {return m_outputParams.always_output_node_rank();}
  bool is_nodeset_variable() const {return m_outputParams.is_nodeset_variable();}
  bool apply_to_entity(Ioss::GroupingEntity *entity) const {
      return m_outputParams.apply_to_entity(entity);
  }
private:
  stk::mesh::FieldBase *m_field;
  std::string m_dbName;
  stk::mesh::EntityRank m_variableType;
  bool m_useAlias;
  OutputVariableParams m_outputParams;
public:
  bool m_wasFound;
  // Field is not defined on UNIVERSAL part, but we still want to output it on the nodeblock.
  // This is done to output, for example, nodal fields that exist on an element block without
  // creating a nodeset for the nodes of the element block.
  mutable bool m_forceNodeblockOutput;
};


struct UserDataAndName
{
public:
    UserDataAndName(const std::vector<std::string>& my_parts, const std::string& my_db_name, DataLocation /*loc*/) :
        m_partNames(my_parts),
        m_dbName(my_db_name)
        {}

    const std::vector<std::string>& get_parts() const {return m_partNames;}
    const std::string& db_name() const {return m_dbName;}
    void set_db_name(const std::string &name) {m_dbName = name;}
private:
    std::vector<std::string> m_partNames;
    std::string  m_dbName;
};

}//namespace io
}//namespace stk


#endif /* STK_IO_FIELDANDNAME_HPP_ */
