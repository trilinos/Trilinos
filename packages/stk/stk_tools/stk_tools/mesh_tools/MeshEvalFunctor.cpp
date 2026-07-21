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

#include <stk_tools/mesh_tools/MeshEvalFunctor.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldBase.hpp>

namespace stk::tools {

MeshEvalFunctor::MeshEvalFunctor(std::shared_ptr<stk::expreval::Eval> eval,
                                 stk::mesh::MetaData& meta,
                                 stk::mesh::EntityRank rank)
 : m_eval(eval)
{
  if (!m_eval->getParseStatus()) {
    m_eval->parse();
  }
  STK_ThrowRequireMsg(m_eval->getParseStatus(),"MeshEvalFunctor: expreval::Eval::getParseStatus == false. Can't continue.");

  for(auto&& varName : m_eval->get_variable_names()) {
    const bool fieldAlreadyExists = meta.get_field(rank, varName) != nullptr;
    if (!fieldAlreadyExists) {
      if (stk::equal_case(varName, "x") ||
          stk::equal_case(varName, "y") ||
          stk::equal_case(varName, "z"))
      {
        auto* coordField = meta.coordinate_field();
        STK_ThrowRequireMsg(coordField != nullptr,"MeshEvalFunctor: coordinate field is null");
        
        m_fieldData.push_back(coordField->data<double,stk::mesh::ReadWrite>());
      }
      else {
        if (meta.is_commit() && !meta.are_late_fields_enabled()) {
          meta.enable_late_fields();
        }
        auto& field = meta.declare_field<double>(rank, varName);
        stk::mesh::put_field_on_mesh(field, meta.universal_part(), m_eval->getVariable(varName).getLength(), nullptr);
        m_fieldData.push_back(field.data<stk::mesh::ReadWrite>());
      }
    }
  }
}

void MeshEvalFunctor::bind_vars(stk::mesh::Entity entity) const
{
  unsigned varIdx = 0;
  for(auto&& varName : m_eval->get_variable_names()) {
    if (stk::equal_case(varName, "x")) {
      m_eval->bindVariable(varName, m_fieldData[varIdx].entity_values(entity)(0_comp), 1);
    }
    else if (stk::equal_case(varName, "y")) {
      m_eval->bindVariable(varName, m_fieldData[varIdx].entity_values(entity)(1_comp), 1);
    }
    else if (stk::equal_case(varName, "z")) {
      m_eval->bindVariable(varName, m_fieldData[varIdx].entity_values(entity)(2_comp), 1);
    }
    else {
      m_eval->bindVariable(varName, m_fieldData[varIdx].entity_values(entity)(0_comp), 1);
    }
    ++varIdx;
  }
}

void MeshEvalFunctor::operator()(const stk::mesh::BulkData& /*mesh*/,
                                 stk::mesh::Entity entity) const
{
  bind_vars(entity);
  m_eval->evaluate();
}

} // namespace stk::tools

