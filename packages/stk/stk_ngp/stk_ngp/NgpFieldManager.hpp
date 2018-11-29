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

#ifndef NGPFIELDMANAGER_H
#define NGPFIELDMANAGER_H

#include "stk_ngp/Ngp.hpp"
#include "stk_ngp/NgpField.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include <vector>

namespace ngp {

class FieldManager
{
public:
  FieldManager(stk::mesh::BulkData & bulk)
    : m_bulk(bulk)
  {
    m_fields.resize(m_bulk.mesh_meta_data().get_fields().size(), nullptr);
  }

  ~FieldManager() {
    for ( ngp::FieldBase * field : m_fields) {
      delete field;
    }
  }

  void add_fields( const std::vector<stk::mesh::FieldBase*> & fields ) {
    if (m_bulk.mesh_meta_data().get_fields().size() != m_fields.size()) {
      m_fields.resize(m_bulk.mesh_meta_data().get_fields().size());
    }

    for (stk::mesh::FieldBase * field : fields) {
      const unsigned fieldOrdinal = field->mesh_meta_data_ordinal();

      if (field->type_is<double>()) {
        m_fields[fieldOrdinal] = new ngp::Field<double>(m_bulk, *field);
      }
      else if (field->type_is<int>()) {
        m_fields[fieldOrdinal] = new ngp::Field<int>(m_bulk, *field);
      }
      else if (field->type_is<unsigned int>()) {
        m_fields[fieldOrdinal] = new ngp::Field<unsigned int>(m_bulk, *field);
      }
      else if (field->type_is<long>()) {
        m_fields[fieldOrdinal] = new ngp::Field<long>(m_bulk, *field);
      }
      else if (field->type_is<unsigned long>()) {
        m_fields[fieldOrdinal] = new ngp::Field<unsigned long>(m_bulk, *field);
      }
      else if (field->type_is<long long>()) {
        m_fields[fieldOrdinal] = new ngp::Field<long long>(m_bulk, *field);
      }
      else if (field->type_is<unsigned long long>()) {
        m_fields[fieldOrdinal] = new ngp::Field<unsigned long long>(m_bulk, *field);
      }
      else if (field->type_is<short>()) {
        m_fields[fieldOrdinal] = new ngp::Field<short>(m_bulk, *field);
      }
      else if (field->type_is<unsigned short>()) {
        m_fields[fieldOrdinal] = new ngp::Field<unsigned short>(m_bulk, *field);
      }
      else if (field->type_is<float>()) {
        m_fields[fieldOrdinal] = new ngp::Field<float>(m_bulk, *field);
      }
      else {
        ThrowErrorMsg("Unsupported Field type: '" << field->data_traits().name << "'");
      }
    }
  }

  template <typename T>
  ngp::Field<T> & get_field(unsigned fieldOrdinal) const {
    return *dynamic_cast<ngp::Field<T>*>(m_fields[fieldOrdinal]);
  }


private:
  stk::mesh::BulkData & m_bulk;
  std::vector<ngp::FieldBase*> m_fields;

};

}


#endif
