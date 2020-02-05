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

#ifndef STK_MESH_NGPFIELDMANAGER_HPP
#define STK_MESH_NGPFIELDMANAGER_HPP

#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/NgpField.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include <vector>
#include <memory>

namespace stk {
namespace mesh {

class NgpFieldManager
{
public:
  NgpFieldManager()
    : m_bulk(nullptr)
  {
  }

  NgpFieldManager(const stk::mesh::BulkData & bulk)
    : m_bulk(&bulk)
  {
    m_fields.resize(m_bulk->mesh_meta_data().get_fields().size());
  }

  ~NgpFieldManager() {
    clear_fields();
  }

  NgpFieldManager & operator=(const NgpFieldManager & rhs) {
    if (this != &rhs) {
      clear_fields();
      m_bulk = rhs.m_bulk;
      m_fields = rhs.m_fields;
    }
    return *this;
  }

  NgpFieldManager(const NgpFieldManager & rhs) = default;
  NgpFieldManager(NgpFieldManager && rhs) = default;

  template <typename T>
  stk::mesh::NgpField<T> & get_field(unsigned fieldOrdinal) const {
    ThrowRequireMsg(m_bulk != nullptr, "NgpFieldManager ERROR, m_bulk is null.");
    if (m_bulk->mesh_meta_data().get_fields().size() != m_fields.size()) {
      m_fields.resize(m_bulk->mesh_meta_data().get_fields().size());
    }
    ThrowRequireMsg(fieldOrdinal < m_fields.size(), "Invalid field ordinal (" << fieldOrdinal << ").  Only have " << m_fields.size() << " fields.");

    if (m_fields[fieldOrdinal].get() == nullptr) {
      construct_ngp_field(fieldOrdinal);
    }

    std::shared_ptr<stk::mesh::NgpFieldBase> sharedFieldPtr = m_fields[fieldOrdinal];
    stk::mesh::NgpField<T> * returnField = dynamic_cast<stk::mesh::NgpField<T>*>(sharedFieldPtr.get());
    ThrowRequireMsg(returnField != nullptr, "Calling get_field() with the wrong data type for field '"
                    << m_bulk->mesh_meta_data().get_fields()[fieldOrdinal]->name() << "'");
    return *returnField;
  }

  const stk::mesh::BulkData& get_bulk() const {return *m_bulk;}

  std::vector<std::shared_ptr<stk::mesh::NgpFieldBase> >& get_fields() { return m_fields; }
  const std::vector<std::shared_ptr<stk::mesh::NgpFieldBase> >& get_fields() const { return m_fields; }

private:
  void construct_ngp_field( stk::mesh::Ordinal fieldOrdinal ) const {
    stk::mesh::FieldBase & field = *(m_bulk->mesh_meta_data().get_fields()[fieldOrdinal]);

    if (field.type_is<double>()) {
      m_fields[fieldOrdinal] = std::make_shared< stk::mesh::NgpField<double> >(*m_bulk, field);
    }
    else if (field.type_is<int>()) {
      m_fields[fieldOrdinal] = std::make_shared< stk::mesh::NgpField<int> >(*m_bulk, field);
    }
    else if (field.type_is<unsigned int>()) {
      m_fields[fieldOrdinal] = std::make_shared< stk::mesh::NgpField<unsigned int> >(*m_bulk, field);
    }
    else if (field.type_is<long>()) {
      m_fields[fieldOrdinal] = std::make_shared< stk::mesh::NgpField<long> >(*m_bulk, field);
    }
    else if (field.type_is<unsigned long>()) {
      m_fields[fieldOrdinal] = std::make_shared< stk::mesh::NgpField<unsigned long> >(*m_bulk, field);
    }
    else if (field.type_is<long long>()) {
      m_fields[fieldOrdinal] = std::make_shared< stk::mesh::NgpField<long long> >(*m_bulk, field);
    }
    else if (field.type_is<unsigned long long>()) {
      m_fields[fieldOrdinal] = std::make_shared< stk::mesh::NgpField<unsigned long long> >(*m_bulk, field);
    }
    else if (field.type_is<short>()) {
      m_fields[fieldOrdinal] = std::make_shared< stk::mesh::NgpField<short> >(*m_bulk, field);
    }
    else if (field.type_is<unsigned short>()) {
      m_fields[fieldOrdinal] = std::make_shared<stk::mesh::NgpField<unsigned short> >(*m_bulk, field);
    }
    else if (field.type_is<float>()) {
      m_fields[fieldOrdinal] = std::make_shared<stk::mesh::NgpField<float> >(*m_bulk, field);
    }
    else {
      ThrowErrorMsg("Unsupported Field type: '" << field.data_traits().name << "'");
    }
  }

  void clear_fields() {
    for ( std::shared_ptr<stk::mesh::NgpFieldBase> & field : m_fields) {
      if (nullptr != field.get()) {
        field->sync_to_host(); //cheap/no-op if not needed
      }
      field.reset();
    }
    m_fields.clear();
  }

  const stk::mesh::BulkData * m_bulk;
  mutable std::vector<std::shared_ptr<stk::mesh::NgpFieldBase> > m_fields;

};

}
}


#endif
