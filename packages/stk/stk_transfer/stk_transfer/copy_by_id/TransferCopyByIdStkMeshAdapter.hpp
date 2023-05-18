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

#ifndef  STK_STKMESHADAPTER_HPP
#define  STK_STKMESHADAPTER_HPP

#include "TransferCopyByIdMeshAdapter.hpp"

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>


namespace stk {
namespace transfer {


class TransferCopyByIdStkMeshAdapter : public TransferCopyByIdMeshAdapter {
public :
  typedef stk::mesh::Entity                  Entity;
  typedef stk::mesh::EntityKey               EntityKey;
  typedef stk::mesh::EntityVector            EntityVector;
  typedef std::vector<stk::mesh::FieldBase*> FieldVector;
  typedef std::vector<EntityKey>             EntityKeyVector;

  TransferCopyByIdStkMeshAdapter(
      stk::mesh::BulkData &   mesh,
      const EntityVector &    entities,
      const FieldVector &     fields)
  : TransferCopyByIdStkMeshAdapter(mesh, entities, fields, mesh.parallel()) {}
  TransferCopyByIdStkMeshAdapter(
      stk::mesh::BulkData &   mesh,
      const EntityVector &    entities,
      const FieldVector &     fields,
      stk::ParallelMachine global_comm)
    :m_mesh              (mesh)
    ,m_comm(global_comm)
    ,m_coordinates_field (m_mesh.mesh_meta_data().coordinate_field())
    ,m_transfer_fields   (fields)
    ,m_lastKey(stk::mesh::EntityKey::INVALID)
    ,m_lastEntity()
  {
    m_entity_keys.reserve(entities.size());
    for (size_t i=0 ; i<entities.size() ; ++i) {
      STK_ThrowRequireMsg(m_mesh.is_valid(entities[i]),
                      "P" << m_mesh.parallel_rank() <<
                      " stk::transfer::StkMeshAdapter Error, entities[" << i <<
                      "]=" << entities[i] << " in constructor is invalid!");
      const mesh::EntityKey key = m_mesh.entity_key(entities[i]);
      m_entity_keys.push_back(key);
    }
    m_ids.resize(m_entity_keys.size());
    for (size_t i=0 ;i<m_entity_keys.size() ; ++i) {
      m_ids[i] = m_entity_keys[i].m_value;
    }
  }

  virtual ~TransferCopyByIdStkMeshAdapter() override = default;

  stk::mesh::Entity get_cached_entity(stk::mesh::EntityKey key) const
  {
    if (key != m_lastKey) {
      m_lastKey = key;
      m_lastEntity = m_mesh.get_entity(key);
    }
    return m_lastEntity;
  }

  const void* field_data(const Mesh_ID & id, const unsigned field_index) const override
  {
    EntityKey key(static_cast<EntityKey::entity_key_t>(id));
    const mesh::Entity entity = get_cached_entity(key);
    stk::mesh::FieldBase* field=m_transfer_fields[field_index];
    return reinterpret_cast<const void*>(stk::mesh::field_data(*field, entity));
  }

  void* field_data(const Mesh_ID & id, const unsigned field_index) override
  {
    EntityKey key(static_cast<EntityKey::entity_key_t>(id));
    const mesh::Entity entity = get_cached_entity(key);
    stk::mesh::FieldBase* field=m_transfer_fields[field_index];
    return reinterpret_cast<void*>(stk::mesh::field_data(*field, entity));
  }

  std::string field_name(const unsigned field_index) const override
  {
    STK_ThrowRequireMsg(field_index < m_transfer_fields.size(),
                    "P" << m_mesh.parallel_rank() <<
                    " stk::transfer::StkMeshAdapter Error, attempt to access invalid field index [" << field_index <<
                    "] in get_field_name(const unsigned field_index) is invalid!");

    stk::mesh::FieldBase* field=m_transfer_fields[field_index];
    return field->name();
  }

  stk::mesh::FieldBase* get_field(const unsigned field_index)
  {
    STK_ThrowRequireMsg(field_index < m_transfer_fields.size(),
                    "P" << m_mesh.parallel_rank() <<
                    " stk::transfer::StkMeshAdapter Error, attempt to access invalid field index [" << field_index <<
                    "] in get_field(const unsigned field_index) is invalid!");

    stk::mesh::FieldBase* field=m_transfer_fields[field_index];
    return field;
  }

  unsigned field_data_size(const Mesh_ID & id, const unsigned field_index) const override
  {
    EntityKey key(static_cast<EntityKey::entity_key_t>(id));
    const mesh::Entity entity = get_cached_entity(key);
    const mesh::FieldBase &field = *m_transfer_fields[field_index];

    return stk::mesh::field_bytes_per_entity(field, entity);
  }

  unsigned num_fields() const override
  {
    return m_transfer_fields.size();
  }

  ParallelMachine comm() const override
  {
    return m_comm;
  }

  const MeshIDVector & get_mesh_ids() const override
  {
    return m_ids;
  }

  bool is_locally_owned(const Mesh_ID & k) const override
  {
    EntityKey key(static_cast<EntityKey::entity_key_t>(k));
    const mesh::Entity entity = get_cached_entity(key);
    return m_mesh.is_valid(entity) && m_mesh.bucket(entity).owned();
  }

  std::string print_mesh_id(const Mesh_ID& k) const override
  {
    EntityKey key(static_cast<EntityKey::entity_key_t>(k));
    std::ostringstream oss;
    oss << key;
    return oss.str();
  }

  void centroid(const Mesh_ID& k, double coords[3]) const override
  {
    for (int i=0 ; i<3 ; ++i) { coords[i] = 0.0; }
    const unsigned dimension = m_mesh.mesh_meta_data().spatial_dimension();
    EntityKey key(static_cast<EntityKey::entity_key_t>(k));
    const mesh::Entity entity = get_cached_entity(key);

    if (key.rank() == stk::topology::NODE_RANK) {
      const double* c = static_cast<const double*>(stk::mesh::field_data(*m_coordinates_field, entity));
      for (unsigned j=0; j<dimension; ++j) {
        coords[j] = c[j];
      }
    } else {
      const int num_nodes = m_mesh.num_nodes(entity);
      Entity const * node_it = m_mesh.begin_nodes(entity);
      for (int node_i=0 ; node_i<num_nodes ; ++node_i) {
        Entity node = node_it[node_i];
        const double* coord_field =
          static_cast<const double*>(stk::mesh::field_data(*m_coordinates_field, node));
        for (unsigned i=0 ; i<dimension ; ++i) {
          coords[i] += coord_field[i];
        }
      }
      for (unsigned i=0 ; i<dimension ; ++i) {
        coords[i] = coords[i] / num_nodes;
      }
    }
  }

  DataTypeKey::data_t get_field_type(const unsigned fieldIndex) const override
  { 
    STK_ThrowRequireMsg(fieldIndex < m_transfer_fields.size(),
                    "P" << m_mesh.parallel_rank() <<
                    " stk::transfer::StkMeshAdapter Error, attempt to access invalid field index [" << fieldIndex <<
                    "] in get_field_type(const unsigned fieldIndex) is invalid!");

    stk::mesh::FieldBase* field=m_transfer_fields[fieldIndex];

    DataTypeKey::data_t dataType( DataTypeKey::data_t::INVALID_TYPE );

 
    if(field->type_is<unsigned>()) {
      dataType = DataTypeKey::data_t::UNSIGNED_INTEGER;
    }
    else if(field->type_is<int>()) {
      dataType = DataTypeKey::data_t::INTEGER;
    }
    else if(field->type_is<int64_t>()) {
      dataType = DataTypeKey::data_t::LONG_INTEGER;
    }
    else if(field->type_is<uint64_t>()) {
      dataType = DataTypeKey::data_t::UNSIGNED_INTEGER_64;
    }
    else if(field->type_is<double>()) {
      dataType = DataTypeKey::data_t::DOUBLE;
    }
    else if(field->type_is<long double>()) {
      dataType = DataTypeKey::data_t::LONG_DOUBLE;
    } else {
      STK_ThrowRequireMsg(false, "Unsupported data type");
    }
    return dataType;
  }
   

protected:
  stk::mesh::BulkData & m_mesh;
  stk::ParallelMachine m_comm;
  const stk::mesh::FieldBase* m_coordinates_field;
  FieldVector m_transfer_fields;
  EntityKeyVector m_entity_keys;
  MeshIDVector m_ids;
  mutable stk::mesh::EntityKey m_lastKey;
  mutable stk::mesh::Entity m_lastEntity;
};

} // transfer
} // stk
#endif //  STK_STKMESHADAPTER_HPP
