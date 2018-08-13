// Copyright (c) 2015, Sandia Corporation.
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
  {
    m_entity_keys.reserve(entities.size());
    for (size_t i=0 ; i<entities.size() ; ++i) {
      ThrowRequireMsg(m_mesh.is_valid(entities[i]),
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

  virtual ~TransferCopyByIdStkMeshAdapter() {}

  const double* field_data(const Mesh_ID & id, const unsigned field_index) const
  {
    EntityKey key(static_cast<EntityKey::entity_key_t>(id));
    const mesh::Entity entity = m_mesh.get_entity(key);
    stk::mesh::FieldBase* field=m_transfer_fields[field_index];
    return static_cast<const double*>(stk::mesh::field_data(*field, entity));
  }

  double* field_data(const Mesh_ID & id, const unsigned field_index)
  {
    EntityKey key(static_cast<EntityKey::entity_key_t>(id));
    const mesh::Entity entity = m_mesh.get_entity(key);
    stk::mesh::FieldBase* field=m_transfer_fields[field_index];
    return static_cast<double*>(stk::mesh::field_data(*field, entity));
  }

  stk::mesh::FieldBase* get_field(const unsigned field_index)
  {
    ThrowRequireMsg(field_index < m_transfer_fields.size(),
                    "P" << m_mesh.parallel_rank() <<
                    " stk::transfer::StkMeshAdapter Error, attempt to access invalid field index [" << field_index <<
                    "] in get_field(const unsigned field_index) is invalid!");

    stk::mesh::FieldBase* field=m_transfer_fields[field_index];
    return field;
  }

  unsigned field_data_size(const Mesh_ID & id, const unsigned field_index) const
  {
    EntityKey key(static_cast<EntityKey::entity_key_t>(id));
    const mesh::Entity entity    = m_mesh.get_entity(key);
    const mesh::FieldBase &field = *m_transfer_fields[field_index];

    return stk::mesh::field_scalars_per_entity(field,entity);
  }

  unsigned num_fields() const
  {
    return m_transfer_fields.size();
  }

  ParallelMachine comm() const
  {
    return m_comm;
  }

  const MeshIDVector & get_mesh_ids() const
  {
    return m_ids;
  }

  bool is_locally_owned(const Mesh_ID & k) const
  {
    EntityKey key(static_cast<EntityKey::entity_key_t>(k));
    Entity entity = m_mesh.get_entity(key);
    return m_mesh.is_valid(entity) && m_mesh.bucket(entity).owned();
  }

  std::string print_mesh_id(const Mesh_ID& k) const
  {
    EntityKey key(static_cast<EntityKey::entity_key_t>(k));
    std::ostringstream oss;
    oss << key;
    return oss.str();
  }

  void centroid(const Mesh_ID& k, double coords[3]) const
  {
    for (int i=0 ; i<3 ; ++i) { coords[i] = 0.0; }
    const unsigned dimension = m_mesh.mesh_meta_data().spatial_dimension();
    EntityKey key(static_cast<EntityKey::entity_key_t>(k));
    Entity entity = m_mesh.get_entity(key);

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

private:
  stk::mesh::BulkData & m_mesh;
  stk::ParallelMachine m_comm;
  const stk::mesh::FieldBase* m_coordinates_field;
  FieldVector m_transfer_fields;
  EntityKeyVector m_entity_keys;
  MeshIDVector m_ids;
};

} // transfer
} // stk
#endif //  STK_STKMESHADAPTER_HPP
