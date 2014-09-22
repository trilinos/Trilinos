// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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


#ifndef  STK_STKNODE_HPP
#define  STK_STKNODE_HPP

#include <boost/shared_ptr.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#include <stk_util/environment/ReportHandler.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>


namespace stk {
namespace transfer {


class STKNode {
public :
  typedef mesh:: Entity                           Entity;
  typedef std::vector<Entity>                     EntityVec;
  typedef mesh:: EntityKey                        EntityKey;
  typedef std::set<EntityKey>                     EntityKeySet;
  typedef search::IdentProc<EntityKey, unsigned>  EntityProc;
  typedef std::vector<EntityProc>                 EntityProcVec;

  typedef search::Point<float>  Point;
  typedef search::Sphere<float> Sphere;

  typedef std::pair<Sphere,EntityProc> BoundingBox;


  STKNode(
      const EntityVec                    &entities,
      const mesh::FieldBase              &coord,
      const std::vector<mesh::FieldBase*> &val,
      const double                        initial_radius) :
    m_bulk_data         (coord.get_mesh()),
    m_mesh_modified     (false),
    m_comm              (m_bulk_data.parallel()),
    m_sphere_rad        (initial_radius),
    m_entity_keys       (entity_keys(m_bulk_data, entities)),
    m_coordinates_field (coord),
    m_values_field      (val),
    m_entities_currently_ghosted()
  {
    const std::string name = "Transfer Ghosting";
    m_bulk_data.modification_begin();
    m_transfer_entity_ghosting = &m_bulk_data.create_ghosting(name);
    m_bulk_data.modification_end();
  }

  // Needed for STK Transfer
  ParallelMachine comm() const {return m_comm;}

  void bounding_boxes (std::vector< std::pair<Sphere,EntityProc> > &v) const
  {
    const unsigned dimension = m_bulk_data.mesh_meta_data().spatial_dimension();
    const float r = m_sphere_rad;
    const int proc_id = parallel_machine_rank(m_comm);

    v.clear();

    for (EntityKeySet::const_iterator k=m_entity_keys.begin(); k!=m_entity_keys.end(); ++k) {
      const EntityKey id = *k;
      Point center;
      const double *c = coord(id);
      for (unsigned j=0; j<dimension; ++j) {
        center[j] = c[j];
      }
      v.push_back( std::make_pair( Sphere(center,r), EntityProc(id, proc_id)));
    }
  }

  void copy_entities(const EntityProcVec  & keys_to_copy, const std::string  & transfer_name)
  {
    m_bulk_data.modification_begin();
    {
      mesh::EntityProcVec new_entities_to_copy(keys_to_copy.size());
      for (size_t i=0; i<keys_to_copy.size(); ++i) {
        // convert from EntityProc based on EntityKey to EntityProc based on raw Entity.
        const EntityProc key_proc = keys_to_copy[i];
        const EntityKey       key = key_proc.id();
        const unsigned       proc = key_proc.proc();
        const Entity            e = entity(key);
        const mesh::EntityProc ep( e, proc);
        new_entities_to_copy[i] = ep;
      }
      m_entities_currently_ghosted.insert(m_entities_currently_ghosted.end(),
          new_entities_to_copy.begin(),
          new_entities_to_copy.end());

      std::sort(m_entities_currently_ghosted.begin(), m_entities_currently_ghosted.end());
      mesh::EntityProcVec::iterator del = std::unique(m_entities_currently_ghosted.begin(), m_entities_currently_ghosted.end());
      m_entities_currently_ghosted.resize(std::distance(m_entities_currently_ghosted.begin(), del));
    }
    {
      m_bulk_data.change_ghosting(*m_transfer_entity_ghosting,
          m_entities_currently_ghosted);

      std::vector<mesh::EntityKey> receive;
      std::vector<mesh::EntityProc> send;
      m_transfer_entity_ghosting->receive_list( receive );
      m_transfer_entity_ghosting->send_list( send );
    }
    m_mesh_modified = true;
    m_bulk_data.modification_end();
  }

  void update_values()
  {
    std::vector<const mesh::FieldBase *> fields(m_values_field.begin(), m_values_field.end());
    if (m_mesh_modified) {
      // Copy coordinates to the newly ghosted nodes
      m_mesh_modified = false;
      fields.push_back(&m_coordinates_field);
    }
    mesh::communicate_field_data( *m_transfer_entity_ghosting , fields);
    mesh::copy_owned_to_shared  (  m_bulk_data, fields );
  }

  // Needed for LinearInterpoate and FEInterpolation
  const double *coord(const EntityKey k) const
  {
    const mesh::Entity e = entity(k);
    return static_cast<const double*>(stk::mesh::field_data(m_coordinates_field, e));
  }

  const double *value(const EntityKey k, const unsigned i=0) const
  {
    const mesh::Entity  e = entity(k);
    mesh::FieldBase *val=m_values_field[i];
    return static_cast<const double*>(stk::mesh::field_data(*val, e));
  }

  double *value(const EntityKey k, const unsigned i=0)
  {
    const mesh::Entity  e = entity(k);
    mesh::FieldBase *val=m_values_field[i];
    return static_cast<double*>(stk::mesh::field_data(*val, e));
  }

  unsigned  value_size(const EntityKey k, const unsigned i=0) const
  {
    const mesh::Entity         e = entity(k);
    const mesh::FieldBase &field = *m_values_field[i];
    const mesh::Bucket    &bucket= m_bulk_data.bucket(e);

    const unsigned bytes = field_bytes_per_entity(field, bucket);
    const unsigned bytes_per_entry = field.data_traits().size_of;
    const unsigned num_entry = bytes/bytes_per_entry;

    ThrowRequireMsg (bytes == num_entry * bytes_per_entry,
        __FILE__<<":"<<__LINE__<<" Error:" <<"  bytes:" <<bytes<<"  num_entry:" <<num_entry
        <<"  bytes_per_entry:" <<bytes_per_entry);
    return  num_entry;
  }

  unsigned      num_values() const
  { return m_values_field.size(); }

  struct Record { virtual ~Record(){} };

private :
  STKNode ();
  STKNode(const STKNode &M);
  STKNode &operator=(const STKNode&);

  mesh::BulkData                        &m_bulk_data;
  bool                               m_mesh_modified;
  const ParallelMachine                       m_comm;
  const double                          m_sphere_rad;
  const EntityKeySet                   m_entity_keys;
  const mesh::FieldBase         &m_coordinates_field;
  const std::vector<mesh::FieldBase*> m_values_field;

  mesh::Ghosting       *m_transfer_entity_ghosting;
  mesh::EntityProcVec   m_entities_currently_ghosted;

  typedef  boost::shared_ptr<Record>      RecordPtr;
  typedef  std::map<EntityKey,RecordPtr>  RecordMap;
  RecordMap                               m_record_map;

  Entity entity(const EntityKey k) const
  {
    return m_bulk_data.get_entity(k);
  }

  static EntityKeySet entity_keys (const mesh::BulkData &bulk_data, const EntityVec &entities)
  {
    EntityKeySet tmp_keys;
    for (EntityVec::const_iterator e=entities.begin(); e!=entities.end(); ++e) {
      const mesh::EntityKey k = bulk_data.entity_key(*e);
      tmp_keys.insert(k);
    }
    return tmp_keys;
  }
};

} // transfer
} // stk
#endif
