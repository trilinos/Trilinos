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
#ifndef STK_MD_MESH_H
#define STK_MD_MESH_H

#include <memory>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

namespace stk {
namespace transfer {

class MDMesh {
public :
  typedef Intrepid::FieldContainer<double>        MDArray;
  typedef unsigned                                Entity;
  typedef std::vector<Entity>                     EntityVec;
  typedef unsigned                                EntityKey;
  typedef std::set   <EntityKey>                  EntityKeySet;
  typedef search::IdentProc<EntityKey, unsigned>  EntityProc;
  typedef std::vector<EntityProc>                 EntityProcVec;

  typedef search::Point<float>   Point;
  typedef search::Sphere<float>  Sphere;

  typedef std::pair<Sphere,EntityProc> BoundingBox;

  MDMesh(
      MDArray               &val,
      const MDArray               &coord,
      const double                 initial_radius,
      const stk::ParallelMachine   comm) :
    m_num_nodes         (coord.dimension(0)),
    m_spatial_dim       (coord.dimension(1)),
    m_num_values        (val  .dimension(1)),
    m_sphere_rad        (initial_radius    ),
    m_coordinates_field (coord),
    m_values_field      (val)  ,
    m_comm              (comm)
  {}

  // Needed for STK Transfer
  ParallelMachine comm() const {return m_comm;}

  void bounding_boxes (std::vector< std::pair<Sphere,EntityProc> > &v) const
  {
    const float r = m_sphere_rad;

    v.clear();

    const int proc_id = parallel_machine_rank(comm());

    for (unsigned i=0; i!=m_num_nodes; ++i) {

      Point center;
      const double *c = &m_coordinates_field(i,0);
      for (unsigned j=0; j<m_spatial_dim; ++j) {
        center[j] = c[j];
      }
      v.emplace_back( Sphere( center, r), EntityProc(i,proc_id));
    }
  }

  void update_values() {}

  // Needed for LinearInterpoate
  const double *coord(const EntityKey k) const
  {
    return &m_coordinates_field(k,0);
  }

  const double *value(const EntityKey k, const unsigned i=0) const
  {
     return &m_values_field(k, i);
  }

  double *value(const EntityKey k, const unsigned i=0)
  {
     return &m_values_field(k, i);
  }

  unsigned      value_size(const EntityKey e, const unsigned i=0) const
  { return 1; }

  unsigned      num_values() const
  {  return m_num_values;  }

  struct Record { virtual ~Record(){} };
  template <class T> T* database(const EntityKey k) {
    typename RecordMap::const_iterator i = m_record_map.find(k);
    if (i == m_record_map.end()) {
      RecordPtr record(new T());
      typename RecordMap::value_type v(k,record);
      i = m_record_map.insert(v).first;
    }
    T *record = dynamic_cast<T*>(i->second.get());
    ThrowRequireMsg (record,__FILE__<<":"<<__LINE__<<" Dynamic Cast failed in MDMesh::database ");
    return record;
  }


private :
  MDMesh ();
  MDMesh(const MDMesh &M);
  MDMesh &operator=(const MDMesh&);

  const unsigned                         m_num_nodes;
  const unsigned                       m_spatial_dim;
  const unsigned                        m_num_values;
  const double                          m_sphere_rad;
  const MDArray                 &m_coordinates_field;
  MDArray                      &m_values_field;
  const ParallelMachine                       m_comm;

  typedef std::shared_ptr<Record>           RecordPtr;
  typedef std::map<EntityKey,RecordPtr>     RecordMap;
  RecordMap                                 m_record_map;
};

} // transfer
} // stk

#endif

