// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef FromMesh_h
#define FromMesh_h

#include <string>
#include <vector>
#include <utility>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>

namespace percept{

//--------------------------------------------------------------------
/**
 *  @class  FromMesh
 *  @author Stefan P. Domino
 *  @date   July, 2013
 *  @brief  Class that holds source mesh for transfer.
 *
 */
//--------------------------------------------------------------------

class FromMeshAdaptor;

template<class F = FromMeshAdaptor>
class FromMesh {
public :

  typedef F Adaptor;

  typedef stk::mesh:: Entity                                 Entity;
  typedef std::vector<Entity>                                EntityVec;
  typedef stk::mesh:: EntityKey                              EntityKey;
  typedef std::set   <EntityKey>                             EntityKeySet;
  typedef stk::search::IdentProc<EntityKey, unsigned> EntityProc;
  typedef std::vector<EntityProc>                            EntityProcVec;

  typedef std::pair< stk::search::Box<double>,EntityProc> BoundingBox;

  FromMesh(stk::mesh::BulkData    &bulkData,
           stk::mesh::FieldBase * coordinates,
           stk::mesh::FieldBase * field,
           const stk::ParallelMachine   comm) :
    fromMetaData_   (bulkData.mesh_meta_data()),
    fromBulkData_   (bulkData),
    fromcoordinates_(coordinates),
    fromFields_      (1,field),
    comm_           (comm),
    mesh_modified_  (false),
    ghosting_       (0),
    ghosting_map_   ()
  {}

  ~FromMesh(){};

  stk::ParallelMachine comm() const {return comm_;}

  void update_ghosting(const EntityProcVec &entity_keys);

  void update_values();

  void bounding_boxes (std::vector<BoundingBox> &v) const;

  Entity entity(const EntityKey k) const;

  const stk::mesh::MetaData    &fromMetaData_;
  stk::mesh::BulkData    &fromBulkData_;
  stk::mesh::FieldBase *fromcoordinates_;
  std::vector<stk::mesh::FieldBase *> fromFields_;
  const stk::ParallelMachine   comm_;

  bool                          mesh_modified_;
  stk::mesh::Ghosting          *ghosting_;
  stk::mesh::EntityProcVec      ghosting_map_;
};

// default adaptor
class FromMeshAdaptor {
public:

  template<class FM>
  static void
  get_entities(const FM& fromMesh, stk::mesh::EntityRank rank, std::vector<stk::mesh::Entity>& vecEntity)
  {
    stk::mesh::get_selected_entities(stk::mesh::Selector(fromMesh.fromMetaData_.locally_owned_part()), fromMesh.fromBulkData_.buckets(rank), vecEntity);
  }

  template<class FM>
  static void
  modify_bounding_box(const FM& /*fromMesh*/,
                      stk::search::Point<double>& /*min_corner*/,
                      stk::search::Point<double>& /*max_corner*/)
  {
  }

};


} // namespace percept

#include "FromMeshDef.hpp"

#endif
