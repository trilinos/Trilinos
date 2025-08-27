// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef ToMesh_h
#define ToMesh_h

#include <string>
#include <vector>
#include <utility>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>

#include <percept/xfer/TransferHelper.hpp>
#include <percept/xfer/FromMesh.hpp>

#include <stk_mesh/base/FieldParallel.hpp>

namespace percept{

class ToMeshAdaptor;

template<class T = ToMeshAdaptor>
class ToMesh {
public :
  typedef T Adaptor;

  typedef stk::mesh:: Entity                                      Entity;
  typedef std::vector<Entity>                                EntityVec;
  typedef stk::mesh:: EntityKey                                   EntityKey;
  typedef std::set   <EntityKey>                             EntityKeySet;
  typedef stk::search::IdentProc<EntityKey, unsigned>      EntityProc;
  typedef std::vector<EntityProc>                            EntityProcVec;
  typedef std::pair< stk::search::Sphere<double>, EntityProc> BoundingBox;

  ToMesh(stk::mesh::BulkData    &bulkData,
   stk::mesh::FieldBase * coordinates,
	 stk::mesh::FieldBase * field,
	 const stk::ParallelMachine   comm,
	 TransferType transferType,
         SrcFieldType srcFieldType=SRC_FIELD,
         const double radius=.0001) :
    toMetaData_(bulkData.mesh_meta_data()),
    toBulkData_(bulkData),
    tocoordinates_(coordinates),
    toFields_   (1,field),
    comm_(comm),
    transferType_(transferType),
    srcFieldType_(srcFieldType),
    radius_(radius)
  {}

  ~ToMesh(){};

  stk::ParallelMachine comm() const {return comm_;}

  void update_values();

  void bounding_boxes (std::vector<BoundingBox> &v) const;

  void set_entity_parametric_coords(const EntityKey& key, const std::vector<double>& coords) {
    par_coords_[key] = coords;
  }

  [[nodiscard]] const std::vector<double>& get_entity_parametric_coords(const EntityKey& key) const {
    if( par_coords_.count(key) != 1 ) {
      throw std::runtime_error("Key not found in database");
    }
    return par_coords_.at(key);
  }

  const stk::mesh::MetaData    &toMetaData_;
  stk::mesh::BulkData    &toBulkData_;
  stk::mesh::FieldBase         *tocoordinates_;
  std::vector<stk::mesh::FieldBase *> toFields_;
  const stk::ParallelMachine            comm_;
  TransferType transferType_;
  SrcFieldType srcFieldType_;
  const double                 radius_;

  typedef std::map<stk::mesh::EntityKey, std::vector<double> > TransferInfo;
  TransferInfo par_coords_;
};

class ToMeshAdaptor {
public:

  typedef std::multimap<stk::mesh::EntityKey, stk::mesh::EntityKey> EntityKeyMap;

  // return true to skip the standard processing
  template<class FromMesh, class ToMesh>
  bool
  apply(ToMesh &/*ToPoints*/,
        const FromMesh  &/*FromElem*/,
        const EntityKeyMap &/*RangeToDomain*/)
  {
    return false;
  }

  // return true to skip the standard processing
  template<class FromMesh, class ToMesh>
  bool
  filter_to_nearest(const EntityKeyMap &/*RangeToDomain*/,
                    const FromMesh  &/*FromElem*/,
                    ToMesh &/*ToPoints*/)
  {
    return false;
  }

  void
  modify_bounding_box(stk::search::Point<double>& /*min_corner*/,
                      stk::search::Point<double>& /*max_corner*/) {}

  template<class ToMesh>
  void
  modify_selector(const ToMesh& /*ToPoints*/, stk::mesh::Selector &/*sel*/) {}

};

} // namespace percept

#include "ToMeshDef.hpp"

#endif
