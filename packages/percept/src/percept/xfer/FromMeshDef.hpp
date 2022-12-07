// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef FromMeshDef_hpp
#define FromMeshDef_hpp

#include <percept/PerceptMesh.hpp>

namespace percept
{

template<class F>
stk::mesh::Entity
FromMesh<F>::entity(const stk::mesh::EntityKey k) const {
  const stk::mesh::Entity  e = fromBulkData_.get_entity(k);
  return e;
}


template<class F>
void
FromMesh<F>::update_ghosting(const EntityProcVec &entity_keys) {
  ghosting_map_.resize(0);
  for (size_t i=0; i<entity_keys.size(); ++i) {
    //convert from EntityProc based on EntityKey to EntityProc based on raw Entity.
    const EntityProc&              key_proc = entity_keys[i];
    const stk::mesh::EntityKey          key = key_proc.id();
    const unsigned                     proc = key_proc.proc();
    const stk::mesh::Entity               e = entity(key);
    const stk::mesh::EntityProc ep( e, proc);
    //VERIFY_OP_ON(fromBulkData_.bucket(e).owned(), ==, true, "not owned");
    if (fromBulkData_.bucket(e).owned())
      ghosting_map_.push_back( ep);
  }

  unsigned s = !ghosting_map_.empty();
  stk::all_reduce( comm_, stk::ReduceSum<1>(&s));

  if (s) {
    std::sort(ghosting_map_.begin(), ghosting_map_.end());
    stk::mesh::EntityProcVec::iterator del = std::unique(ghosting_map_.begin(), ghosting_map_.end());
    ghosting_map_.resize(std::distance(ghosting_map_.begin(), del));

    std::string theGhostName = "percept_transfer_ghosting";
    ghosting_ = &fromBulkData_.create_ghosting( theGhostName );

    fromBulkData_.change_ghosting( *ghosting_, ghosting_map_);
    mesh_modified_ = true;
  }
}

template<class F>
void
FromMesh<F>::update_values () {
  if (ghosting_) {
    std::vector<const stk::mesh::FieldBase *> fields;
    for (unsigned ii=0; ii < fromFields_.size(); ++ii)
      fields.push_back(fromFields_[ii]);
    if (mesh_modified_) {
      // Copy coordinates to the newly ghosted nodes
      mesh_modified_ = false;
      fields.push_back(fromcoordinates_);
    }
    stk::mesh::communicate_field_data( *ghosting_ ,    fields);
    stk::mesh::copy_owned_to_shared  (  fromBulkData_, fields);
  }
}

template<class F>
void
FromMesh<F>::bounding_boxes(std::vector<FromMesh<F>::BoundingBox> &v_box) const
{
  const unsigned nDim = fromMetaData_.spatial_dimension();

  stk::search::Point<double> min_corner, max_corner;

  stk::mesh::EntityRank rank = stk::topology::ELEMENT_RANK;

  std::vector<stk::mesh::Entity> vecEntity;
  Adaptor adaptor;
  adaptor.get_entities(*this, rank, vecEntity);

  for (size_t ii = 0; ii < vecEntity.size(); ++ii)
    {
      stk::mesh::Entity entity = vecEntity[ii];

      // initialize max and min
      for (unsigned j = 0; j < nDim; ++j ) {
        min_corner[j] = +std::numeric_limits<double>::max();
        max_corner[j] = -std::numeric_limits<double>::max();
      }

      stk::mesh::Entity const * entity_node_rels = fromBulkData_.begin_nodes(entity);
      int num_entity_nodes = fromBulkData_.num_nodes(entity);
      for ( int ni = 0; ni < num_entity_nodes; ++ni ) {
        stk::mesh::Entity node = entity_node_rels[ni];

        double * coords = static_cast<double*>(stk::mesh::field_data(*fromcoordinates_, node));
        for ( unsigned j=0; j < nDim; ++j ) {
          min_corner[j] = std::min(min_corner[j], coords[j]);
          max_corner[j] = std::max(max_corner[j], coords[j]);
        }
      }

      adaptor.modify_bounding_box(*this, min_corner, max_corner);

      // setup ident
      FromMesh<F>::EntityProc theIdent( fromBulkData_.entity_key(entity), fromBulkData_.parallel_rank());

      FromMesh<F>::BoundingBox theBox( stk::search::Box<double>(min_corner,max_corner), theIdent );
      v_box.push_back(theBox);
    }
  }

} // namespace percept

#endif
