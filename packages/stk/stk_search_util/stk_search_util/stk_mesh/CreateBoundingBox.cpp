/*------------------------------------------------------------------------*/
/*                 Copyright 2010 - 2011 Sandia Corporation.              */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <assert.h>

#include <stk_search_util/stk_mesh/CreateBoundingBox.hpp>

#include <stk_search/BoundingBox.hpp>

#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/EntityKey.hpp>

static const size_t NODE_RANK = stk::mesh::fem::FEMMetaData::NODE_RANK;

#define ct_assert(e) extern char (*ct_assert(void)) [sizeof(char[1 - 2*!(e)])]

typedef stk::search::ident::IdentProc<stk::mesh::EntityKey, unsigned> IdentProc;
typedef stk::search::box::AxisAlignedBoundingBox<IdentProc, double, 3> AxisAlignedBoundingBox3D;
typedef stk::search::box::PointBoundingBox<IdentProc, double, 3> PointBoundingBox3D;

namespace {
void build_node_axis_bbox(stk::mesh::Part &part,
                          stk::mesh::BulkData &bulk_data,
                          stk::mesh::EntityRank type,
                          CartesianField *coordinates,
                          std::vector<AxisAlignedBoundingBox3D> &box_vector,
                          const stk::search_util::Op &op);

void build_axis_bbox(stk::mesh::Part &part,
                     stk::mesh::BulkData &bulk_data,
                     stk::mesh::EntityRank type,
                     CartesianField *coordinates,
                     std::vector<AxisAlignedBoundingBox3D> &box_vector,
                     const stk::search_util::Op &op);

void build_node_cent_bbox(stk::mesh::Part &part,
                          stk::mesh::BulkData &bulk_data,
                          stk::mesh::EntityRank type,
                          CartesianField *coordinates,
                          std::vector<PointBoundingBox3D> &box_vector);

void build_cent_bbox(stk::mesh::Part &part,
                     stk::mesh::BulkData &bulk_data,
                     stk::mesh::EntityRank type,
                     CartesianField *coordinates,
                     std::vector<PointBoundingBox3D> &box_vector);
}

namespace stk {
namespace search_util {

// The adjust_box should be implemented as an Op rather than a scale and offset parameter.
void build_axis_aligned_bbox(stk::mesh::BulkData &bulk_data, stk::mesh::EntityRank type,
                             CartesianField *coordinates,
                             std::vector<AxisAlignedBoundingBox3D> &box_vector,
                             bool use_universal_part,
                             const stk::search_util::Op &op)
{
  // Build an axis-aligned bounding box of each entity in the mesh...

  // NOTE: Unless 'use_universal_part' is true, then a
  // 'type' of 'Node' will search for a nodeset, not all nodes in
  // the model.

  const stk::mesh::MetaData& meta_data = stk::mesh::MetaData::get(bulk_data);
  const stk::mesh::fem::FEMMetaData & fem = stk::mesh::fem::FEMMetaData::get(meta_data);
  const stk::mesh::EntityRank side_rank = fem.side_rank();

  if (use_universal_part && type == NODE_RANK) {
    stk::mesh::Part &universal = meta_data.universal_part();
    build_node_axis_bbox(universal, bulk_data, type, coordinates, box_vector, op);
  } else if (use_universal_part && type == side_rank) {
    stk::mesh::Part *skin = meta_data.get_part("skin");
    build_axis_bbox(*skin, bulk_data, type, coordinates, box_vector, op);
  } else {
    const stk::mesh::PartVector & all_parts = meta_data.get_parts();
    for ( stk::mesh::PartVector::const_iterator ip = all_parts.begin();
          ip != all_parts.end(); ++ip ) {
      stk::mesh::Part * const part = *ip;
      if ( part->primary_entity_rank() == type ) {
        if (type == NODE_RANK) {
          build_node_axis_bbox(*part, bulk_data, type, coordinates, box_vector, op);
        } else {
          if (type == side_rank && part->name() == "skin")
            continue;
          build_axis_bbox(*part, bulk_data, type, coordinates, box_vector, op);
        }
      }
    }
  }
}

void build_centroid_bbox(stk::mesh::BulkData &bulk_data,  stk::mesh::EntityRank type,
                         CartesianField *coordinates,
                         std::vector<PointBoundingBox3D> &box_vector,
                         bool use_universal_part)
{
  // NOTE: Unless 'nodes_use_universal_part' is true, then a
  // 'type' of 'Node' will search for a nodeset, not all nodes in
  // the model.

  // The box for this case is the centroid of each entity in the mesh...
  const stk::mesh::MetaData& meta_data = stk::mesh::MetaData::get(bulk_data);

  const stk::mesh::fem::FEMMetaData & fem = stk::mesh::fem::FEMMetaData::get(meta_data);
  const stk::mesh::EntityRank side_rank = fem.side_rank();

  if (use_universal_part && type == NODE_RANK) {
    stk::mesh::Part &universal = meta_data.universal_part();
    build_node_cent_bbox(universal, bulk_data, type, coordinates, box_vector);
  } else if (use_universal_part && type == side_rank) {
    stk::mesh::Part *skin = meta_data.get_part("skin");
    build_cent_bbox(*skin, bulk_data, type, coordinates, box_vector);
  } else {
    const stk::mesh::PartVector & all_parts = meta_data.get_parts();
    for ( stk::mesh::PartVector::const_iterator ip = all_parts.begin();
          ip != all_parts.end(); ++ip ) {
      stk::mesh::Part * const part = *ip;
      if ( part->primary_entity_rank() == type ) {
        if (type == NODE_RANK) {
          build_node_cent_bbox(*part, bulk_data, type, coordinates, box_vector);
        } else {
          if (type == side_rank && part->name() == "skin")
            continue;
          build_cent_bbox(*part, bulk_data, type, coordinates, box_vector);
        }
      }
    }
  }
}

void OffsetScaleOp::operator()(Data &xmin, Data &ymin, Data &zmin,
                               Data &xmax, Data &ymax, Data &zmax) const {
  if (m_offset != 0.0 || m_scale != 0.0) {
    const double dx = (xmax - xmin);
    const double dy = (ymax - ymin);
    const double dz = (zmax - zmin);

    double delta = dx;
    if (delta < dy)
      delta = dy;
    if (delta < dz)
      delta = dz;

    delta = delta * m_scale + m_offset;

    xmin -= delta;
    ymin -= delta;
    zmin -= delta;
    xmax += delta;
    ymax += delta;
    zmax += delta;
  }
}

void OffsetScaleOp::operator()(AxisAlignedBoundingBox3D &box) const {
  if (m_offset != 0.0 || m_scale != 0.0) {
    double dx = box.upper(0) - box.lower(0);
    double dy = box.upper(1) - box.lower(1);
    double dz = box.upper(2) - box.lower(2);

    double delta = dx;
    if (delta < dy)
      delta = dy;
    if (delta < dz)
      delta = dz;

    delta = delta * m_scale + m_offset;

    box.expand(delta);
  }
}

} // namespace search_util
} // namespace stk

namespace {

void build_node_axis_bbox(stk::mesh::Part &part,
                          stk::mesh::BulkData &bulk_data,
                          stk::mesh::EntityRank type,
                          CartesianField *coordinates,
                          std::vector<AxisAlignedBoundingBox3D> &box_vector,
                          const stk::search_util::Op &op)
{
  const stk::mesh::MetaData& meta_data = stk::mesh::MetaData::get(bulk_data);

  std::vector<stk::mesh::Entity *> entities;
  stk::mesh::Selector selector = part & ( meta_data.locally_owned_part() | meta_data.globally_shared_part() );
  get_selected_entities(selector, bulk_data.buckets(type), entities);
  size_t num_entities = entities.size();

  for (size_t i = 0; i < num_entities; ++i) {
    AxisAlignedBoundingBox3D   domain;
    ct_assert(sizeof(domain.key.ident) >= sizeof(stk::mesh::EntityKey));
    domain.key.ident = entities[i]->key();

    double *fld_data = (double*)stk::mesh::field_data(*coordinates, *entities[i]);
    assert(fld_data != NULL);

    domain.set_box(fld_data);
    op(domain);
    box_vector.push_back(domain);
  }
}

void build_axis_bbox(stk::mesh::Part &part,
                     stk::mesh::BulkData &bulk_data,
                     stk::mesh::EntityRank type,
                     CartesianField *coordinates,
                     std::vector<AxisAlignedBoundingBox3D> &box_vector,
                     const stk::search_util::Op &op)
{
  const stk::mesh::MetaData& meta_data = stk::mesh::MetaData::get(bulk_data);
  stk::mesh::fem::FEMMetaData * fem_meta = const_cast<stk::mesh::fem::FEMMetaData *>(meta_data.get_attribute<stk::mesh::fem::FEMMetaData>());

  const CellTopologyData * cell_topo = NULL;
  cell_topo = fem_meta->get_cell_topology(part).getCellTopologyData();
  if (fem_meta && !cell_topo) cell_topo = fem_meta->get_cell_topology(part).getCellTopologyData();
  if (cell_topo == NULL)
    return;

  const int nodes_per_entity = cell_topo->node_count;

  std::vector<stk::mesh::Entity *> entities;
  stk::mesh::Selector selector = part & ( meta_data.locally_owned_part() | meta_data.globally_shared_part() );
  get_selected_entities(selector, bulk_data.buckets(type), entities);
  size_t num_entities = entities.size();

  for (size_t i = 0; i < num_entities; ++i) {
    AxisAlignedBoundingBox3D   domain;
    ct_assert(sizeof(domain.key.ident) >= sizeof(stk::mesh::EntityKey));
    domain.key.ident = entities[i]->key();

    const stk::mesh::PairIterRelation entity_nodes = entities[i]->relations(NODE_RANK);
    assert(static_cast<int>(entity_nodes.size()) == nodes_per_entity);

    double *fld_data = (double*)stk::mesh::field_data(*coordinates, *entity_nodes[0].entity());
    assert(fld_data != NULL);

    AxisAlignedBoundingBox3D::Data bbox[6];
    bbox[0]= fld_data[0];
    bbox[1] = fld_data[1];
    bbox[2] = fld_data[2];

    bbox[3] = fld_data[0];
    bbox[4] = fld_data[1];
    bbox[5] = fld_data[2];

    for (int j = 1; j < nodes_per_entity; ++j) {
      fld_data = (double*)stk::mesh::field_data(*coordinates, *entity_nodes[j].entity());
      assert(fld_data != NULL);
      bbox[0] = fld_data[0] < bbox[0] ? fld_data[0] : bbox[0];
      bbox[1] = fld_data[1] < bbox[1] ? fld_data[1] : bbox[1];
      bbox[2] = fld_data[2] < bbox[2] ? fld_data[2] : bbox[2];

      bbox[3] = fld_data[0] > bbox[3] ? fld_data[0] : bbox[3];
      bbox[4] = fld_data[1] > bbox[4] ? fld_data[1] : bbox[4];
      bbox[5] = fld_data[2] > bbox[5] ? fld_data[2] : bbox[5];
    }

    domain.set_box(bbox);
    op(domain);
    box_vector.push_back(domain);
  }
}

void build_node_cent_bbox(stk::mesh::Part &part,
                          stk::mesh::BulkData &bulk_data,
                          stk::mesh::EntityRank type,
                          CartesianField *coordinates,
                          std::vector<PointBoundingBox3D> &box_vector)
{
  const stk::mesh::MetaData& meta_data = stk::mesh::MetaData::get(bulk_data);

  std::vector<stk::mesh::Entity *> entities;
  stk::mesh::Selector selector = part & ( meta_data.locally_owned_part() | meta_data.globally_shared_part() );
  get_selected_entities(selector, bulk_data.buckets(type), entities);
  size_t num_entities = entities.size();

  for (size_t i = 0; i < num_entities; ++i) {
    PointBoundingBox3D   p;
    ct_assert(sizeof(p.key.ident) >= sizeof(stk::mesh::EntityKey));
    p.key.ident = entities[i]->key();

    double *fld_data = (double*)stk::mesh::field_data(*coordinates, *entities[i]);
    assert(fld_data != NULL);

    p.set_center(fld_data);
    box_vector.push_back(p);
  }
}

void build_cent_bbox(stk::mesh::Part &part,
                     stk::mesh::BulkData &bulk_data,
                     stk::mesh::EntityRank type,
                     CartesianField *coordinates,
                     std::vector<PointBoundingBox3D> &box_vector)
{
  const stk::mesh::MetaData& meta_data = stk::mesh::MetaData::get(bulk_data);

  std::vector<stk::mesh::Entity *> entities;
  stk::mesh::Selector selector = part & ( meta_data.locally_owned_part() | meta_data.globally_shared_part() );
  get_selected_entities(selector, bulk_data.buckets(type), entities);
  size_t num_entities = entities.size();

  for (size_t i = 0; i < num_entities; ++i) {
    PointBoundingBox3D   p;
    ct_assert(sizeof(p.key.ident) >= sizeof(stk::mesh::EntityKey));
    p.key.ident = entities[i]->key();

    p.center[0] = 0;
    p.center[1] = 0;
    p.center[2] = 0;

    const stk::mesh::PairIterRelation entity_nodes = entities[i]->relations(NODE_RANK);
    const size_t nodes_per_entity = entity_nodes.size();
    for (size_t j = 0; j < nodes_per_entity; ++j) {
      double *fld_data = (double*)stk::mesh::field_data(*coordinates, *entity_nodes[j].entity());
      assert(fld_data != NULL);
      p.center[0] += fld_data[0];
      p.center[1] += fld_data[1];
      p.center[2] += fld_data[2];
    }
    p.center[0] /= nodes_per_entity;
    p.center[1] /= nodes_per_entity;
    p.center[2] /= nodes_per_entity;
    box_vector.push_back(p);
  }
}

}
