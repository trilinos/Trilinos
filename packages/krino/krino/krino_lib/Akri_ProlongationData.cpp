// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CDFEM_Support.hpp>
#include <Akri_ProlongationData.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_FieldRef.hpp>

#include <stk_mesh/base/CommunicateMeshTypes.hpp>
#include <stk_util/parallel/CommSparse.hpp>

#include <cmath>

namespace krino{

ProlongationElementData::ProlongationElementData(const stk::mesh::BulkData& stk_mesh, const std::vector<const ProlongationElementData *> & children_data, const std::vector< std::vector<double> > & children_intg_wts)
{
  ThrowAssert(!children_data.empty() && children_data.size() == children_intg_wts.size());

  const unsigned num_intg_pts = children_intg_wts[0].size();

  // homogenize fields
  const stk::mesh::FieldVector & all_fields = stk_mesh.mesh_meta_data().get_fields();
  my_field_indices.resize(all_fields.size(), -1);
  for ( stk::mesh::FieldVector::const_iterator it = all_fields.begin(); it != all_fields.end() ; ++it )
  {
    const FieldRef field = **it;

    if( field.entity_rank()!=stk::topology::ELEMENT_RANK || !field.type_is<double>() ) continue;

    bool any_child_has_field = false;
    for (auto && child_data : children_data)
    {
      if (NULL != child_data->get_field_data(field))
      {
        any_child_has_field = true;
        break;
      }
    }

    if (any_child_has_field)
    {
      const unsigned field_length = field.length();
      const unsigned field_data_index = my_field_data.size();
      my_field_indices[field.field().mesh_meta_data_ordinal()] = field_data_index;
      my_field_data.resize(field_data_index+field_length, 0.0);

      // TODO: Add a method to distinguish between vector fields and gauss point fields.
      const bool data_is_gauss_pt_field = (num_intg_pts == field_length);

      if (data_is_gauss_pt_field)
      {
        double tot_child_sum = 0.;
        double tot_child_vol = 0.;

        for (unsigned n = 0; n < children_data.size(); ++n)
        {
          const double * subelem_field_data = children_data[n]->get_field_data(field);
          if(NULL == subelem_field_data)
          {
            continue;
          }

          const std::vector<double> & child_intg_wts = children_intg_wts[n];
          ThrowAssertMsg(child_intg_wts.size() == num_intg_pts, "Children have different integration rules.");

          for (unsigned j=0; j<num_intg_pts; ++j)
          {
            tot_child_sum += subelem_field_data[j] * child_intg_wts[j];
            tot_child_vol += child_intg_wts[j];
          }
        }

        const double tot_child_avg = tot_child_sum / tot_child_vol;
        for (unsigned i=0; i<field_length; ++i)
        {
          my_field_data[field_data_index+i] = tot_child_avg;
        }
      }
      else // vector field (includes scalar case)
      {
        double tot_child_vol = 0.;

        for (unsigned n = 0; n < children_data.size(); ++n )
        {
          const double * subelem_field_data = children_data[n]->get_field_data(field);
          if(NULL == subelem_field_data)
          {
            continue;
          }

          const std::vector<double> & child_intg_wts = children_intg_wts[n];
          // We could relax this assertion if we had another way to distinguish gauss point fields from vector fields
          ThrowAssertMsg(child_intg_wts.size() == num_intg_pts, "Children have different integration rules.");

          double child_vol = 0.;
          for (unsigned j=0; j<num_intg_pts; ++j)
          {
            child_vol += child_intg_wts[j];
          }
          tot_child_vol += child_vol;

          for (unsigned i=0; i<field_length; ++i)
          {
            my_field_data[field_data_index+i] += child_vol * subelem_field_data[i];
          }
        }

        for (unsigned i=0; i<field_length; ++i)
        {
          my_field_data[field_data_index+i] /= tot_child_vol;
        }
      }
    }
  }
}

void
ProlongationData::save_fields(const stk::mesh::BulkData& stk_mesh, stk::mesh::Entity entity)
{
  const stk::mesh::EntityRank entity_rank = stk_mesh.entity_rank(entity);
  const stk::mesh::FieldVector & all_fields = stk_mesh.mesh_meta_data().get_fields();
  my_field_indices.resize(all_fields.size(), -1);
  for ( auto&& field_ptr : all_fields )
  {
    const FieldRef field = *field_ptr;

    if( field.entity_rank()!=entity_rank || !field.type_is<double>() ) continue;

    double * val = field_data<double>(field, entity);
    const bool has_field = (NULL != val);

    if (has_field)
    {
      const unsigned field_length = field.length();
      my_field_indices[field.field().mesh_meta_data_ordinal()] = my_field_data.size();
      for (unsigned i=0; i<field_length; ++i)
      {
        my_field_data.push_back(val[i]);
      }
    }
  }
}

void
ProlongationData::restore_fields(const stk::mesh::BulkData& stk_mesh, stk::mesh::Entity entity) const
{
  stk::mesh::EntityRank entity_rank = stk_mesh.entity_rank(entity);
  const stk::mesh::FieldVector & all_fields = stk_mesh.mesh_meta_data().get_fields();
  for ( stk::mesh::FieldVector::const_iterator it = all_fields.begin(); it != all_fields.end() ; ++it )
  {
    const FieldRef field = **it;

    if( field.entity_rank()!=entity_rank || !field.type_is<double>() ) continue;
    const unsigned field_length = field.length();

    double * val = field_data<double>(field, entity);
    const double * prolong_field = get_field_data(field);

    if (nullptr == val) continue;
    if(nullptr == prolong_field)
    {
      std::stringstream err_msg;
      err_msg << "Missing prolongation field data when restoring fields on entity:\n";
      err_msg << stk_mesh.identifier(entity) << " of rank " << stk_mesh.entity_rank(entity);
      err_msg << " with parts:\n  ";
      for(auto && part : stk_mesh.bucket(entity).supersets())
      {
        err_msg << part->name() << ", ";
      }
      err_msg << "\n";
      err_msg << "Missing field data for field " << field.name() << "\n";
      throw std::runtime_error(err_msg.str());
    }

    for (unsigned i=0; i<field_length; ++i) val[i] = prolong_field[i];
  }
}

std::vector<unsigned>
ProlongationNodeData::get_fields_on_node(const stk::mesh::BulkData& mesh, stk::mesh::Entity entity)
{
  const stk::mesh::FieldVector & all_fields = mesh.mesh_meta_data().get_fields();
  std::vector<unsigned> entity_fields;
  entity_fields.reserve(all_fields.size());
  for ( auto && fieldPtr : all_fields )
  {
    const FieldRef field = *fieldPtr;
    if (field.field().entity_rank() == stk::topology::NODE_RANK && nullptr != field_data<double>(field, entity))
      entity_fields.push_back(field.field().mesh_meta_data_ordinal());
  }

  std::sort(entity_fields.begin(), entity_fields.end());
  return entity_fields;
}

Vector3d
ProlongationNodeData::get_node_coordinates(const CDMesh & mesh, stk::mesh::Entity node)
{
  const double * coordsPtr = field_data<double>(mesh.get_coords_field(), node);
  ThrowAssert(coordsPtr);
  Vector3d coords(coordsPtr, mesh.spatial_dim());
  FieldRef cdfemSnapDispField = mesh.get_cdfem_support().get_cdfem_snap_displacements_field();
  if (cdfemSnapDispField.valid())
  {
    FieldRef oldCdfemSnapDispField = cdfemSnapDispField.field_state(stk::mesh::StateOld);
    double * cdfemSnapDispPtr = field_data<double>(cdfemSnapDispField, node);
    double * oldCdfemSnapDispPtr = field_data<double>(oldCdfemSnapDispField, node);
    if (nullptr != cdfemSnapDispPtr)
      coords += Vector3d(oldCdfemSnapDispPtr, mesh.spatial_dim()) - Vector3d(cdfemSnapDispPtr, mesh.spatial_dim());
  }
  return coords;
}

std::vector<unsigned>
ProlongationNodeData::get_node_io_parts(const stk::mesh::BulkData& stk_mesh, stk::mesh::Entity entity)
{
  // This list of parts is used to determine if a node needs to be ALE prolonged
  stk::mesh::PartVector node_parts = filter_non_io_parts(stk_mesh.bucket(entity).supersets());
  std::vector<unsigned> node_part_ids;
  node_part_ids.reserve(node_parts.size());
  for (auto * node_part : node_parts)
    node_part_ids.push_back(node_part->mesh_meta_data_ordinal());
  std::sort(node_part_ids.begin(), node_part_ids.end());
  return node_part_ids;
}

ProlongationNodeData::ProlongationNodeData(const CDMesh & mesh, stk::mesh::Entity node, bool communicate_me_to_all_sharers)
  : ProlongationPointData(get_node_coordinates(mesh, node)),
    my_entityId(mesh.stk_bulk().identifier(node)),
    myCommunicateMeToAllSharersFlag(communicate_me_to_all_sharers)
{
  const stk::mesh::BulkData& stk_mesh = mesh.stk_bulk();

  my_fields = get_fields_on_node(stk_mesh, node);
  my_ioparts = get_node_io_parts(stk_mesh, node);

  save_fields(stk_mesh, node);
}

ProlongationPointData::ProlongationPointData(const CDMesh & mesh, const FacetDistanceQuery & facet_dist_query,
    const std::vector<const ProlongationNodeData *> & facet_nodes)
  : my_coordinates(facet_dist_query.closest_point())
{
  const Vector3d node_wts = facet_dist_query.closest_point_weights();

  // interpolate fields
  const stk::mesh::FieldVector & all_fields = mesh.stk_meta().get_fields();
  my_field_indices.resize(all_fields.size(), -1);
  for ( stk::mesh::FieldVector::const_iterator it = all_fields.begin(); it != all_fields.end() ; ++it )
  {
    const FieldRef field = **it;

    if( field.entity_rank()!=stk::topology::NODE_RANK || !field.type_is<double>() ) continue;

    const unsigned field_length = field.length();

    bool any_node_has_field = false;
    for (auto && facet_node : facet_nodes)
    {
      if (NULL != facet_node->get_field_data(field))
      {
        any_node_has_field = true;
        break;
      }
    }

    if (any_node_has_field)
    {
      const unsigned field_data_index = my_field_data.size();
      my_field_indices[field.field().mesh_meta_data_ordinal()] = field_data_index;
      my_field_data.resize(field_data_index+field_length, 0.0);

      double node_wt_sum = 0.0;
      for (unsigned n = 0; n < facet_nodes.size(); ++n)
      {
        const double * node_field_data = facet_nodes[n]->get_field_data(field);

        if(node_field_data)
        {
          node_wt_sum += node_wts[n];
          for (unsigned i=0; i<field_length; ++i)
          {
            my_field_data[field_data_index+i] += node_wts[n] * node_field_data[i];
          }
        }
      }

      for(unsigned i=0; i<field_length; ++i)
      {
        my_field_data[field_data_index+i] /= node_wt_sum;
      }
    }
  }
}

static
void pack_nodes_that_need_to_be_communicated_to_sharers(const stk::mesh::BulkData & mesh, const EntityProlongationNodeMap & proc_prolong_nodes, stk::CommSparse &commSparse)
{
  std::vector<int> nodeSharedProcs;
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (auto entry : proc_prolong_nodes)
    {
      const ProlongationNodeData & prolongNode = *entry.second;
      if (prolongNode.communicate_me_to_all_sharers())
      {
        stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, prolongNode.entityId());
        mesh.comm_shared_procs(node, nodeSharedProcs);
        for (int procId : nodeSharedProcs)
          prolongNode.pack_into_buffer(commSparse.send_buffer(procId));
      }
    }
  });
}

static
void pack_facet_prolong_nodes(const stk::mesh::BulkData & mesh, const ProlongFacetVec & proc_prolong_facets, const std::vector<BoundingBox> & proc_target_bboxes, stk::CommSparse &commSparse)
{
  stk::pack_and_communicate(commSparse,[&]()
  {
    for ( int procId=0; procId<commSparse.parallel_size(); ++procId )
    {
      if ( commSparse.parallel_rank() == procId ) continue;  // Don't talk to yourself, it's embarrassing
      stk::CommBuffer & buffer = commSparse.send_buffer(procId);

      for ( auto && node : ProlongationFacet::get_facet_nodes_to_communicate(proc_prolong_facets, proc_target_bboxes[procId]) )
        node->pack_into_buffer(buffer);
    }
  });
}

static
void receive_and_build_prolong_nodes(const CDMesh & mesh, EntityProlongationNodeMap & proc_prolong_nodes, stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      ProlongationNodeData * node = ProlongationNodeData::unpack_from_buffer( buffer, mesh.stk_meta() ); // This calls new to create a new ProlongationNodeData
      EntityProlongationNodeMap::iterator it = proc_prolong_nodes.find(node->entityId());
      if( it == proc_prolong_nodes.end() || it->second == nullptr )
      {
        proc_prolong_nodes[node->entityId()] = node;
      }
      else
      {
        delete node;
      }
    }
  });
}

void ProlongationFacet::communicate_shared_nodes( const CDMesh & mesh, EntityProlongationNodeMap & proc_prolong_nodes )
{
  stk::CommSparse commSparse(mesh.stk_bulk().parallel());
  pack_nodes_that_need_to_be_communicated_to_sharers(mesh.stk_bulk(), proc_prolong_nodes, commSparse);
  receive_and_build_prolong_nodes(mesh, proc_prolong_nodes, commSparse);
}


void
ProlongationFacet::communicate_facet_nodes( const CDMesh & mesh, const ProlongFacetVec & proc_prolong_facets, EntityProlongationNodeMap & proc_prolong_nodes, const std::vector<BoundingBox> & proc_target_bboxes )
{ /* %TRACE[ON]% */ Trace trace__("krino:ProlongationFacet::communicate_facet_nodes()"); /* %TRACE% */
  const int num_procs = mesh.stk_bulk().parallel_size();
  if ( num_procs == 1 ) return;  // Don't talk to yourself, it's embarrassing

  stk::CommSparse commSparse(mesh.stk_bulk().parallel());
  pack_facet_prolong_nodes(mesh.stk_bulk(), proc_prolong_facets, proc_target_bboxes, commSparse);
  receive_and_build_prolong_nodes(mesh, proc_prolong_nodes, commSparse);
}

void
ProlongationNodeData::pack_into_buffer(stk::CommBuffer & b) const
{
  b.pack(my_entityId);

  b.pack(my_coordinates.data(),3);

  const size_t num_fields = my_fields.size();
  b.pack(num_fields);
  b.pack(my_fields.data(), my_fields.size());

  const size_t num_parts = my_ioparts.size();
  b.pack(num_parts);
  b.pack(my_ioparts.data(), my_ioparts.size());

  b.pack(my_field_indices.data(),my_field_indices.size());

  const size_t num_field_data = my_field_data.size();
  b.pack(num_field_data);
  b.pack(my_field_data.data(),my_field_data.size());
}

ProlongationNodeData *
ProlongationNodeData::unpack_from_buffer( stk::CommBuffer & b, const stk::mesh::MetaData & stk_meta )
{
  stk::mesh::EntityId global_id;
  b.unpack(global_id);

  Vector3d coords;
  b.unpack(coords.data(),3);

  size_t num_fields = 0;
  b.unpack(num_fields);
  std::vector<unsigned> node_fields(num_fields);
  b.unpack(node_fields.data(), num_fields);

  size_t num_parts = 0;
  b.unpack(num_parts);
  std::vector<unsigned> node_ioparts(num_parts);
  b.unpack(node_ioparts.data(), num_parts);

  ProlongationNodeData * node = new ProlongationNodeData(global_id, coords, node_fields, node_ioparts);

  const size_t len_field_indices = stk_meta.get_fields().size();
  std::vector<int> & field_indices = node->get_field_indices();
  field_indices.resize(len_field_indices);
  b.unpack(field_indices.data(), len_field_indices);

  size_t num_field_data = 0;
  b.unpack(num_field_data);
  std::vector<double> & field_data = node->get_field_data();
  field_data.resize(num_field_data);
  b.unpack(field_data.data(), num_field_data);

  return node;
}

std::string
ProlongationData::missing_prolongation_fields_for_entity( const CDMesh & mesh, const stk::mesh::Entity dst ) const
{
  std::string missing_fields;
  const FieldSet & ale_prolongation_fields = mesh.get_ale_prolongation_fields();
  for ( auto&& field : ale_prolongation_fields )
  {
    if( !field.type_is<double>() || field.entity_rank() != mesh.stk_bulk().entity_rank(dst) ) continue;

    double * val = field_data<double>(field, dst);
    if (NULL != val && NULL == get_field_data(field))
    {
      missing_fields = missing_fields + " " + field.name();
    }
  }

  return missing_fields;
}

void ProlongationFacet::compute_common_fields()
{
  my_common_fields.clear();
  for (unsigned side_node_index=0; side_node_index<my_prolong_nodes.size(); ++side_node_index)
  {
    const std::vector<unsigned> & node_fields = my_prolong_nodes[side_node_index]->get_fields();
    if (0 == side_node_index)
    {
      my_common_fields = node_fields;
    }
    else
    {
      std::vector<unsigned> working_set;
      working_set.swap(my_common_fields);
      std::set_intersection(working_set.begin(),working_set.end(),node_fields.begin(),node_fields.end(),std::back_inserter(my_common_fields));
    }
  }
}

ProlongationFacet::ProlongationFacet(const CDMesh & mesh, stk::mesh::Entity side)
: my_mesh(mesh)
{
  const stk::mesh::BulkData & stk_mesh = my_mesh.stk_bulk();

  ThrowAssert(stk_mesh.num_elements(side) > 0);
  stk::mesh::Entity elem0 = stk_mesh.begin_elements(side)[0];
  const PhaseTag elem0_phase = mesh.determine_entity_phase(elem0);

  const unsigned num_side_nodes = stk_mesh.bucket(side).topology().base().num_nodes();
  const stk::mesh::Entity* side_nodes = stk_mesh.begin_nodes(side);
  my_prolong_nodes.resize(num_side_nodes);

  for (unsigned side_node_index=0; side_node_index<num_side_nodes; ++side_node_index)
  {
    stk::mesh::Entity side_node = side_nodes[side_node_index];
    const ProlongationNodeData * prolong_node = mesh.fetch_prolong_node(stk_mesh.identifier(side_node));

    if (NULL == prolong_node)
    {
      krinolog << stk::diag::dendl;
      krinolog << "Failed to find node " << stk_mesh.identifier(side_node) << " while stashing facet " << debug_entity(mesh.stk_bulk(), side);
    }
    ThrowRequire(NULL != prolong_node);
    my_prolong_nodes[side_node_index] = prolong_node;
  }

  compute_common_fields();

  ThrowAssert((int)my_prolong_nodes.size() == my_mesh.spatial_dim());
  if (2 == my_prolong_nodes.size())
  {
    my_facet = std::make_unique<Facet2d>( my_prolong_nodes[0]->get_coordinates(), my_prolong_nodes[1]->get_coordinates());
  }
  else
  {
    ThrowAssert(3 == my_prolong_nodes.size());
    my_facet = std::make_unique<Facet3d>( my_prolong_nodes[0]->get_coordinates(), my_prolong_nodes[1]->get_coordinates(), my_prolong_nodes[2]->get_coordinates());
  }
}

ProlongationFacet::ProlongationFacet(const CDMesh & mesh, const std::vector<const ProlongationNodeData *> & prolong_nodes, const std::vector<unsigned> & common_fields)
: my_mesh (mesh),
  my_prolong_nodes(prolong_nodes),
  my_common_fields(common_fields)
{
  ThrowAssert((int)my_prolong_nodes.size() == my_mesh.spatial_dim());
  if (2 == my_prolong_nodes.size())
  {
    my_facet = std::make_unique<Facet2d>( my_prolong_nodes[0]->get_coordinates(), my_prolong_nodes[1]->get_coordinates());
  }
  else
  {
    ThrowAssert(3 == my_prolong_nodes.size());
    my_facet = std::make_unique<Facet3d>( my_prolong_nodes[0]->get_coordinates(), my_prolong_nodes[1]->get_coordinates(), my_prolong_nodes[2]->get_coordinates());
  }
}

void ProlongationFacet::update_prolongation_point_data(const FacetDistanceQuery & dist_query) const
{
  ThrowAssert(&dist_query.facet() == my_facet.get());
  my_prolongation_point_data = std::make_unique<ProlongationPointData>(my_mesh, dist_query, my_prolong_nodes);
}

bool ProlongationFacet::communicate_me(const BoundingBox & proc_target_bbox) const
{
  for(auto && prolong_node : my_prolong_nodes)
  {
    if (!proc_target_bbox.contains(prolong_node->get_coordinates()))
    {
      return false;
    }
  }
  return true;
}

std::set<const ProlongationNodeData *>
ProlongationFacet::get_facet_nodes_to_communicate( const ProlongFacetVec & proc_prolong_facets, const BoundingBox & proc_target_bbox )
{
  std::set<const ProlongationNodeData *> procFacetNodes;

  for ( auto && facet : proc_prolong_facets )
    if (facet->communicate_me(proc_target_bbox))
      for (auto node : facet->my_prolong_nodes)
        procFacetNodes.insert(node);

  return procFacetNodes;
}

void
ProlongationFacet::communicate( const CDMesh & mesh, ProlongFacetVec & proc_prolong_facets, EntityProlongationNodeMap & proc_prolong_nodes, const std::vector<BoundingBox> & proc_target_bboxes )
{ /* %TRACE[ON]% */ Trace trace__("krino:ProlongationFacet::communicate()"); /* %TRACE% */
  communicate_shared_nodes(mesh, proc_prolong_nodes);
  communicate_facet_nodes(mesh, proc_prolong_facets, proc_prolong_nodes, proc_target_bboxes);
  communicate_facets(mesh, proc_prolong_facets, proc_target_bboxes);
}

static
void pack_facets_within_proc_bboxes(const stk::mesh::BulkData & mesh, const ProlongFacetVec & proc_prolong_facets, const std::vector<BoundingBox> & proc_target_bboxes, stk::CommSparse &commSparse)
{
  stk::pack_and_communicate(commSparse,[&]()
  {
    for ( int procId=0; procId<commSparse.parallel_size(); ++procId )
    {
      if ( commSparse.parallel_rank() == procId ) continue;  // Don't talk to yourself, it's embarrassing
      const BoundingBox & proc_target_bbox = proc_target_bboxes[procId];
      stk::CommBuffer & buffer = commSparse.send_buffer(procId);

      for ( auto && f_i : proc_prolong_facets )
      {
        const ProlongationFacet & facet = *f_i;
        if (facet.communicate_me(proc_target_bbox))
          facet.pack_into_buffer(buffer);
      }
    }
  });
}

static
void receive_and_build_prolong_facets(const CDMesh & mesh, ProlongFacetVec & proc_prolong_facets, stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      ProlongationFacet * facet = ProlongationFacet::unpack_from_buffer( mesh, buffer );
      proc_prolong_facets.push_back( facet );
    }
  });
}

void ProlongationFacet::communicate_facets( const CDMesh & mesh, ProlongFacetVec & proc_prolong_facets, const std::vector<BoundingBox> & proc_target_bboxes )
{
  stk::CommSparse commSparse(mesh.stk_bulk().parallel());
  pack_facets_within_proc_bboxes(mesh.stk_bulk(), proc_prolong_facets, proc_target_bboxes, commSparse);
  receive_and_build_prolong_facets(mesh, proc_prolong_facets, commSparse);
}

void
ProlongationFacet::pack_into_buffer(stk::CommBuffer & b) const
{
  const size_t num_prolong_nodes = my_prolong_nodes.size();
  b.pack(num_prolong_nodes);
  for(auto && prolong_node : my_prolong_nodes)
  {
    b.pack(prolong_node->entityId());
  }

  b.pack(my_common_fields.size());
  for (unsigned field : my_common_fields)
    b.pack(field);
}

ProlongationFacet *
ProlongationFacet::unpack_from_buffer(const CDMesh & mesh, stk::CommBuffer & b )
{
  size_t num_prolong_nodes = 0;
  b.unpack(num_prolong_nodes);
  std::vector<const ProlongationNodeData *> prolong_nodes(num_prolong_nodes);

  for(auto && prolong_node : prolong_nodes)
  {
    stk::mesh::EntityId node_id;
    b.unpack(node_id);
    prolong_node = mesh.fetch_prolong_node(node_id);
    ThrowRequireMsg(prolong_node, "Communication error, missing prolongation node " << node_id << " on processor " << mesh.stk_bulk().parallel_rank());
  }

  size_t num_common_fields = 0;
  b.unpack(num_common_fields);

  std::vector<unsigned> common_fields(num_common_fields);
  for (size_t ifield=0; ifield<num_common_fields; ++ifield)
    b.unpack(common_fields[ifield]);

  ProlongationFacet * facet = new ProlongationFacet(mesh, prolong_nodes, common_fields);

  return facet;
}

} // namespace krino
