// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_tools/mesh_clone/ReplaceBulkData.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_io/IossBridge.hpp>

#include <Akri_MeshClone.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_MeshHelpers.hpp>

namespace krino{

bool
MeshClone::stash_or_restore_mesh(stk::mesh::BulkData & mesh, const unsigned step_count)
{ /* %TRACE[ON]% */ Trace trace__("krino::MeshClone::stash_or_restore_mesh(tftk::mesh::Mesh &mesh)"); /* %TRACE% */

  const auto empty_function = []() {};
  return stash_or_restore_mesh(mesh, step_count, empty_function, empty_function);
}

bool
MeshClone::stash_or_restore_mesh(stk::mesh::BulkData & mesh, const unsigned step_count,
    const std::function<void()> & notify_of_pre_mesh_modification,
    const std::function<void()> & notify_of_post_mesh_modification)
{ /* %TRACE[ON]% */ Trace trace__("krino::MeshClone::stash_or_restore_mesh(tftk::mesh::Mesh &mesh)"); /* %TRACE% */
  bool restored_mesh = false;
  MeshClone * clone = const_cast<MeshClone*>(mesh.mesh_meta_data().get_attribute<MeshClone>());
  if (!clone)
  {
    clone = new MeshClone(mesh, sierra::Diag::sierraTimer(), step_count);
    mesh.mesh_meta_data().declare_attribute_with_delete<MeshClone>(clone);
  }
  else
  {
    if (step_count > clone->my_step_count)
    {
      if (!clone->mesh_is_up_to_date())
      {
        clone->update(step_count);
      }
    }
    else if (!clone->mesh_is_up_to_date())
    {
      notify_of_pre_mesh_modification();
      clone->restore(step_count);
      notify_of_post_mesh_modification();
      restored_mesh = true;
    }
  }
  return restored_mesh;
}

bool
MeshClone::exists(const stk::mesh::BulkData & mesh)
{ /* %TRACE[ON]% */ Trace trace__("krino::MeshClone::exists(stk::mesh::BulkData & mesh)"); /* %TRACE% */
  MeshClone * clone = const_cast<MeshClone*>(mesh.mesh_meta_data().get_attribute<MeshClone>());
  return (clone != nullptr);
}

MeshClone &
MeshClone::get(const stk::mesh::BulkData & mesh)
{ /* %TRACE[ON]% */ Trace trace__("krino::MeshClone::update_cloned_fields(stk::mesh::BulkData & mesh)"); /* %TRACE% */
  MeshClone * clone = const_cast<MeshClone*>(mesh.mesh_meta_data().get_attribute<MeshClone>());
  STK_ThrowRequireMsg(clone != nullptr, "Could not find MeshClone.");
  return *clone;
}

MeshClone::MeshClone( stk::mesh::BulkData & orig_mesh, stk::diag::Timer parent_timer, const unsigned step_count )
: my_orig_mesh(&orig_mesh),
  my_timer("Clone mesh", parent_timer),
  my_step_count(step_count),
  my_synchronized_count(orig_mesh.synchronized_count())
{
  stk::diag::TimeBlock timer__(my_timer);
  const stk::mesh::MetaData & in_meta = my_orig_mesh->mesh_meta_data();
  my_meta = stk::mesh::MeshBuilder().create_meta_data();
  clone_meta_data_parts_and_fields(in_meta, *my_meta);

  my_mesh = stk::mesh::MeshBuilder(my_orig_mesh->parallel())
                       .set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA)
                       .set_add_fmwk_data(my_orig_mesh->add_fmwk_data())
                       .create(my_meta);

  my_meta->commit();

  my_mesh->modification_begin();
  clone_bulk_data_entities(*my_orig_mesh, *my_mesh, false);
  my_mesh->modification_end();

  copy_field_data(*my_orig_mesh, *my_mesh);
}

void MeshClone::update(const unsigned step_count)
{ /* %TRACE[ON]% */ Trace trace__("krino::MeshClone::update(const unsigned step_count)"); /* %TRACE% */
  stk::diag::TimeBlock timer__(my_timer);

  const stk::mesh::BulkData & in_mesh = *my_orig_mesh;
  stk::mesh::BulkData & out_mesh = *my_mesh;

  clone_mesh(in_mesh, out_mesh, false);

  my_step_count = step_count;
  mark_mesh_as_up_to_date();
}

void MeshClone::restore(const unsigned step_count)
{ /* %TRACE[ON]% */ Trace trace__("krino::MeshClone::restore(const unsigned step_count)"); /* %TRACE% */
  stk::diag::TimeBlock timer__(my_timer);
  krinolog << "Restoring mesh from clone." << stk::diag::dendl;

  const stk::mesh::BulkData & in_mesh = *my_mesh;
  stk::mesh::BulkData & out_mesh = *my_orig_mesh;

  clone_mesh(in_mesh, out_mesh, false);

  STK_ThrowRequire(step_count == my_step_count);
  mark_mesh_as_up_to_date();
}

void MeshClone::clone_mesh(const stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh, const bool full_overwrite)
{ /* %TRACE[ON]% */ Trace trace__("krino::MeshClone::clone_mesh(const stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh, const bool full_overwrite)"); /* %TRACE% */
  if (full_overwrite)
  {
//    std::function<void(stk::mesh::BulkData& outMesh_)> op = [](stk::mesh::BulkData& outMesh_) {};
    stk::tools::replace_bulk_data(in_mesh, out_mesh/*, op*/);
  }
  else
  {
    // I don't think you can start with clone_parallel because you can't know where something goes unless you know about it locally.
    // Therefore I don't think that creation can go before deletion because you can't know if it already exists.
    // Best to just delete everything wrong and then rebuild.
    STK_ThrowRequireMsg(out_mesh.modification_begin(), "MeshClone::restore must not be called within a modification cycle.");
    delete_extraneous_entities(in_mesh, out_mesh);
    // Occasionally, there is an issue with deleting and recreating the same entity within the same modification cycle,
    // so we need to end here and begin again for the creation.
    out_mesh.modification_end();

    //Delete graph because the following mod cycle creates elements and faces in the same mode cycle.
    out_mesh.delete_face_adjacent_element_graph();

    out_mesh.modification_begin();
    clone_bulk_data_entities(in_mesh, out_mesh, true);
    out_mesh.modification_end();
    copy_field_data(in_mesh, out_mesh);
  }

}


void
MeshClone::clone_meta_data_parts_and_fields(const stk::mesh::MetaData & in_meta, stk::mesh::MetaData & out_meta)
{
  /* %TRACE[ON]% */ Trace trace__("void MeshClone::clone_meta_data(stk::mesh::MetaData & in_meta)"); /* %TRACE% */

  // This is pretty nasty.  We want the part.mesh_meta_data_ordinal to be the same for the in_meta and out_meta.
  // To accomplish this, we must be careful about the order of the part creation.
  // Specifically, we must clone the parts that were created before in_meta.initialize() were called, then
  // call out_meta.initialize(), then clone the rest.

  const stk::mesh::PartVector & in_parts = in_meta.get_parts();

  unsigned ipart = 0;
  bool more_to_do = ipart < in_parts.size();
  while (more_to_do)
  {
    stk::mesh::Part * in_part = in_parts[ipart++];
    if (stk::mesh::is_topology_root_part(*in_part))
    {
      more_to_do = false;
    }
    else
    {
      more_to_do = ipart < in_parts.size();
      STK_ThrowRequire(in_part->primary_entity_rank() == stk::topology::INVALID_RANK);
      stk::mesh::Part & out_part = out_meta.declare_part(in_part->name());
      STK_ThrowRequire(out_part.mesh_meta_data_ordinal() == in_part->mesh_meta_data_ordinal());
    }
  }

  out_meta.initialize(in_meta.spatial_dimension(), in_meta.entity_rank_names());

  for ( auto&& in_part : in_parts )
  {
    stk::mesh::Part & out_part =
        (in_part->primary_entity_rank() == stk::topology::INVALID_RANK) ?
        out_meta.declare_part(in_part->name()) :
        out_meta.declare_part(in_part->name(), in_part->primary_entity_rank(), in_part->force_no_induce());
    STK_ThrowRequire(out_part.mesh_meta_data_ordinal() == in_part->mesh_meta_data_ordinal());
    if (stk::io::is_part_io_part(*in_part))
    {
      stk::io::put_io_part_attribute(out_part);
    }
  }

  for ( auto&& in_part : in_parts )
  {
    stk::mesh::Part & out_part = out_meta.get_part(in_part->mesh_meta_data_ordinal());
    const stk::mesh::PartVector & in_subsets = in_part->subsets();
    for (auto && in_subset : in_subsets)
    {
      stk::mesh::Part & out_subset = out_meta.get_part(in_subset->mesh_meta_data_ordinal());
      out_meta.declare_part_subset(out_part, out_subset);
    }
  }

  const stk::mesh::FieldVector & in_fields = in_meta.get_fields();
  for ( auto&& in_field : in_fields )
  {
    if (in_field->state() == stk::mesh::StateNone)
    {
      stk::mesh::FieldBase * out_field = in_field->clone(out_meta.get_field_repository());

      for ( auto&& in_restriction : in_field->restrictions() )
      {
        const stk::mesh::Selector out_selector = MeshClone::translate_selector(in_restriction.selector(), out_meta);

        out_meta.declare_field_restriction( *out_field, out_selector, in_restriction.num_scalars_per_entity(), in_restriction.dimension() );
      }
    }
  }
}

void
MeshClone::delete_extraneous_entities(const stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh)
{
  /* %TRACE[ON]% */ Trace trace__("void MeshClone::delete_extraneous_entities(stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh)"); /* %TRACE% */

  std::vector<stk::mesh::Entity> entities;
  std::vector<stk::mesh::Entity> relatives;
  std::vector<stk::mesh::ConnectivityOrdinal> relative_ordinals;
  std::vector<int> in_mesh_comm_procs;
  std::vector<int> out_mesh_comm_procs;

  const stk::mesh::EntityRank highest_entity_rank = static_cast<stk::mesh::EntityRank>(in_mesh.mesh_meta_data().entity_rank_count()-1);
  for (stk::mesh::EntityRank entity_rank = highest_entity_rank; entity_rank >= stk::topology::NODE_RANK; --entity_rank)
  {
    stk::mesh::get_entities( out_mesh, entity_rank, entities );

    for (auto&& out_entity : entities)
    {
      // Delete all entities for which any of the following are true:
      // 1. There is no corresponding entity in in_mesh
      // 2. The parallel owner changed.
      // 3. The comm procs or comm shared procs differ. (Ugh.  This is to work around strange corner case issues.)
      // 4. The entity has invalid nodes (due to the nodes violating 1 or 2).

      stk::mesh::Entity in_entity = get_entity_on_other_mesh(out_mesh, out_entity, in_mesh);
      bool delete_entity = !in_mesh.is_valid(in_entity);

      if (!delete_entity)
      {
        delete_entity =
            in_mesh.parallel_owner_rank(in_entity) != out_mesh.parallel_owner_rank(out_entity);
      }

      if (!delete_entity)
      {
        in_mesh.comm_procs(in_entity, in_mesh_comm_procs);
        out_mesh.comm_procs(out_entity, out_mesh_comm_procs);
        delete_entity = in_mesh_comm_procs != out_mesh_comm_procs;
      }

      if (!delete_entity)
      {
        in_mesh.comm_shared_procs(in_entity, in_mesh_comm_procs);
        out_mesh.comm_shared_procs(out_entity, out_mesh_comm_procs);
        delete_entity = in_mesh_comm_procs != out_mesh_comm_procs;
      }

      if (delete_entity)
      {
        STK_ThrowRequireMsg(disconnect_and_destroy_entity(out_mesh, out_entity), "Could not destroy entity " << out_mesh.entity_key(out_entity) << debug_entity(out_mesh, out_entity));
      }
    }
  }

  for (stk::mesh::EntityRank entity_rank = stk::topology::ELEMENT_RANK; entity_rank > stk::topology::NODE_RANK; --entity_rank)
  {
    stk::mesh::get_entities( out_mesh, entity_rank, entities );

    for (auto&& out_entity : entities)
    {
      const unsigned num_nodes = out_mesh.bucket(out_entity).topology().num_nodes();
      if (out_mesh.count_valid_connectivity(out_entity, stk::topology::NODE_RANK) != num_nodes)
      {
        STK_ThrowRequireMsg(disconnect_and_destroy_entity(out_mesh, out_entity), "Could not destroy entity " << out_mesh.entity_key(out_entity));
      }
    }
  }
}

void
MeshClone::delete_all_entities(stk::mesh::BulkData & mesh)
{
  /* %TRACE[ON]% */ Trace trace__("void MeshClone::delete_extraneous_entities(stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh)"); /* %TRACE% */

  std::vector<stk::mesh::Entity> entities;

  const stk::mesh::EntityRank highest_entity_rank = static_cast<stk::mesh::EntityRank>(mesh.mesh_meta_data().entity_rank_count()-1);
  for (stk::mesh::EntityRank entity_rank = highest_entity_rank; entity_rank >= stk::topology::NODE_RANK; --entity_rank)
  {
    stk::mesh::get_entities( mesh, entity_rank, entities );

    for (auto && entity : entities)
    {
      STK_ThrowRequireMsg(mesh.destroy_entity(entity), "Could not destroy entity " << mesh.entity_key(entity));
    }
  }
}

void
MeshClone::clone_bulk_data_entities(const stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh, const bool search_for_existing)
{
  /* %TRACE[ON]% */ Trace trace__("void MeshClone::clone_bulk_data(stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh)"); /* %TRACE% */

  const bool debug = false;

  const stk::mesh::MetaData & in_meta = in_mesh.mesh_meta_data();
  const stk::mesh::MetaData & out_meta = out_mesh.mesh_meta_data();

  std::vector<int> sharing;

  const stk::mesh::EntityRank highest_entity_rank = static_cast<stk::mesh::EntityRank>(in_meta.entity_rank_count()-1);
  for (stk::mesh::EntityRank entity_rank = stk::topology::NODE_RANK; entity_rank <= highest_entity_rank; ++entity_rank)
  {
    stk::mesh::Selector not_ghost = in_meta.locally_owned_part() | in_meta.globally_shared_part();
    const stk::mesh::BucketVector & buckets = in_mesh.get_buckets(entity_rank, not_ghost);

    for ( auto&& bucket_ptr : buckets )
    {
      const stk::mesh::Bucket & b = *bucket_ptr;

      stk::mesh::PartVector in_bucket_parts, out_bucket_parts;
      get_bucket_parts(b,in_bucket_parts);

      translate_parts(in_bucket_parts, out_meta, out_bucket_parts);

      const bool bucket_is_locally_owned = b.member(in_meta.locally_owned_part());


      const int length = b.size();
      for (int i = 0; i < length; ++i)
      {
        stk::mesh::Entity in_entity = b[i];

        stk::mesh::Entity out_entity;
        if (search_for_existing)
        {
          out_entity = out_mesh.get_entity( entity_rank, in_mesh.identifier(in_entity) );
          if (out_mesh.is_valid(out_entity) &&
              out_mesh.bucket(out_entity).member(out_meta.locally_owned_part()))
          {
            // make sure parts are right
            stk::mesh::PartVector current_parts;
            get_bucket_parts(out_mesh.bucket(out_entity), current_parts);
            out_mesh.change_entity_parts(out_entity, out_bucket_parts, current_parts);
          }
        }

        // Unfortunately, there is a subtle bug that may be present in the input mesh that can cause a face to be present on
        // this processor even when there are no locally owned elements using that face.  We will skip these faces.
        if (!bucket_is_locally_owned && entity_rank == in_meta.side_rank())
        {
          bool found_locally_owned_elem = false;
          const unsigned num_elems = in_mesh.num_elements(in_entity);
          const stk::mesh::Entity* elems = in_mesh.begin_elements(in_entity);
          for (unsigned ielem=0; ielem<num_elems; ++ielem)
          {
            if (in_mesh.bucket(elems[ielem]).member(in_meta.locally_owned_part()))
            {
              found_locally_owned_elem = true;
              break;
            }
          }
          if (!found_locally_owned_elem)
          {
            continue;
          }
        }

        if (!out_mesh.is_valid(out_entity))
        {
          out_entity = out_mesh.declare_entity( entity_rank, in_mesh.identifier(in_entity), out_bucket_parts );
        }

        if (stk::topology::NODE_RANK == entity_rank && out_mesh.parallel_size() > 1)
        {
          in_mesh.comm_shared_procs(in_entity,sharing);
          for ( size_t k = 0 ; k < sharing.size() ; ++k )
          {
            out_mesh.add_node_sharing(out_entity, sharing[k]);
          }
        }

        for (stk::mesh::EntityRank relative_rank = stk::topology::NODE_RANK; relative_rank < entity_rank; ++relative_rank)
        {
          const unsigned num_relatives = in_mesh.num_connectivity(in_entity, relative_rank);
          const stk::mesh::Entity* in_relatives = in_mesh.begin(in_entity, relative_rank);
          const stk::mesh::ConnectivityOrdinal * relative_ordinals = in_mesh.begin_ordinals(in_entity, relative_rank);
          const stk::mesh::Permutation * relative_permutations = in_mesh.begin_permutations(in_entity, relative_rank);

          for (unsigned relative_index=0; relative_index<num_relatives; ++relative_index)
          {
            stk::mesh::Entity in_relative = in_relatives[relative_index];

            stk::mesh::Entity out_relative = out_mesh.get_entity(relative_rank, in_mesh.identifier(in_relative));
            if (!out_mesh.is_valid(out_relative))
            {
              // This relative might live only on another processor
              continue;
            }

            if (in_mesh.has_permutation(in_entity, relative_rank))
            {
              out_mesh.declare_relation( out_entity, out_relative, relative_ordinals[relative_index], relative_permutations[relative_index] );
            }
            else
            {
              out_mesh.declare_relation( out_entity, out_relative, relative_ordinals[relative_index] );
            }
          }
        }
        if (debug)
        {
          krinolog << "Created missing entity "  << debug_entity(out_mesh, out_entity);
          krinolog << "By cloning entity      "  << debug_entity(in_mesh, in_entity);
        }
      }
    }
  }
}

void
MeshClone::copy_field_data(const stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh, const stk::mesh::FieldBase & in_field, const stk::mesh::FieldBase & out_field, const bool out_mesh_aura_from_communication)
{
  /* %TRACE[ON]% */ Trace trace__("void MeshClone::copy_field_data(stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh)"); /* %TRACE% */

  if ( !in_field.data_traits().is_pod )
  {
    return;
  }

  in_field.sync_to_host();
  out_field.sync_to_host();
  out_field.modify_on_host();

  stk::mesh::EntityRank entity_rank = out_field.entity_rank();
  STK_ThrowRequire(in_field.entity_rank() == entity_rank);

  stk::mesh::MetaData & out_meta = out_mesh.mesh_meta_data();
  stk::mesh::Selector field_selector =
      out_mesh_aura_from_communication ?
      (stk::mesh::selectField(out_field) & !out_meta.aura_part()) :
      stk::mesh::selectField(out_field);

  const stk::mesh::BucketVector & buckets = out_mesh.get_buckets(entity_rank, field_selector);

  for ( auto && bucket_ptr : buckets )
  {
    const stk::mesh::Bucket & b = *bucket_ptr;
    const unsigned bucket_length = b.size();
    const unsigned out_length = stk::mesh::field_bytes_per_entity(out_field, b);
    unsigned char * const out_bucket_data = reinterpret_cast<unsigned char *>( stk::mesh::field_data( out_field, b ) );

    for (unsigned ib=0; ib<bucket_length; ++ib)
    {
      stk::mesh::Entity out_entity = b[ib];
      stk::mesh::Entity in_entity = in_mesh.get_entity( entity_rank, out_mesh.identifier(out_entity) );
      const auto in_length = stk::mesh::field_bytes_per_entity(in_field, in_entity);
      STK_ThrowRequireMsg(in_mesh.is_valid(in_entity), "Missing entity " << out_mesh.entity_key(out_entity));
      STK_ThrowRequireMsg(in_length == out_length,
          "Mismatched field size for field " << in_field.name() << " in_length = " << in_length << " out_length = " << out_length << "\n"
          << " for input entity " << debug_entity(in_mesh, in_entity) << " on " << in_mesh.parallel_owner_rank(in_entity)
          << " and output entity " << debug_entity(out_mesh, out_entity) << " on " << out_mesh.parallel_owner_rank(out_entity) );

      unsigned char * const out_data = out_bucket_data + ib*out_length;
      unsigned char * const in_data = reinterpret_cast<unsigned char *>( stk::mesh::field_data( in_field, in_entity ) );
      for ( unsigned i = 0; i < out_length; i++ ) out_data[i] = in_data[i];
    }
  }
}

void
MeshClone::copy_field_data(const stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh)
{
  /* %TRACE[ON]% */ Trace trace__("void MeshClone::copy_field_data(stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh)"); /* %TRACE% */

  const stk::mesh::MetaData & in_meta = in_mesh.mesh_meta_data();
  stk::mesh::MetaData & out_meta = out_mesh.mesh_meta_data();

  const stk::mesh::FieldVector & in_fields = in_meta.get_fields();
  const stk::mesh::FieldVector & out_fields = out_meta.get_fields();
  STK_ThrowAssert(in_fields.size() == out_fields.size());

  const bool out_mesh_aura_from_communication = out_mesh.is_automatic_aura_on() && !in_mesh.is_automatic_aura_on();

  for ( unsigned field_index=0; field_index < in_fields.size(); ++field_index )
  {
    const stk::mesh::FieldBase & in_field = *in_fields[field_index];
    const stk::mesh::FieldBase & out_field = *out_fields[field_index];
    copy_field_data(in_mesh, out_mesh, in_field, out_field, out_mesh_aura_from_communication);
  }

  if (out_mesh_aura_from_communication)
  {
    const std::vector<stk::mesh::Ghosting *> ghostings = out_mesh.ghostings();
    const std::vector<const stk::mesh::FieldBase *> const_fields(out_fields.begin(), out_fields.end());
    stk::mesh::communicate_field_data(*ghostings[stk::mesh::BulkData::AURA], const_fields);
  }
}

void
MeshClone::translate_parts(const stk::mesh::PartVector & in_parts, const stk::mesh::MetaData & out_meta, stk::mesh::PartVector & out_parts)
{
  out_parts.clear();
  for ( auto&& in_part : in_parts )
  {
    stk::mesh::Part & out_part = out_meta.get_part(in_part->mesh_meta_data_ordinal());
    stk::mesh::insert(out_parts, out_part);
  }
}

stk::mesh::Selector
MeshClone::translate_selector(const stk::mesh::Selector & in_selector, const stk::mesh::MetaData & out_meta)
{
  if (in_selector == stk::mesh::Selector())
  {
    return in_selector;
  }
  return in_selector.clone_for_different_mesh(out_meta);
}

stk::mesh::Part *
MeshClone::translate_part(const stk::mesh::Part & in_part, const stk::mesh::MetaData & out_meta)
{
  return &out_meta.get_part(in_part.mesh_meta_data_ordinal());
}

void
MeshClone::get_bucket_parts(const stk::mesh::Bucket & bucket, stk::mesh::PartVector & parts)
{
  parts.clear();

  stk::mesh::PartVector const& bucket_parts = bucket.supersets();
  for ( stk::mesh::PartVector::const_iterator ip = bucket_parts.begin(); ip != bucket_parts.end(); ++ip )
  {
    stk::mesh::Part & bucket_part = **ip;
    if (bucket_part.primary_entity_rank() != stk::topology::INVALID_RANK && bucket_part.primary_entity_rank() != bucket.entity_rank())
    {
      continue;
    }
    if (stk::mesh::is_auto_declared_part(bucket_part) && !stk::mesh::is_topology_root_part(bucket_part))
    {
      continue;
    }

    parts.push_back(&bucket_part);
  }
  std::sort( parts.begin(), parts.end(), stk::mesh::PartLess() );
}

stk::mesh::Entity
MeshClone::get_entity_on_other_mesh(const stk::mesh::BulkData & mesh, stk::mesh::Entity entity, const stk::mesh::BulkData & other_mesh)
{
  stk::mesh::EntityRank entity_rank = mesh.entity_rank(entity);
  stk::mesh::Entity other_entity = other_mesh.get_entity( entity_rank, mesh.identifier(entity) );

  if (other_mesh.is_valid(other_entity) && entity_rank != stk::topology::NODE_RANK)
  {
    // check if nodes are the same
    const unsigned num_entity_nodes = mesh.num_nodes(entity);
    const unsigned num_other_entity_nodes = other_mesh.num_nodes(other_entity);
    if (num_entity_nodes != num_other_entity_nodes)
    {
      return stk::mesh::Entity();
    }
    const stk::mesh::Entity* entity_nodes = mesh.begin_nodes(entity);
    const stk::mesh::Entity* other_entity_nodes = other_mesh.begin_nodes(other_entity);
    for (unsigned n=0; n<num_entity_nodes; ++n)
    {
      stk::mesh::Entity node = entity_nodes[n];
      stk::mesh::Entity other_node = other_entity_nodes[n];
      if (!other_mesh.is_valid(other_node) || mesh.identifier(node) != other_mesh.identifier(other_node))
      {
        return stk::mesh::Entity();
      }
    }
  }

  return other_entity;
}

} // namespace krino
