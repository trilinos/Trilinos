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
#include <Akri_MeshCloneUtils.hpp>
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

  // The clone preserves the original entity identifiers, so resolve shared
  // side/edge sharing by entity id rather than by node connectivity alone.
  my_mesh->set_use_entity_ids_for_resolving_sharing(true);

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



} // namespace krino
