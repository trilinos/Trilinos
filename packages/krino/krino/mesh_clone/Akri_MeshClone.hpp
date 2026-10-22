// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_MeshClone_h
#define Akri_MeshClone_h

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/diag/Timer.hpp>
#include <memory>

namespace krino {

class MeshClone {
public:

  MeshClone( stk::mesh::BulkData & orig_mesh, stk::diag::Timer parent_timer, const unsigned step_count = 0 );

  static bool exists(const stk::mesh::BulkData & mesh);
  static MeshClone & get(const stk::mesh::BulkData & mesh);
  static bool stash_or_restore_mesh(stk::mesh::BulkData & mesh, const unsigned step_count);
  static bool stash_or_restore_mesh(stk::mesh::BulkData & mesh,
    const unsigned step_count,
    const std::function<void()> & notify_of_pre_mesh_modification,
    const std::function<void()> & notify_of_post_mesh_modification);
  bool mesh_is_up_to_date() const {return my_orig_mesh->synchronized_count() == my_synchronized_count;}
  void mark_mesh_as_up_to_date() {my_synchronized_count = my_orig_mesh->synchronized_count();}
  void mark_mesh_as_out_of_date() {--my_synchronized_count;}

private:
  void update(const unsigned step_count = 0);
  void restore(const unsigned step_count = 0);

  static void clone_meta_data_parts_and_fields(const stk::mesh::MetaData & in_meta, stk::mesh::MetaData & out_meta);
  static void translate_parts(const stk::mesh::PartVector & in_parts, const stk::mesh::MetaData & out_meta, stk::mesh::PartVector & out_parts);
  static stk::mesh::Part * translate_part(const stk::mesh::Part & in_part, const stk::mesh::MetaData & out_meta);
  static stk::mesh::Selector translate_selector(const stk::mesh::Selector & in_selector, const stk::mesh::MetaData & out_meta);
  static void get_bucket_parts(const stk::mesh::Bucket & bucket, stk::mesh::PartVector & parts);

  static void clone_mesh(const stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh, const bool full_overwrite=false);
  static void clone_bulk_data(const stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh);
  static void copy_field_data(const stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh, const stk::mesh::FieldBase & in_field, const stk::mesh::FieldBase & out_field, const bool out_mesh_aura_from_communication);
  static void copy_field_data(const stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh);

  static void delete_extraneous_entities(const stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh);
  static void delete_all_entities(stk::mesh::BulkData & mesh);
  static void clone_bulk_data_entities(const stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh, const bool search_for_existing);

  static stk::mesh::Entity get_entity_on_other_mesh(const stk::mesh::BulkData & mesh, stk::mesh::Entity entity, const stk::mesh::BulkData & other_mesh);

  stk::mesh::BulkData* my_orig_mesh;
  std::shared_ptr<stk::mesh::MetaData> my_meta;
  std::unique_ptr<stk::mesh::BulkData> my_mesh;

  mutable stk::diag::Timer my_timer;
  unsigned my_step_count;
  size_t my_synchronized_count;
};

} // namespace krino

#endif // Akri_MeshClone_h
