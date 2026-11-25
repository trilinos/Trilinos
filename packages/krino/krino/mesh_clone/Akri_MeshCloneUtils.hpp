// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_MESH_CLONE_AKRI_MESHCLONEUTILS_HPP_
#define KRINO_KRINO_MESH_CLONE_AKRI_MESHCLONEUTILS_HPP_

namespace stk
{
namespace mesh
{
class BulkData;
class MetaData;
} // namespace mesh
} // namespace stk

namespace krino
{

// Clones the parts and fields of in_meta into out_meta, which must not be initialized yet, keeping
// the parts' and fields' ordinals, so that more parts and fields can then be declared on out_meta
// before it is committed.
void clone_meta_data_parts_and_fields(
    const stk::mesh::MetaData & in_meta, stk::mesh::MetaData & out_meta);

// Copies the IO attributes of in_meta's parts and fields (part attributes, distribution factor
// fields, field roles and output types) and its surface-to-block map to out_meta, a clone from
// clone_meta_data_parts_and_fields that has no fields of its own yet and is not committed.
void clone_io_attributes_and_surface_to_block_mapping(
    const stk::mesh::MetaData & in_meta, stk::mesh::MetaData & out_meta);

// Clones the locally owned and shared entities of in_mesh into out_mesh, whose MetaData is a
// committed clone of in_mesh's (possibly with more parts and fields), with their parts, relations,
// and node sharing, then copies the field data of in_mesh's fields.  Entity ids are kept.  out_mesh
// must have no entities and not be in a modification cycle.
void clone_bulk_data(const stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh);

// The steps of clone_bulk_data, for updating an existing clone.  clone_bulk_data_entities must be
// called within a modification cycle; if search_for_existing, entities already in out_mesh are kept
// and their parts are updated.  delete_extraneous_entities destroys the entities of out_mesh that
// are not in in_mesh, or whose owner or sharing differs.  copy_field_data copies the data of
// in_mesh's fields to the fields of out_mesh with the same ordinals.
void clone_bulk_data_entities(const stk::mesh::BulkData & in_mesh,
    stk::mesh::BulkData & out_mesh,
    const bool search_for_existing);
void delete_extraneous_entities(
    const stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh);
void copy_field_data(const stk::mesh::BulkData & in_mesh, stk::mesh::BulkData & out_mesh);

} // namespace krino

#endif /* KRINO_KRINO_MESH_CLONE_AKRI_MESHCLONEUTILS_HPP_ */
