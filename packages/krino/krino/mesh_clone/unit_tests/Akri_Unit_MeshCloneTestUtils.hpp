// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_MESH_CLONE_UNIT_TESTS_AKRI_UNIT_MESHCLONETESTUTILS_HPP_
#define KRINO_KRINO_MESH_CLONE_UNIT_TESTS_AKRI_UNIT_MESHCLONETESTUTILS_HPP_

#include <memory>

#include <Akri_StkMeshFixture.hpp>

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

using TetMeshAndBuilder = StkMeshAndBuilder<stk::topology::TETRAHEDRON_4>;

// Four tets around the z axis, the first two in block_1 and the others in block_2, on processors
// 0,0,1,1 for 2 processors and one each for 4, with their exposed and block boundary sides. Fields:
// a two-state nodal field, an int element field on block_2 with an initial value, and a nodal field
// on subPart, a subset of superPart, which holds the center node.  Field values are set from the
// entity ids.  Returns false, building nothing, on more than 4 processors.
bool build_two_block_tet_mesh(TetMeshAndBuilder & fixture);

// A clone of mesh, by clone_meta_data_parts_and_fields and clone_bulk_data, for comparing mesh to
// later.
std::unique_ptr<stk::mesh::BulkData> clone_mesh(const stk::mesh::BulkData & mesh);

// in's parts and fields have the same ordinals in out, with the same names, ranks, topologies,
// subsets, states, data types, initial values and restrictions.  out may have more parts and
// fields.
void expect_equal_meta_data_parts_and_fields(
    const stk::mesh::MetaData & in, const stk::mesh::MetaData & out);

// in's fields have the same roles and output types in out, its parts the same distribution factor
// fields, and its surfaces in the surface to block map touch the same blocks, by name.
void expect_equal_io_attributes_and_surface_to_block_mapping(
    const stk::mesh::MetaData & in, const stk::mesh::MetaData & out);

// out has the same global number of entities of each rank as in, and each owned or shared entity of
// in is in out with the same owner, parts, sharing, downward connectivity (ids, ordinals and
// permutations), and data of in's fields.
void expect_equal_meshes(const stk::mesh::BulkData & in, const stk::mesh::BulkData & out);

} // namespace krino

#endif /* KRINO_KRINO_MESH_CLONE_UNIT_TESTS_AKRI_UNIT_MESHCLONETESTUTILS_HPP_ */
