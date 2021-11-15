// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Unit_Single_Element_Fixtures.hpp>

#include <gtest/gtest.h>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <Akri_CDFEM_Support.hpp>
#include <Akri_AuxMetaData.hpp>

namespace krino {

void SimpleStkFixture::write_results(const std::string & filename, stk::mesh::BulkData & mesh, const bool use64bitIds)
{
  stk::io::StkMeshIoBroker io(mesh.parallel());
  io.set_bulk_data(mesh);

  Ioss::PropertyManager properties;
  if (use64bitIds)
  {
    properties.add(Ioss::Property("INTEGER_SIZE_API", 8));
    properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
  }

  auto index = io.create_output_mesh(filename, stk::io::WRITE_RESULTS, properties);
  io.write_output_mesh(index);

  for(auto && field : mesh.mesh_meta_data().get_fields())
  {
    io.add_field(index, *field);
  }

  io.begin_output_step(index, 0.);
  io.write_defined_output_fields(index);
  io.end_output_step(index);
}

SingleElementFixture::SingleElementFixture(const stk::topology & topology)
  : my_topology(topology),
    stk_fixture(my_topology.dimension())
{
  AuxMetaData & aux_meta = AuxMetaData::get(stk_fixture.meta_data());
  block_part = &stk_fixture.meta_data().declare_part_with_topology("block_1", my_topology);
  const FieldType & vec_type = (my_topology.dimension() == 3) ? FieldType::VECTOR_3D : FieldType::VECTOR_2D;
  coord_field = aux_meta.register_field("coordinates", vec_type, stk::topology::NODE_RANK, 1u, 1u, *block_part);
  CDFEM_Support::get(stk_fixture.meta_data()).set_coords_field( coord_field );
  scalar_field = aux_meta.register_field("scalar_field", FieldType::REAL, stk::topology::NODE_RANK, 1u, 1u, *block_part);
  stk_fixture.commit();
}

void
SingleElementFixture::generate_mesh()
{
  const AuxMetaData & aux_meta = AuxMetaData::get(stk_fixture.meta_data());
  stk_fixture.bulk_data().modification_begin();
  {
    stk::mesh::PartVector elem_parts;
    elem_parts.push_back(block_part);
    elem_parts.push_back(&aux_meta.active_part());
    elem_parts.push_back(&stk_fixture.meta_data().locally_owned_part());
    stk::mesh::EntityIdVector node_ids(my_topology.num_nodes());
    for(unsigned i=0; i < node_ids.size(); ++i)
    {
      node_ids[i] = i+1;
    }
    const stk::mesh::EntityId elem_id = 1;
    my_elem = stk::mesh::declare_element( stk_fixture.bulk_data(), elem_parts, elem_id, node_ids );
    const stk::mesh::Entity * const nodes = stk_fixture.bulk_data().begin_nodes(my_elem);
    for(unsigned i=0; i < node_ids.size(); ++i)
    {
      EXPECT_EQ(node_ids[i], stk_fixture.bulk_data().identifier(nodes[i]));
    }
  }
  stk_fixture.bulk_data().modification_end();

  // set node coordinates
  stk::mesh::field_fill(0., coord_field);
  const stk::mesh::Entity * const nodes = stk_fixture.bulk_data().begin_nodes(my_elem);
  for(unsigned i=0; i < my_topology.num_nodes(); ++i)
  {
    if (i > 0)
    {
      double * node_coords = field_data<double>(coord_field, nodes[i]);
      node_coords[i-1] = 1.0;
    }
  }
}

TEST(SingleElementFixture, tri3)
{
  stk::topology tri3 = stk::topology::TRIANGLE_3_2D;
  SingleElementFixture test_fixture(tri3);

  test_fixture.generate_mesh();

  EXPECT_EQ(2u, test_fixture.stk_fixture.meta_data().spatial_dimension());

  const stk::mesh::BucketVector & elem_buckets = test_fixture.stk_fixture.bulk_data().buckets(stk::topology::ELEMENT_RANK);
  ASSERT_EQ(1u, elem_buckets.size());
  EXPECT_EQ(1u, (*elem_buckets.begin())->size());
  EXPECT_EQ(3u, test_fixture.stk_fixture.bulk_data().num_nodes((*elem_buckets[0])[0]));
  const stk::mesh::BucketVector & node_buckets = test_fixture.stk_fixture.bulk_data().buckets(stk::topology::NODE_RANK);
  ASSERT_EQ(1u, node_buckets.size());
  EXPECT_EQ(3u, (*node_buckets.begin())->size());
}

TEST(SingleElementFixture, tet4)
{
  stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  SingleElementFixture test_fixture(tet4);

  test_fixture.generate_mesh();

  EXPECT_EQ(3u, test_fixture.stk_fixture.meta_data().spatial_dimension());

  const stk::mesh::BucketVector & elem_buckets = test_fixture.stk_fixture.bulk_data().buckets(stk::topology::ELEMENT_RANK);
  ASSERT_EQ(1u, elem_buckets.size());
  EXPECT_EQ(1u, (*elem_buckets.begin())->size());
  EXPECT_EQ(4u, test_fixture.stk_fixture.bulk_data().num_nodes((*elem_buckets[0])[0]));
  const stk::mesh::BucketVector & node_buckets = test_fixture.stk_fixture.bulk_data().buckets(stk::topology::NODE_RANK);
  ASSERT_EQ(1u, node_buckets.size());
  EXPECT_EQ(4u, (*node_buckets.begin())->size());
}

}
