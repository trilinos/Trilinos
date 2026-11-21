// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Unit_MeshCloneTestUtils.hpp>

#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/IossBridge.hpp>

#include <Akri_FieldRef.hpp>
#include <Akri_MeshCloneUtils.hpp>
#include <Akri_MeshSpecs.hpp>

namespace krino
{

bool build_two_block_tet_mesh(TetMeshAndBuilder & fixture)
{
  if (fixture.parallel_size() > 4) return false;

  stk::mesh::BulkData & mesh = fixture.get_mesh();
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();

  meta.declare_part_with_topology("block_1", stk::topology::TETRAHEDRON_4);
  stk::mesh::Part & block2 =
      meta.declare_part_with_topology("block_2", stk::topology::TETRAHEDRON_4);
  stk::mesh::Part & superPart = meta.declare_part("superPart", stk::topology::NODE_RANK);
  stk::mesh::Part & subPart = meta.declare_part("subPart", stk::topology::NODE_RANK);
  meta.declare_part_subset(superPart, subPart);

  auto & twoStateField = meta.declare_field<double>(stk::topology::NODE_RANK, "two_state", 2);
  stk::mesh::put_field_on_mesh(twoStateField, meta.universal_part(), nullptr);
  auto & elemIntField = meta.declare_field<int>(stk::topology::ELEMENT_RANK, "elem_int");
  const int elemIntInitialValue = 7;
  stk::mesh::put_field_on_mesh(elemIntField, block2, &elemIntInitialValue);
  auto & subsetField = meta.declare_field<double>(stk::topology::NODE_RANK, "subset_field");
  stk::mesh::put_field_on_mesh(subsetField, subPart, nullptr);

  fixture.build_mesh_from_spec_in_blocks(FourRightTets(), {1, 1, 2, 2}, {{2, {0, 0, 1, 1}}});

  std::vector<stk::mesh::Entity> nodesToAdd;
  const stk::mesh::Entity centerNode = fixture.get_assigned_node_for_index(0);
  if (fixture.is_owned_node(0)) nodesToAdd.push_back(centerNode);
  mesh.batch_change_entity_parts(nodesToAdd, {&subPart}, {});

  stk::mesh::FieldBase & oldStateField = *twoStateField.field_state(stk::mesh::StateOld);
  for (auto * bucketPtr : mesh.buckets(stk::topology::NODE_RANK))
  {
    for (auto node : *bucketPtr)
    {
      const double id = static_cast<double>(mesh.identifier(node));
      *field_data<double>(twoStateField, node) = id;
      *field_data<double>(oldStateField, node) = 10. * id;
      if (bucketPtr->member(subPart)) *field_data<double>(subsetField, node) = id + 0.5;
    }
  }
  for (auto * bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, block2))
    for (auto elem : *bucketPtr)
      *field_data<int>(elemIntField, elem) = static_cast<int>(mesh.identifier(elem));

  return true;
}

std::unique_ptr<stk::mesh::BulkData> clone_mesh(const stk::mesh::BulkData & mesh)
{
  std::shared_ptr<stk::mesh::MetaData> meta =
      stk::mesh::MeshBuilder(mesh.parallel()).create_meta_data();
  clone_meta_data_parts_and_fields(mesh.mesh_meta_data(), *meta);
  std::unique_ptr<stk::mesh::BulkData> clone = stk::mesh::MeshBuilder(mesh.parallel()).create(meta);
  meta->commit();
  clone_bulk_data(mesh, *clone);
  return clone;
}

static std::string to_string(const stk::mesh::Selector & selector)
{
  std::ostringstream os;
  os << selector;
  return os.str();
}

// The ordinals of the subsets of part, of those below numParts.
static std::vector<stk::mesh::PartOrdinal>
subset_ordinals(const stk::mesh::Part & part, const size_t numParts)
{
  std::vector<stk::mesh::PartOrdinal> ordinals;
  for (auto * subset : part.subsets())
    if (subset->mesh_meta_data_ordinal() < numParts)
      ordinals.push_back(subset->mesh_meta_data_ordinal());
  return ordinals;
}

void expect_equal_meta_data_parts_and_fields(
    const stk::mesh::MetaData & in, const stk::mesh::MetaData & out)
{
  EXPECT_EQ(in.spatial_dimension(), out.spatial_dimension());
  EXPECT_EQ(in.coordinate_field_name(), out.coordinate_field_name());

  ASSERT_LE(in.get_parts().size(), out.get_parts().size());
  for (auto * inPart : in.get_parts())
  {
    const stk::mesh::Part & outPart = *out.get_parts()[inPart->mesh_meta_data_ordinal()];
    EXPECT_EQ(inPart->name(), outPart.name());
    EXPECT_EQ(inPart->primary_entity_rank(), outPart.primary_entity_rank()) << inPart->name();
    EXPECT_EQ(inPart->topology(), outPart.topology()) << inPart->name();
    EXPECT_EQ(inPart->force_no_induce(), outPart.force_no_induce()) << inPart->name();
    EXPECT_EQ(stk::io::is_part_io_part(*inPart), stk::io::is_part_io_part(outPart))
        << inPart->name();
    const size_t numInParts = in.get_parts().size();
    EXPECT_EQ(subset_ordinals(*inPart, numInParts), subset_ordinals(outPart, numInParts))
        << inPart->name();
  }

  ASSERT_LE(in.get_fields().size(), out.get_fields().size());
  for (auto * inField : in.get_fields())
  {
    const stk::mesh::FieldBase & outField = *out.get_fields()[inField->mesh_meta_data_ordinal()];
    EXPECT_EQ(inField->name(), outField.name());
    EXPECT_EQ(inField->entity_rank(), outField.entity_rank()) << inField->name();
    EXPECT_EQ(inField->number_of_states(), outField.number_of_states()) << inField->name();
    EXPECT_EQ(inField->state(), outField.state()) << inField->name();
    EXPECT_EQ(inField->data_traits().type_info, outField.data_traits().type_info)
        << inField->name();

    const auto & inInitialValue = inField->get_initial_value_bytes();
    const auto & outInitialValue = outField.get_initial_value_bytes();
    ASSERT_EQ(inInitialValue.extent(0), outInitialValue.extent(0)) << inField->name();
    for (size_t i = 0; i < inInitialValue.extent(0); ++i)
      EXPECT_EQ(inInitialValue.data()[i], outInitialValue.data()[i]) << inField->name();

    const stk::mesh::FieldRestrictionVector & inRestrictions = inField->restrictions();
    const stk::mesh::FieldRestrictionVector & outRestrictions = outField.restrictions();
    ASSERT_EQ(inRestrictions.size(), outRestrictions.size()) << inField->name();
    for (size_t i = 0; i < inRestrictions.size(); ++i)
    {
      EXPECT_EQ(to_string(inRestrictions[i].selector()), to_string(outRestrictions[i].selector()))
          << inField->name();
      EXPECT_EQ(
          inRestrictions[i].num_scalars_per_entity(), outRestrictions[i].num_scalars_per_entity())
          << inField->name();
      EXPECT_EQ(inRestrictions[i].dimension(), outRestrictions[i].dimension()) << inField->name();
    }
  }
}

static int field_ordinal_or_minus_one(const stk::mesh::FieldBase * field)
{
  return (field == nullptr) ? -1 : static_cast<int>(field->mesh_meta_data_ordinal());
}

static std::set<std::string>
names_of_blocks_touching_surface(const stk::mesh::MetaData & meta, const stk::mesh::Part & surface)
{
  std::set<std::string> blockNames;
  for (auto * block : meta.get_blocks_touching_surface(&surface))
    blockNames.insert(block->name());
  return blockNames;
}

void expect_equal_io_attributes_and_surface_to_block_mapping(
    const stk::mesh::MetaData & in, const stk::mesh::MetaData & out)
{
  for (auto * inField : in.get_fields())
  {
    const stk::mesh::FieldBase & outField = *out.get_fields()[inField->mesh_meta_data_ordinal()];
    const Ioss::Field::RoleType * inRole = stk::io::get_field_role(*inField);
    const Ioss::Field::RoleType * outRole = stk::io::get_field_role(outField);
    ASSERT_EQ(inRole == nullptr, outRole == nullptr) << inField->name();
    if (inRole != nullptr)
    {
      EXPECT_EQ(*inRole, *outRole) << inField->name();
    }
    ASSERT_EQ(stk::io::has_field_output_type(*inField), stk::io::has_field_output_type(outField))
        << inField->name();
    if (stk::io::has_field_output_type(*inField))
    {
      EXPECT_EQ(stk::io::get_field_output_type(*inField), stk::io::get_field_output_type(outField))
          << inField->name();
    }
  }

  for (auto * inPart : in.get_parts())
  {
    const stk::mesh::Part & outPart = *out.get_parts()[inPart->mesh_meta_data_ordinal()];
    EXPECT_EQ(field_ordinal_or_minus_one(stk::io::get_distribution_factor_field(*inPart)),
        field_ordinal_or_minus_one(stk::io::get_distribution_factor_field(outPart)))
        << inPart->name();
  }

  const std::vector<const stk::mesh::Part *> inSurfaces = in.get_surfaces_in_surface_to_block_map();
  EXPECT_EQ(inSurfaces.size(), out.get_surfaces_in_surface_to_block_map().size());
  for (auto * inSurface : inSurfaces)
  {
    const stk::mesh::Part & outSurface = *out.get_parts()[inSurface->mesh_meta_data_ordinal()];
    EXPECT_EQ(names_of_blocks_touching_surface(in, *inSurface),
        names_of_blocks_touching_surface(out, outSurface))
        << inSurface->name();
  }
}

static std::vector<stk::mesh::PartOrdinal>
part_ordinals(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity)
{
  std::vector<stk::mesh::PartOrdinal> ordinals;
  for (auto * part : mesh.bucket(entity).supersets())
    ordinals.push_back(part->mesh_meta_data_ordinal());
  return ordinals;
}

static void expect_equal_downward_connectivity(const stk::mesh::BulkData & in,
    const stk::mesh::Entity inEntity,
    const stk::mesh::BulkData & out,
    const stk::mesh::Entity outEntity)
{
  for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < in.entity_rank(inEntity);
      ++rank)
  {
    const unsigned numConnected = in.num_connectivity(inEntity, rank);
    ASSERT_EQ(numConnected, out.num_connectivity(outEntity, rank))
        << in.entity_key(inEntity) << " rank " << rank;
    const stk::mesh::Entity * inConnected = in.begin(inEntity, rank);
    const stk::mesh::Entity * outConnected = out.begin(outEntity, rank);
    const stk::mesh::ConnectivityOrdinal * inOrdinals = in.begin_ordinals(inEntity, rank);
    const stk::mesh::ConnectivityOrdinal * outOrdinals = out.begin_ordinals(outEntity, rank);
    const stk::mesh::Permutation * inPermutations = in.begin_permutations(inEntity, rank);
    const stk::mesh::Permutation * outPermutations = out.begin_permutations(outEntity, rank);
    EXPECT_EQ(inPermutations == nullptr, outPermutations == nullptr)
        << in.entity_key(inEntity) << " rank " << rank;
    for (unsigned i = 0; i < numConnected; ++i)
    {
      EXPECT_EQ(in.identifier(inConnected[i]), out.identifier(outConnected[i]))
          << in.entity_key(inEntity);
      EXPECT_EQ(inOrdinals[i], outOrdinals[i]) << in.entity_key(inEntity);
      if (inPermutations != nullptr && outPermutations != nullptr)
      {
        EXPECT_EQ(inPermutations[i], outPermutations[i]) << in.entity_key(inEntity);
      }
    }
  }
}

static void expect_equal_field_data(const stk::mesh::BulkData & in,
    const stk::mesh::Entity inEntity,
    const stk::mesh::BulkData & out,
    const stk::mesh::Entity outEntity)
{
  const stk::mesh::FieldVector & outFields = out.mesh_meta_data().get_fields();
  for (auto * inField : in.mesh_meta_data().get_fields(in.entity_rank(inEntity)))
  {
    if (!inField->data_traits().is_pod) continue;
    const stk::mesh::FieldBase & outField = *outFields[inField->mesh_meta_data_ordinal()];
    auto inBytes = inField->data_bytes<const std::byte>().entity_bytes(inEntity);
    auto outBytes = outField.data_bytes<const std::byte>().entity_bytes(outEntity);
    ASSERT_EQ(inBytes.num_bytes(), outBytes.num_bytes())
        << inField->name() << " " << in.entity_key(inEntity);
    for (stk::mesh::ByteIdx i : inBytes.bytes())
      EXPECT_EQ(inBytes(i), outBytes(i)) << inField->name() << " " << in.entity_key(inEntity);
  }
}

void expect_equal_meshes(const stk::mesh::BulkData & in, const stk::mesh::BulkData & out)
{
  std::vector<size_t> inCounts, outCounts;
  stk::mesh::comm_mesh_counts(in, inCounts);
  stk::mesh::comm_mesh_counts(out, outCounts);
  EXPECT_EQ(inCounts, outCounts);

  const stk::mesh::MetaData & inMeta = in.mesh_meta_data();
  const stk::mesh::Selector ownedOrShared =
      inMeta.locally_owned_part() | inMeta.globally_shared_part();
  std::vector<int> inSharing, outSharing;
  for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < inMeta.entity_rank_count();
      ++rank)
  {
    for (auto * bucketPtr : in.get_buckets(rank, ownedOrShared))
    {
      for (auto inEntity : *bucketPtr)
      {
        const stk::mesh::Entity outEntity = out.get_entity(in.entity_key(inEntity));
        ASSERT_TRUE(out.is_valid(outEntity)) << in.entity_key(inEntity);
        EXPECT_EQ(in.parallel_owner_rank(inEntity), out.parallel_owner_rank(outEntity))
            << in.entity_key(inEntity);
        EXPECT_EQ(part_ordinals(in, inEntity), part_ordinals(out, outEntity))
            << in.entity_key(inEntity);
        in.comm_shared_procs(inEntity, inSharing);
        out.comm_shared_procs(outEntity, outSharing);
        EXPECT_EQ(inSharing, outSharing) << in.entity_key(inEntity);
        expect_equal_downward_connectivity(in, inEntity, out, outEntity);
        expect_equal_field_data(in, inEntity, out, outEntity);
      }
    }
  }
}

} // namespace krino
