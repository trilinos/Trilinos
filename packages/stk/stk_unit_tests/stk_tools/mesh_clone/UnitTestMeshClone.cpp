// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_tools/mesh_clone/MeshClone.hpp>
#include <stk_tools/mesh_clone/MeshCloneUtils.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>

namespace
{
using stk::unit_test_util::build_mesh;

class NoDeleteAttribute {
public:
  NoDeleteAttribute()
    : m_value(123)
  {}
  int value() const { return m_value; }

private:
  int m_value;
};

class WithDeleteAttribute {
public:
  WithDeleteAttribute()
    : m_value(456)
  {}
  int value() const { return m_value; }

private:
  int m_value;
};

void expect_equal_parts(const stk::mesh::Part &oldPart, const stk::mesh::Part &newPart)
{
  EXPECT_EQ(oldPart.name(), newPart.name());
  EXPECT_EQ(oldPart.id(), newPart.id()) << oldPart.name();
  EXPECT_EQ(oldPart.topology(), newPart.topology()) << oldPart.name();
  EXPECT_EQ(oldPart.primary_entity_rank(), newPart.primary_entity_rank()) << oldPart.name();
  EXPECT_EQ(oldPart.force_no_induce(), newPart.force_no_induce());

  EXPECT_EQ(stk::io::is_part_io_part(oldPart), stk::io::is_part_io_part(newPart));
  if (oldPart.attribute<NoDeleteAttribute>()) {
    EXPECT_EQ(oldPart.attribute<NoDeleteAttribute>()->value(), newPart.attribute<NoDeleteAttribute>()->value());
  }
  if (oldPart.attribute<WithDeleteAttribute>()) {
    EXPECT_EQ(oldPart.attribute<WithDeleteAttribute>()->value(), newPart.attribute<WithDeleteAttribute>()->value());
  }

  const stk::mesh::FieldBase *oldDistFactField = stk::io::get_distribution_factor_field(oldPart);
  if (oldDistFactField != nullptr) {
    const stk::mesh::FieldBase *newDistFactField = stk::io::get_distribution_factor_field(newPart);
    ASSERT_TRUE(newDistFactField != nullptr);
    EXPECT_EQ(oldDistFactField->mesh_meta_data_ordinal(), newDistFactField->mesh_meta_data_ordinal());
    EXPECT_NE(oldDistFactField, newDistFactField) << "Distribution factor fields must be different for different meshes";
  }

  stk::mesh::OrdinalVector oldSupersets = stk::tools::get_part_supersets(oldPart);
  stk::mesh::OrdinalVector newSupersets = stk::tools::get_part_supersets(newPart);
  EXPECT_TRUE(oldSupersets == newSupersets) << oldPart.name();
}

void expect_all_parts_equal(stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
  const stk::mesh::PartVector &oldParts = oldMeta.get_parts();
  const stk::mesh::PartVector &newParts = newMeta.get_parts();
  ASSERT_EQ(oldParts.size(), newParts.size());
  for(size_t i=0; i<oldParts.size(); i++)
    expect_equal_parts(*oldParts[i], *newParts[i]);
}

void expect_part_pointers_in_new_meta(stk::mesh::MetaData &newMeta, const stk::mesh::FieldRestriction &newRestriction)
{
  stk::mesh::PartVector selectorParts;
  newRestriction.selector().get_parts(selectorParts);
  for(stk::mesh::Part *part : selectorParts)
  {
    stk::mesh::Part &metaPart = newMeta.get_part(part->mesh_meta_data_ordinal());
    EXPECT_EQ(&metaPart, part);
  }
}

void expect_equal_selectors(stk::mesh::Selector oldSel, stk::mesh::Selector newSel)
{
  std::ostringstream oldS;
  oldS << oldSel;
  std::ostringstream newS;
  newS << newSel;
  EXPECT_EQ(oldS.str(), newS.str());
}

void expect_equal_field_restrictions(stk::mesh::MetaData &newMeta, stk::mesh::FieldBase &oldField, stk::mesh::FieldBase &newField)
{
  const stk::mesh::FieldRestrictionVector &oldRestrictions = oldField.restrictions();
  const stk::mesh::FieldRestrictionVector &newRestrictions = newField.restrictions();
  ASSERT_EQ(oldRestrictions.size(), newRestrictions.size());
  for(unsigned j=0; j<oldRestrictions.size(); j++)
  {
    expect_equal_selectors(oldRestrictions[j].selector(), newRestrictions[j].selector());
    expect_part_pointers_in_new_meta(newMeta, newRestrictions[j]);

    EXPECT_EQ(oldRestrictions[j].num_scalars_per_entity(), newRestrictions[j].num_scalars_per_entity());
    EXPECT_EQ(oldRestrictions[j].dimension(), newRestrictions[j].dimension());
  }
}

void expect_equal_field_initial_values(stk::mesh::FieldBase &oldField, stk::mesh::FieldBase &newField)
{
  ASSERT_EQ(oldField.get_initial_value_num_bytes(), newField.get_initial_value_num_bytes());
  const char *oldInitValPtr = static_cast<const char *>(oldField.get_initial_value());
  const char *newInitValPtr = static_cast<const char *>(newField.get_initial_value());
  for(unsigned j=0; j<oldField.get_initial_value_num_bytes(); j++)
    EXPECT_EQ(oldInitValPtr[j], newInitValPtr[j]);
}

void expect_equal_data_traits(const stk::mesh::DataTraits &oldTraits, const stk::mesh::DataTraits &newTraits)
{
  EXPECT_EQ(oldTraits.type_info, newTraits.type_info);
  EXPECT_EQ(oldTraits.size_of, newTraits.size_of);
  EXPECT_EQ(oldTraits.is_void, newTraits.is_void);
  EXPECT_EQ(oldTraits.is_integral, newTraits.is_integral);
  EXPECT_EQ(oldTraits.is_floating_point, newTraits.is_floating_point);
  EXPECT_EQ(oldTraits.is_array, newTraits.is_array);
  EXPECT_EQ(oldTraits.is_pointer, newTraits.is_pointer);
  EXPECT_EQ(oldTraits.is_enum, newTraits.is_enum);
  EXPECT_EQ(oldTraits.is_class, newTraits.is_class);
  EXPECT_EQ(oldTraits.is_pod, newTraits.is_pod);
  EXPECT_EQ(oldTraits.is_signed, newTraits.is_signed);
  EXPECT_EQ(oldTraits.is_unsigned, newTraits.is_unsigned);
  EXPECT_EQ(oldTraits.alignment_of, newTraits.alignment_of);
  EXPECT_EQ(oldTraits.stride_of, newTraits.stride_of);
  EXPECT_EQ(oldTraits.name, newTraits.name);
}

void expect_equal_fields(stk::mesh::MetaData &newMeta, stk::mesh::FieldBase &oldField, stk::mesh::FieldBase &newField)
{
  EXPECT_EQ(oldField.mesh_meta_data_ordinal(), newField.mesh_meta_data_ordinal());
  EXPECT_EQ(oldField.name(), newField.name());
  EXPECT_EQ(oldField.entity_rank(), newField.entity_rank());
  EXPECT_EQ(oldField.number_of_states(), newField.number_of_states());
  EXPECT_EQ(oldField.state(), newField.state());
  EXPECT_EQ(oldField.type_is<int>(), newField.type_is<int>());
  EXPECT_EQ(oldField.type_is<double>(), newField.type_is<double>());

  if(oldField.attribute<Ioss::Field::RoleType>() != nullptr)
  {
    EXPECT_EQ(*oldField.attribute<Ioss::Field::RoleType>(), *newField.attribute<Ioss::Field::RoleType>());
  }

  for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<oldField.entity_rank(); rank++)
  {
    EXPECT_EQ(oldField.max_size(), newField.max_size());
  }

  ASSERT_EQ(stk::io::has_field_output_type(oldField), stk::io::has_field_output_type(newField));
  if (stk::io::has_field_output_type(oldField)) {
    EXPECT_EQ(stk::io::get_field_output_type(oldField), stk::io::get_field_output_type(newField));
  }

  expect_equal_data_traits(oldField.data_traits(), newField.data_traits());
  expect_equal_field_initial_values(oldField, newField);
  expect_equal_field_restrictions(newMeta, oldField, newField);
}

void expect_all_fields_equal(stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
  const stk::mesh::FieldVector &oldFields = oldMeta.get_fields();
  const stk::mesh::FieldVector &newFields = newMeta.get_fields();
  ASSERT_EQ(oldFields.size(), newFields.size());
  for(size_t i=0; i<oldFields.size(); i++) {
    expect_equal_fields(newMeta, *oldFields[i], *newFields[i]);
  }
}

void expect_surface_to_block_mappings_equal(stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
  std::vector<const stk::mesh::Part *> oldSurfacesInMap = oldMeta.get_surfaces_in_surface_to_block_map();
  std::vector<const stk::mesh::Part *> newSurfacesInMap = newMeta.get_surfaces_in_surface_to_block_map();

  ASSERT_EQ(oldSurfacesInMap.size(), newSurfacesInMap.size());

  for(size_t i=0;i<oldSurfacesInMap.size();++i)
  {
    std::vector<const stk::mesh::Part*> oldBlocks = oldMeta.get_blocks_touching_surface(oldSurfacesInMap[i]);
    std::vector<const stk::mesh::Part*> newBlocks = newMeta.get_blocks_touching_surface(newSurfacesInMap[i]);
    ASSERT_EQ(oldBlocks.size(), newBlocks.size());
    for(size_t j=0; j<oldBlocks.size(); ++j)
    {
      EXPECT_EQ(oldBlocks[j]->name(), newBlocks[j]->name());
    }
  }
}

void expect_equal_meta_datas(stk::mesh::MetaData& oldMeta, stk::mesh::MetaData& newMeta)
{
  EXPECT_TRUE(newMeta.is_initialized());
  EXPECT_EQ(oldMeta.spatial_dimension(), newMeta.spatial_dimension());
  EXPECT_EQ(oldMeta.entity_rank_count(), newMeta.entity_rank_count());
  EXPECT_EQ(oldMeta.entity_rank_names(), newMeta.entity_rank_names());
  EXPECT_EQ(oldMeta.coordinate_field()->name(), newMeta.coordinate_field()->name());
  expect_all_parts_equal(oldMeta, newMeta);
  expect_all_fields_equal(oldMeta, newMeta);
  expect_surface_to_block_mappings_equal(oldMeta, newMeta);
}

void expect_equal_entity_counts(stk::mesh::BulkData& oldBulk, stk::mesh::BulkData& newBulk)
{
  std::vector<unsigned> oldCount;
  std::vector<unsigned> newCount;

  const size_t oldNranks = oldBulk.mesh_meta_data().entity_rank_count();
  const size_t newNranks = newBulk.mesh_meta_data().entity_rank_count();

  EXPECT_EQ(oldNranks, newNranks);

  oldCount.resize( oldNranks );
  newCount.resize( newNranks );

  for ( size_t i = 0 ; i < oldNranks ; ++i )
  {
    oldCount[i] = stk::mesh::count_selected_entities(oldBulk.mesh_meta_data().locally_owned_part(), oldBulk.buckets( static_cast<stk::mesh::EntityRank>(i) ));
    newCount[i] = stk::mesh::count_selected_entities(newBulk.mesh_meta_data().locally_owned_part(), newBulk.buckets( static_cast<stk::mesh::EntityRank>(i) ));
  }

  EXPECT_EQ(oldCount, newCount);
}

class MeshClone : public stk::unit_test_util::MeshFixture
{
public:
  MeshClone()
    : m_coordsField(nullptr),
      m_noDeleteAttribute()
  {}

protected:
  const char *get_mesh_spec() const {return "generated:1x1x8|sideset:x";}
  stk::mesh::BulkData::AutomaticAuraOption get_aura_option() const {return stk::mesh::BulkData::NO_AUTO_AURA;}

  void initialize_mesh_with_parts_and_fields()
  {
    setup_empty_mesh(get_aura_option());
    setup_parts();
    setup_fields();
    stk::io::fill_mesh(get_mesh_spec(), get_bulk());
    get_meta().set_coordinate_field(m_coordsField);
    add_distribution_factor();
    add_part_attributes({"block_1", "surface_1"});
  }

  void initialize_multiproc_multiblock_mesh()
  {
    const std::string meshSpec = "textmesh:0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                                 "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2\n";
    setup_mesh(meshSpec, stk::mesh::BulkData::AUTO_AURA);
  }

  void add_part_attributes(const std::vector<std::string> & partNames)
  {
    for (const std::string & partName : partNames) {
      stk::mesh::Part * part = get_meta().get_part(partName);
      STK_ThrowRequire(part != nullptr);
      get_meta().declare_attribute_no_delete(*part, &m_noDeleteAttribute);
      get_meta().declare_attribute_with_delete(*part, new WithDeleteAttribute);
    }
  }

  void add_distribution_factor() {
    stk::io::set_distribution_factor_field(*m_block1Part, *m_distFactField);
  }

  void setup_parts()
  {
    m_superPart = &get_meta().declare_part("superPart");
    m_subPart = &get_meta().declare_part("subPart");
    get_meta().declare_part_subset(*m_superPart, *m_subPart);

    m_block1Part = &get_meta().declare_part("block_1");
    m_surface1Part = &get_meta().declare_part("surface_1");
  }

  void setup_fields()
  {
    m_coordsField = &get_meta().declare_field<double>(stk::topology::NODE_RANK, "Coordinates", 1);
    m_surfaceField = &get_meta().declare_field<double>(stk::topology::NODE_RANK, "surfaceField", 2);
    m_distFactField = &get_meta().declare_field<double>(stk::topology::NODE_RANK, "dist_fact", 1);

    const double vectInitValue[] = {13, 14, 15};
    const double scalarInitValue = 3;
    stk::mesh::put_field_on_mesh(*m_coordsField, *m_block1Part, 3, vectInitValue);
    stk::mesh::put_field_on_mesh(*m_surfaceField, *m_surface1Part, 3, vectInitValue);
    stk::mesh::put_field_on_mesh(*m_distFactField, *m_block1Part, 1, &scalarInitValue);
    stk::io::set_field_output_type(*m_coordsField, stk::io::FieldOutputType::VECTOR_3D);
    stk::io::set_field_output_type(*m_surfaceField, stk::io::FieldOutputType::VECTOR_3D);
  }

  void add_orphan_nodes(const unsigned numOrphansPerProc)
  {
    std::vector<stk::mesh::EntityId> newIds;
    get_bulk().generate_new_ids(stk::topology::NODE_RANK, numOrphansPerProc, newIds);

    get_bulk().modification_begin();
    for(unsigned i=0; i<numOrphansPerProc; i++)
    {
      stk::mesh::EntityId id = newIds[i];
      get_bulk().declare_entity(stk::topology::NODE_RANK, id, *m_block1Part);
    }
    get_bulk().modification_end();
  }

  stk::mesh::Part* m_superPart;
  stk::mesh::Part* m_subPart;
  stk::mesh::Part* m_block1Part;
  stk::mesh::Part* m_surface1Part;
  stk::mesh::Field<double>* m_coordsField;
  stk::mesh::Field<double>* m_surfaceField;
  stk::mesh::Field<double>* m_distFactField;
  NoDeleteAttribute m_noDeleteAttribute;
};

TEST_F(MeshClone, ifSameMesh_throws)
{
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  ASSERT_THROW(stk::tools::copy_mesh(get_bulk(), get_meta().universal_part(), get_bulk()), std::exception);
}

TEST_F(MeshClone, copyOnlyMeta)
{
  initialize_mesh_with_parts_and_fields();

  stk::mesh::MetaData newMeta;
  stk::tools::copy_meta_with_io_attributes(get_meta(), newMeta);

  expect_equal_meta_datas(get_meta(), newMeta);
}

TEST_F(MeshClone, copyMetaAndBulk)
{
  initialize_mesh_with_parts_and_fields();

  std::shared_ptr<stk::mesh::BulkData> newBulkPtr = build_mesh(get_comm());
  stk::mesh::BulkData& newBulk = *newBulkPtr;
  stk::mesh::MetaData& newMeta = newBulk.mesh_meta_data();

  stk::tools::copy_mesh(get_bulk(), get_meta().universal_part(), newBulk);

  expect_equal_meta_datas(get_meta(), newMeta);
  expect_equal_entity_counts(get_bulk(), newBulk);
}

TEST_F(MeshClone, copyMeshWithPreExistingOrphanNodes)
{
  initialize_mesh_with_parts_and_fields();

  unsigned numOrphansPerProc = 5;
  add_orphan_nodes(numOrphansPerProc);

  std::shared_ptr<stk::mesh::BulkData> newBulkPtr = build_mesh(get_comm());
  stk::mesh::BulkData& newBulk = *newBulkPtr;
  stk::mesh::MetaData& newMeta = newBulk.mesh_meta_data();
  stk::tools::copy_mesh(get_bulk(), get_meta().universal_part(), newBulk);

  expect_equal_meta_datas(get_meta(), newMeta);
  expect_equal_entity_counts(get_bulk(), newBulk);
}

TEST_F(MeshClone, copyMeshWithOrphanNodesOnSharedSide)
{
  if (get_parallel_size() != 2) return;

  initialize_multiproc_multiblock_mesh();

  std::shared_ptr<stk::mesh::BulkData> newBulkPtr = build_mesh(get_comm());
  stk::mesh::BulkData& newBulk = *newBulkPtr;

  stk::tools::copy_mesh(get_bulk(), *get_meta().get_part("block_1"), newBulk);

  std::vector<size_t> count;
  stk::mesh::count_entities(get_meta().locally_owned_part() | get_meta().globally_shared_part(), newBulk, count);
  if (get_parallel_rank() == 0) {
    EXPECT_EQ(count[stk::topology::NODE_RANK], 8u);
  }
  else {
    EXPECT_EQ(count[stk::topology::NODE_RANK], 0u);
  }
}

TEST_F(MeshClone, copyMeshWithOrphanNodesOnOwnedSide)
{
  if (get_parallel_size() != 2) return;

  initialize_multiproc_multiblock_mesh();

  std::shared_ptr<stk::mesh::BulkData> newBulkPtr = build_mesh(get_comm());
  stk::mesh::BulkData& newBulk = *newBulkPtr;

  stk::tools::copy_mesh(get_bulk(), *get_meta().get_part("block_2"), newBulk);

  std::vector<size_t> count;
  stk::mesh::count_entities(get_meta().locally_owned_part() | get_meta().globally_shared_part(), newBulk, count);
  if (get_parallel_rank() == 0) {
    EXPECT_EQ(count[stk::topology::NODE_RANK], 0u);
  }
  else {
    EXPECT_EQ(count[stk::topology::NODE_RANK], 8u);
  }
}

#if defined(__GNUC__) && (__GNUC__ > 4) && !defined(__INTEL_COMPILER)
TEST(MetaDataSize, sizeChanges_needToUpdateCopyMesh)
{
  stk::mesh::MetaData meta;
  EXPECT_GE(632u, sizeof(meta)) << "Size of MetaData changed.  Does mesh copying capability need to be updated?";
}
#endif

}
