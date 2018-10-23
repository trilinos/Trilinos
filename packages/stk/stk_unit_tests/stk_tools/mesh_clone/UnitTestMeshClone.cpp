// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_tools/mesh_clone/MeshClone.hpp>
#include <stk_tools/mesh_clone/MeshCloneUtils.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace
{

void expect_equal_parts(const stk::mesh::Part &oldPart, const stk::mesh::Part &newPart)
{
    EXPECT_EQ(oldPart.name(), newPart.name());
    EXPECT_EQ(oldPart.id(), newPart.id()) << oldPart.name();
    EXPECT_EQ(oldPart.topology(), newPart.topology()) << oldPart.name();
    EXPECT_EQ(oldPart.primary_entity_rank(), newPart.primary_entity_rank()) << oldPart.name();
    EXPECT_EQ(oldPart.force_no_induce(), newPart.force_no_induce());
    EXPECT_EQ(stk::io::is_part_io_part(oldPart), stk::io::is_part_io_part(newPart));
    const stk::mesh::FieldBase *oldDistFactField = stk::io::get_distribution_factor_field(oldPart);
    if(oldDistFactField != nullptr)
    {
        const stk::mesh::FieldBase *newDistFactField = stk::io::get_distribution_factor_field(newPart);
        ASSERT_TRUE(newDistFactField != nullptr);
        EXPECT_EQ(oldDistFactField->mesh_meta_data_ordinal(), newDistFactField->mesh_meta_data_ordinal());
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
    EXPECT_EQ(oldField.field_array_rank(), newField.field_array_rank());
    EXPECT_EQ(oldField.number_of_states(), newField.number_of_states());
    EXPECT_EQ(oldField.state(), newField.state());
    EXPECT_EQ(oldField.type_is<int>(), newField.type_is<int>());
    EXPECT_EQ(oldField.type_is<double>(), newField.type_is<double>());
    if(oldField.attribute<Ioss::Field::RoleType>() != nullptr)
    {
        EXPECT_EQ(*oldField.attribute<Ioss::Field::RoleType>(), *newField.attribute<Ioss::Field::RoleType>());
    }
    for(unsigned i = 0; i < oldField.field_array_rank(); ++i)
    {
        EXPECT_EQ(oldField.dimension_tags()[i]->name(), newField.dimension_tags()[i]->name());
    }
    for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<oldField.entity_rank(); rank++)
    {
        EXPECT_EQ(oldField.max_size(rank), newField.max_size(rank));
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


class CopyingMesh : public stk::unit_test_util::MeshFixture
{
protected:
    const char *get_mesh_spec() const {return "generated:1x1x8|sideset:x";}
    stk::mesh::BulkData::AutomaticAuraOption get_aura_option() const {return stk::mesh::BulkData::NO_AUTO_AURA;}

    void initialize_mesh_with_parts_and_fields()
    {
        setup_parts();
        setup_fields();
        setup_mesh(get_mesh_spec(), get_aura_option());
        get_meta().set_coordinate_field(myCoords);
    }

    void setup_parts()
    {
        stk::mesh::Part &superPart = get_meta().declare_part("superPart");
        stk::mesh::Part &subPart = get_meta().declare_part("subPart");
        get_meta().declare_part_subset(superPart, subPart);
    }

    void setup_fields()
    {
        myCoords = &get_meta().declare_field < stk::mesh::Field<double, stk::mesh::Cartesian3d>
                > (stk::topology::NODE_RANK, "myField1", 2);
        auto& myField2 = get_meta().declare_field < stk::mesh::Field<double, stk::mesh::Cartesian3d>
                > (stk::topology::NODE_RANK, "myField2", 1);
        stk::mesh::Part& myPart1 = get_meta().declare_part("myPart1");
        stk::mesh::Part& myPart2 = get_meta().declare_part("myPart2");
        const double initValue[] = {13, 14, 15};
        stk::mesh::put_field_on_mesh(*myCoords, myPart1, 2, initValue);
        stk::mesh::put_field_on_mesh(*myCoords, myPart2, 3, initValue);
        stk::mesh::put_field_on_mesh(myField2, myPart1 & myPart2, 3, initValue);
    }

    void add_orphan_nodes(const unsigned numOrphansPerProc)
    {
        stk::mesh::Part* myPart1 = get_meta().get_part("myPart1");

        std::vector<stk::mesh::EntityId> newIds;
        get_bulk().generate_new_ids(stk::topology::NODE_RANK, numOrphansPerProc, newIds);

        get_bulk().modification_begin();
        for(unsigned i=0; i<numOrphansPerProc; i++)
        {
            stk::mesh::EntityId id = newIds[i];
            get_bulk().declare_entity(stk::topology::NODE_RANK, id, *myPart1);
        }
        get_bulk().modification_end();
    }
    stk::mesh::Field<double, stk::mesh::Cartesian3d> *myCoords;
};

TEST_F(CopyingMesh, ifSameMesh_throws)
{
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    ASSERT_THROW(stk::tools::copy_mesh(get_bulk(), get_meta().universal_part(), get_bulk()), std::exception);
}

TEST_F(CopyingMesh, copyingMetas_same)
{
    initialize_mesh_with_parts_and_fields();

    stk::mesh::MetaData newMeta;
    stk::mesh::BulkData newBulk(newMeta, get_bulk().parallel());

    stk::tools::copy_mesh(get_bulk(), get_meta().universal_part(), newBulk);

    expect_equal_meta_datas(get_meta(), newMeta);
}

TEST_F(CopyingMesh, copyingMeshWithOrphanNodes_same)
{
    initialize_mesh_with_parts_and_fields();

    unsigned numOrphansPerProc = 5;
    add_orphan_nodes(numOrphansPerProc);

    stk::mesh::MetaData newMeta;
    stk::mesh::BulkData newBulk(newMeta, get_bulk().parallel());

    stk::tools::copy_mesh(get_bulk(), get_meta().universal_part(), newBulk);

    expect_equal_meta_datas(get_meta(), newMeta);
    expect_equal_entity_counts(get_bulk(), newBulk);
}

TEST(MetaDataSize, sizeChanges_needToUpdateCopyMesh)
{
    stk::mesh::MetaData meta;
    EXPECT_EQ(520u, sizeof(meta)) << "Size of MetaData changed.  Does mesh copying capability need to be updated?";
}

TEST(MetaData, cloneDoubleField)
{
    std::string fieldName = "dist_fact";
    stk::mesh::MetaData meta(3);
    meta.declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "coordinates");
    meta.declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, fieldName);
    meta.commit();

    stk::mesh::MetaData metaOut;
    stk::tools::copy_meta_with_io_attributes(meta, metaOut);

    stk::mesh::Field<double> *df_field = metaOut.get_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, fieldName);
    EXPECT_TRUE(df_field != nullptr);
}

}
