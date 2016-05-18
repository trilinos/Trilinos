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
#include <stk_io/IossBridge.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/environment/ReportHandler.hpp>

namespace stk {
namespace mesh {

stk::mesh::Part *create_new_part(const std::string &name, stk::topology topo, stk::mesh::MetaData &newMeta)
{
    if(topo != stk::topology::INVALID_TOPOLOGY)
        return &newMeta.declare_part_with_topology(name, topo);
    else
        return &newMeta.declare_part(name);
}

stk::mesh::OrdinalVector get_part_supersets(const stk::mesh::Part &part)
{
    const stk::mesh::PartVector &supersetParts = part.supersets();
    stk::mesh::OrdinalVector supersetOrdinals(supersetParts.size());
    for(size_t i=0; i<supersetParts.size(); i++)
        supersetOrdinals[i] = supersetParts[i]->mesh_meta_data_ordinal();
    return supersetOrdinals;
}

void copy_part_supersets(const stk::mesh::Part &oldPart, stk::mesh::Part &newPart, stk::mesh::MetaData &newMeta)
{
    stk::mesh::OrdinalVector oldSupersets = get_part_supersets(oldPart);
    for(stk::mesh::PartOrdinal partOrd : oldSupersets)
        newMeta.declare_part_subset(newMeta.get_part(partOrd), newPart);
}

void copy_io_part_attributes(const stk::mesh::Part &oldPart, stk::mesh::Part &newPart)
{
    if(stk::io::is_part_io_part(oldPart))
    {
        stk::io::put_io_part_attribute(newPart);
    }
}

void clone_part_to_other_meta(const stk::mesh::Part &oldPart, stk::mesh::MetaData &newMeta)
{
    stk::mesh::Part *newPart = create_new_part(oldPart.name(), oldPart.topology(), newMeta);
    newMeta.set_part_id(*newPart, oldPart.id());
    copy_part_supersets(oldPart, *newPart, newMeta);
    copy_io_part_attributes(oldPart, *newPart);
}

void copy_parts(const stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
    const stk::mesh::PartVector &allParts = oldMeta.get_mesh_parts();
    for(size_t i = 0; i < allParts.size(); i++)
        clone_part_to_other_meta(*allParts[i], newMeta);
}



stk::mesh::FieldBase* clone_field(const stk::mesh::FieldBase& field, stk::mesh::MetaData& newMeta)
{
    return newMeta.declare_field_base(field.name(),
                                      field.entity_rank(),
                                      field.data_traits(),
                                      field.field_array_rank(),
                                      field.dimension_tags(),
                                      field.number_of_states());
}

void copy_field_restrictions(const stk::mesh::FieldBase& field, stk::mesh::MetaData& newMeta, stk::mesh::FieldBase* newField)
{
    const stk::mesh::FieldRestrictionVector& oldRestrictions = field.restrictions();
    for(const stk::mesh::FieldRestriction& res : oldRestrictions)
    {
//        stk::mesh::Selector selectNewParts = res.selector().clone_for_different_mesh(newMeta);
        stk::mesh::Selector selectNewParts = res.selector(); // WRONG: needs to use above function
        newMeta.declare_field_restriction(*newField,
                                          selectNewParts,
                                          res.num_scalars_per_entity(),
                                          res.dimension(),
                                          field.get_initial_value());
    }
}

void copy_fields(const stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
    const stk::mesh::FieldVector &fields = oldMeta.get_fields();
    for(size_t i = 0; i < fields.size(); i++)
    {
        stk::mesh::FieldBase* newField = clone_field(*fields[i], newMeta);
        copy_field_restrictions(*fields[i], newMeta, newField);
    }
}

// TODO: need this in Selector class to clone correctly
//stk::mesh::Selector clone_for_different_mesh(const stk::mesh::MetaData &differentMeta) const;
//
//stk::mesh::Selector Selector::clone_for_different_mesh(const stk::mesh::MetaData &differentMeta) const
//{
//    stk::mesh::Selector newSelector(*this);
//    for(SelectorNode &selectorNode : newSelector.m_expr)
//    {
//        if(selectorNode.m_type == SelectorNodeType::PART)
//        {
//            unsigned ord = selectorNode.m_value.part_ptr->mesh_meta_data_ordinal();
//            ThrowRequireMsg(selectorNode.m_value.part_ptr->name() == differentMeta.get_part(ord).name(),
//                            "Attepting to clone selector into mesh with different parts");
//            selectorNode.m_value.part_ptr = &differentMeta.get_part(ord);
//        }
//        else if(selectorNode.m_type == SelectorNodeType::FIELD)
//        {
//            unsigned ord = selectorNode.m_value.field_ptr->mesh_meta_data_ordinal();
//            ThrowRequireMsg(selectorNode.m_value.field_ptr->name() == differentMeta.get_fields()[ord]->name(),
//                            "Attepting to clone selector into mesh with different parts");
//            selectorNode.m_value.field_ptr = differentMeta.get_fields()[ord];
//        }
//    }
//    return newSelector;
//}

void copy_mesh_subset(stk::mesh::BulkData &inputBulk, stk::mesh::Selector selectedElems, stk::mesh::BulkData &outputBulk)
{
    ThrowRequireMsg(&inputBulk != &outputBulk, "Can't copy to same mesh.");
    stk::mesh::MetaData &inputMeta = inputBulk.mesh_meta_data();
    stk::mesh::MetaData &outputMeta = outputBulk.mesh_meta_data();
    outputMeta.initialize(inputMeta.spatial_dimension(), inputMeta.entity_rank_names());
    copy_parts(inputMeta, outputMeta);
    copy_fields(inputMeta, outputMeta);
}

}
}

namespace
{

class CopyMeshEmpty : public stk::unit_test_util::MeshFixture {};

TEST_F(CopyMeshEmpty, throwIfSameMesh)
{
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    ASSERT_THROW(stk::mesh::copy_mesh_subset(get_bulk(), stk::mesh::Selector(), get_bulk()), std::exception);
}

class CopyMesh : public stk::unit_test_util::MeshFixture
{
protected:
    const char *get_mesh_spec() const {return "generated:1x1x8";}
    stk::mesh::BulkData::AutomaticAuraOption get_aura_option() const {return stk::mesh::BulkData::NO_AUTO_AURA;}
};

TEST_F(CopyMesh, copyParts)
{
    stk::mesh::Part &superPart = get_meta().declare_part("superPart");
    stk::mesh::Part &subPart = get_meta().declare_part("subPart");
    get_meta().declare_part_subset(superPart, subPart);
    setup_mesh(get_mesh_spec(), get_aura_option());

    stk::mesh::MetaData newMeta;
    stk::mesh::BulkData newBulk(newMeta, get_bulk().parallel());
    stk::mesh::copy_mesh_subset(get_bulk(), stk::mesh::Selector(), newBulk);

    const stk::mesh::PartVector &oldParts = get_meta().get_parts();
    const stk::mesh::PartVector &newParts = newMeta.get_parts();
    ASSERT_EQ(oldParts.size(), newParts.size());
    for(size_t i=0; i<oldParts.size(); i++)
    {
        EXPECT_EQ(oldParts[i]->name(), newParts[i]->name());
        EXPECT_EQ(oldParts[i]->id(), newParts[i]->id());
        EXPECT_EQ(oldParts[i]->topology(), newParts[i]->topology());
        EXPECT_EQ(oldParts[i]->primary_entity_rank(), newParts[i]->primary_entity_rank());
        EXPECT_EQ(oldParts[i]->force_no_induce(), newParts[i]->force_no_induce());
        EXPECT_EQ(stk::io::is_part_io_part(*oldParts[i]), stk::io::is_part_io_part(*newParts[i]));
        stk::mesh::OrdinalVector oldSupersets = get_part_supersets(*oldParts[i]);
        stk::mesh::OrdinalVector newSupersets = get_part_supersets(*newParts[i]);
        EXPECT_TRUE(oldSupersets == newSupersets) << oldParts[i]->name();
    }
}

TEST_F(CopyMesh, copyFields)
{
    stk::mesh::Field<double> &myField1 = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "myField1", 1);
    stk::mesh::Field<double> &myField2 = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "myField2", 1);
    stk::mesh::Part &myPart1 = get_meta().declare_part("myPart1");
    stk::mesh::Part &myPart2 = get_meta().declare_part("myPart2");
    const double initValue = 13;
    stk::mesh::put_field(myField1, myPart1, 2, &initValue);
    stk::mesh::put_field(myField1, myPart2, 3, &initValue);
    stk::mesh::put_field(myField2, myPart1 & myPart2, 3, &initValue);

    setup_mesh(get_mesh_spec(), get_aura_option());

    stk::mesh::MetaData newMeta;
    stk::mesh::BulkData newBulk(newMeta, get_bulk().parallel());
    stk::mesh::copy_mesh_subset(get_bulk(), stk::mesh::Selector(), newBulk);

    const stk::mesh::FieldVector &oldFields = get_meta().get_fields();
    const stk::mesh::FieldVector &newFields = newMeta.get_fields();
    ASSERT_EQ(oldFields.size(), newFields.size());
    for(size_t i=0; i<oldFields.size(); i++)
    {
        EXPECT_EQ(oldFields[i]->name(), newFields[i]->name());
        EXPECT_EQ(oldFields[i]->entity_rank(), newFields[i]->entity_rank());
        EXPECT_EQ(oldFields[i]->field_array_rank(), newFields[i]->field_array_rank());
        EXPECT_EQ(oldFields[i]->number_of_states(), newFields[i]->number_of_states());

        ASSERT_EQ(oldFields[i]->get_initial_value_num_bytes(), newFields[i]->get_initial_value_num_bytes());
        const char *oldInitValPtr = static_cast<const char *>(oldFields[i]->get_initial_value());
        const char *newInitValPtr = static_cast<const char *>(newFields[i]->get_initial_value());
        for(unsigned j=0; j<oldFields[i]->get_initial_value_num_bytes(); j++)
            EXPECT_EQ(oldInitValPtr[j], newInitValPtr[j]);

        const stk::mesh::FieldRestrictionVector &oldRestrictions = oldFields[i]->restrictions();
        const stk::mesh::FieldRestrictionVector &newRestrictions = newFields[i]->restrictions();
        ASSERT_EQ(oldRestrictions.size(), newRestrictions.size());
        for(unsigned j=0; j<oldRestrictions.size(); j++)
        {
            std::ostringstream oldS;
            oldS << oldRestrictions[j].selector();
            std::ostringstream newS;
            newS << newRestrictions[j].selector();
            EXPECT_EQ(oldS.str(), newS.str());

            EXPECT_EQ(oldRestrictions[j].num_scalars_per_entity(), newRestrictions[j].num_scalars_per_entity());
            EXPECT_EQ(oldRestrictions[j].dimension(), newRestrictions[j].dimension());
        }
    }
}

}
