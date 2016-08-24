#include <ngp/Ngp.hpp>
#include <ngp/NgpMultistateField.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_util/stk_config.h>


void set_field_on_device(stk::mesh::BulkData &bulk, stk::mesh::EntityRank rank, stk::mesh::Part &quadPart, stk::mesh::Field<double> &quadField)
{
    ngp::Field<double> ngpQuadField(bulk, quadField);
    ngp::Mesh ngpMesh(bulk);
    ngp::for_each_entity_run(ngpMesh, rank, quadPart, KOKKOS_LAMBDA(ngp::Mesh::MeshIndex entity)
    {
        ngpQuadField.get(entity, 0) = 13.0;
    });
    ngpQuadField.copy_device_to_host(bulk, quadField);
}

class NgpHowTo : public stk::unit_test_util::MeshFixture {};

TEST_F(NgpHowTo, loopOverSubsetOfMesh)
{
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::Part &quadPart = get_meta().get_topology_root_part(stk::topology::SHELL_QUAD_4);
    auto &quadField = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::ELEM_RANK, "myField");
    double init = 0.0;
    stk::mesh::put_field(quadField, quadPart, &init);
    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
         0,2,SHELL_QUAD_4,5,6,7,8";
    stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());

    set_field_on_device(get_bulk(), stk::topology::ELEM_RANK, quadPart, quadField);

    for(const stk::mesh::Bucket *bucket : get_bulk().get_buckets(stk::topology::ELEM_RANK, quadPart))
        for(stk::mesh::Entity elem : *bucket)
            EXPECT_EQ(13.0, *stk::mesh::field_data(quadField, elem));
}

TEST_F(NgpHowTo, loopOverNodes)
{
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    auto &field = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "myField");
    double init = 0.0;
    stk::mesh::put_field(field, get_meta().universal_part(), &init);
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
    stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());

    set_field_on_device(get_bulk(), stk::topology::NODE_RANK, get_meta().universal_part(), field);

    for(const stk::mesh::Bucket *bucket : get_bulk().get_buckets(stk::topology::NODE_RANK, get_meta().universal_part()))
        for(stk::mesh::Entity node : *bucket)
            EXPECT_EQ(13.0, *stk::mesh::field_data(field, node));
}

TEST_F(NgpHowTo, loopOverFaces)
{
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::Part &facePart = get_meta().declare_part("facePart", stk::topology::FACE_RANK);
    auto &field = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::FACE_RANK, "myField");
    double init = 0.0;
    stk::mesh::put_field(field, facePart, &init);
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
    stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());

    stk::mesh::create_exposed_block_boundary_sides(get_bulk(), get_meta().universal_part(), {&facePart});

    set_field_on_device(get_bulk(), stk::topology::FACE_RANK, facePart, field);

    for(const stk::mesh::Bucket *bucket : get_bulk().get_buckets(stk::topology::FACE_RANK, get_meta().universal_part()))
        for(stk::mesh::Entity node : *bucket)
            EXPECT_EQ(13.0, *stk::mesh::field_data(field, node));
}

unsigned count_num_elems(ngp::Mesh ngpMesh,
                         ngp::Field<double> ngpField,
                         stk::mesh::EntityRank rank,
                         stk::mesh::Part &part)
{
    Kokkos::View<unsigned*> numElems("numElems", 1);
    ngp::for_each_entity_run(ngpMesh, rank, part, KOKKOS_LAMBDA(ngp::Mesh::MeshIndex entity)
    {
        unsigned fieldValue = static_cast<unsigned>(ngpField.get(entity, 0));
        Kokkos::atomic_add(&numElems(0), fieldValue);
    });
    Kokkos::View<unsigned*>::HostMirror numElemsHost = Kokkos::create_mirror_view(numElems);
    Kokkos::deep_copy(numElemsHost, numElems);
    return numElemsHost(0);
}

void set_num_elems_in_field_on_device(stk::mesh::BulkData &bulk,
                         stk::mesh::Part &part,
                         stk::mesh::Field<double> &field)
{
    ngp::Field<double> ngpField(bulk, field);
    ngp::Mesh ngpMesh(bulk);
    unsigned numElems = count_num_elems(ngpMesh, ngpField, field.entity_rank(), part);
    ngp::for_each_entity_run(ngpMesh, field.entity_rank(), part, KOKKOS_LAMBDA(ngp::Mesh::MeshIndex entity)
    {
        ngpField.get(entity, 0) = numElems;
    });
    ngpField.copy_device_to_host(bulk, field);
}

TEST_F(NgpHowTo, exerciseAura)
{
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    auto &field = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::ELEM_RANK, "myField");
    double init = 1.0;
    stk::mesh::put_field(field, get_meta().universal_part(), &init);
    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
         1,2,HEX_8,5,6,7,8,9,10,11,12";
    stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());

    set_num_elems_in_field_on_device(get_bulk(), get_meta().universal_part(), field);

    for(const stk::mesh::Bucket *bucket : get_bulk().get_buckets(stk::topology::ELEM_RANK, get_meta().universal_part()))
        for(stk::mesh::Entity elem : *bucket)
            EXPECT_EQ(2.0, *stk::mesh::field_data(field, elem));
}

stk::mesh::Field<int> &create_field_with_num_states(stk::mesh::MetaData &meta, int numStates)
{
    auto &field = meta.declare_field<stk::mesh::Field<int>>(stk::topology::ELEM_RANK, "myField", numStates);
    int init = -1;
    stk::mesh::put_field(field, meta.universal_part(), &init);
    return field;
}

stk::mesh::Field<int> &create_field(stk::mesh::MetaData &meta)
{
    return create_field_with_num_states(meta, 1);
}

void verify_states_np1_and_n_have_values(stk::mesh::BulkData &bulk,
                                         stk::mesh::Field<int>& field,
                                         ngp::MultistateField<int> &ngpMultistateField,
                                         int np1Value,
                                         int nValue)
{
    ngpMultistateField.copy_device_to_host(bulk, field);
    for(const stk::mesh::Bucket* bucket : bulk.buckets(stk::topology::ELEM_RANK))
    {
        for(stk::mesh::Entity elem : *bucket)
        {
            EXPECT_EQ(np1Value, *static_cast<int*>(stk::mesh::field_data(*field.field_state(stk::mesh::StateNP1), elem)));
            EXPECT_EQ(nValue, *static_cast<int*>(stk::mesh::field_data(*field.field_state(stk::mesh::StateN), elem)));
        }
    }
}

void set_states_in_field_on_device(ngp::Mesh &ngpMesh,
                         stk::mesh::EntityRank rank,
                         stk::mesh::Selector sel,
                         ngp::MultistateField<int> &ngpMultistateField)
{
    ngp::for_each_entity_run(ngpMesh, rank, sel, KOKKOS_LAMBDA(ngp::Mesh::MeshIndex entity)
    {
        ngp::Field<int> stateNp1Field = ngpMultistateField.get_field_of_state(stk::mesh::StateNP1);
        stateNp1Field.get(entity, 0) = 1;
        ngp::Field<int> stateNField = ngpMultistateField.get_field_of_state(stk::mesh::StateN);
        stateNField.get(entity, 0) = 2;
    });
}

TEST_F(NgpHowTo, useMultistateFields)
{
    stk::mesh::Field<int> &stkField = create_field_with_num_states(get_meta(), 2);
    setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);

    ngp::MultistateField<int> ngpMultistateField(get_bulk(), stkField);
    ngp::Mesh ngpMesh(get_bulk());

    set_states_in_field_on_device(ngpMesh, stk::topology::ELEM_RANK, get_meta().universal_part(), ngpMultistateField);

    verify_states_np1_and_n_have_values(get_bulk(), stkField, ngpMultistateField, 1, 2);
    get_bulk().update_field_data_states();
    ngpMultistateField.increment_state();
    verify_states_np1_and_n_have_values(get_bulk(), stkField, ngpMultistateField, 2, 1);
}

void set_states_in_convenient_field_on_device(ngp::Mesh &ngpMesh,
                         stk::mesh::EntityRank rank,
                         stk::mesh::Selector sel,
                         ngp::ConvenientMultistateField<int> &ngpMultistateField)
{
    ngp::for_each_entity_run(ngpMesh, rank, sel, KOKKOS_LAMBDA(ngp::Mesh::MeshIndex entity)
    {
        ngpMultistateField.get(stk::mesh::StateNP1, entity, 0) = 1;
        ngpMultistateField.get(stk::mesh::StateN, entity, 0) = 2;
    });
}

TEST_F(NgpHowTo, useConvenientMultistateFields)
{
    stk::mesh::Field<int> &stkField = create_field_with_num_states(get_meta(), 2);
    setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);

    ngp::ConvenientMultistateField<int> ngpMultistateField(get_bulk(), stkField);
    ngp::Mesh ngpMesh(get_bulk());

    set_states_in_convenient_field_on_device(ngpMesh, stk::topology::ELEM_RANK, get_meta().universal_part(), ngpMultistateField);

    verify_states_np1_and_n_have_values(get_bulk(), stkField, ngpMultistateField, 1, 2);
    get_bulk().update_field_data_states();
    ngpMultistateField.increment_state();
    verify_states_np1_and_n_have_values(get_bulk(), stkField, ngpMultistateField, 2, 1);
}

class NgpReduceHowTo : public stk::unit_test_util::MeshFixture
{
protected:
    NgpReduceHowTo()
    {
        stk::mesh::Field<int> &elemField = create_field(get_meta());
        setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);
        stk::mesh::EntityVector elems;
        stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, elems);
        for(stk::mesh::Entity elem : elems)
        {
            int *fieldData = stk::mesh::field_data(elemField, elem);
            fieldData[0] = get_bulk().identifier(elem);
        }
        ngpMesh = ngp::Mesh(get_bulk());
        ngpElemField = ngp::Field<int>(get_bulk(), elemField);
    }
    int get_num_elems()
    {
        return stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::ELEM_RANK));
    }
    ngp::Mesh ngpMesh;
    ngp::Field<int> ngpElemField;
};

int get_min_field_value(ngp::Mesh &ngpMesh, ngp::Field<int> &ngpElemField, stk::mesh::Selector sel)
{
    return ngp::get_field_min(ngpMesh, ngpElemField, sel);
}
TEST_F(NgpReduceHowTo, getMinFieldValue)
{
    EXPECT_EQ(1, get_min_field_value(ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part()));
}

int get_max_field_value(ngp::Mesh &ngpMesh, ngp::Field<int> &ngpElemField, stk::mesh::Selector sel)
{
    return ngp::get_field_max(ngpMesh, ngpElemField, sel);
}
TEST_F(NgpReduceHowTo, getMaxFieldValue)
{
    EXPECT_EQ(get_num_elems(), get_max_field_value(ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part()));
}

int get_sum_field_value(ngp::Mesh &ngpMesh, ngp::Field<int> &ngpElemField, stk::mesh::Selector sel)
{
    return ngp::get_field_sum(ngpMesh, ngpElemField, sel);
}
TEST_F(NgpReduceHowTo, getSumFieldValue)
{
    int numElems = get_num_elems();
    EXPECT_EQ(numElems*(numElems+1)/2, get_sum_field_value(ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part()));
}
