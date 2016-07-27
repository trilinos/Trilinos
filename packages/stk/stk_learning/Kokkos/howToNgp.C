#include "StaticMesh.h"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/stk_config.h>


void set_field_on_device(stk::mesh::BulkData &bulk, stk::mesh::Part &quadPart, stk::mesh::Field<double> &quadField)
{
    ngp::StkNgpField ngpQuadField(bulk, quadField);
    ngp::StkNgpMesh ngpMesh(bulk);
    ngp::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, quadPart, KOKKOS_LAMBDA(ngp::StkNgpMesh::MeshIndex elem)
    {
        ngpQuadField.get(elem, 0) = 13.0;
    });
    ngpQuadField.copy_device_to_host(bulk, quadField);
}

class NgpHowTo : public stk::unit_test_util::MeshFixture
{
};

TEST_F(NgpHowTo, loopOverSubsetOfMesh)
{
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::Part &quadPart = get_meta().get_topology_root_part(stk::topology::SHELL_QUAD_4);
    auto &quadField = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::ELEM_RANK, "quadField");
    double init = 0.0;
    stk::mesh::put_field(quadField, quadPart, &init);
    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
         0,2,SHELL_QUAD_4,5,6,7,8";
    stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());

    set_field_on_device(get_bulk(), quadPart, quadField);

    for(const stk::mesh::Bucket *bucket : get_bulk().get_buckets(stk::topology::ELEM_RANK, quadPart))
        for(stk::mesh::Entity elem : *bucket)
            EXPECT_EQ(13.0, *stk::mesh::field_data(quadField, elem));
}
