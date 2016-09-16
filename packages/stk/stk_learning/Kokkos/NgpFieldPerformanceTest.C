#include <gtest/gtest.h>
#include <ngp/Ngp.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/stk_config.h>

namespace
{

void calculate_centroid(const ngp::Mesh &ngpMesh, const ngp::Field<double> &ngpCoords, const stk::mesh::Selector &sel, ngp::Field<double> &ngpCentroid)
{
    ngp::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, sel, KOKKOS_LAMBDA(ngp::Mesh::MeshIndex elem)
    {
       ngp::Mesh::ConnectedNodes nodes = ngpMesh.get_nodes(elem);

       for(unsigned dim = 0; dim < 3; dim++)
       {
           for(size_t i = 0; i < nodes.size(); i++)
           {
               ngpCentroid.get(elem, dim) += ngpCoords.get(ngpMesh, nodes[i], dim);
           }
           ngpCentroid.get(elem, dim) /= nodes.size();
       }
    });
}

class NgpFieldPerf : public stk::unit_test_util::MeshFixture {};

TEST_F(NgpFieldPerf, fieldDataAccess)
{
    stk::mesh::Field<double, stk::mesh::Cartesian3d> &centroid = get_meta().declare_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::ELEM_RANK, "centroid");
    stk::mesh::put_field(centroid, get_meta().universal_part(), 3);

    std::string meshSpec = stk::unit_test_util::get_mesh_spec("-dim");
    setup_mesh(meshSpec, stk::mesh::BulkData::NO_AUTO_AURA);

    const stk::mesh::FieldBase& coords = *get_meta().coordinate_field();
    ngp::Field<double> ngpCoords(get_bulk(), coords);
    ngp::Field<double> ngpCentroid(get_bulk(), centroid);
    ngp::Mesh ngpMesh(get_bulk());

    calculate_centroid(ngpMesh, ngpCoords, get_meta().locally_owned_part(), ngpCentroid);

    //verify_averaged_centroids_are_center_of_mesh(ngpCentroid)
    double average[3] = {0, 0, 0};
    size_t numElems = 0;
    for(const stk::mesh::Bucket *bucket : get_bulk().buckets(stk::topology::ELEM_RANK))
    {
        for(stk::mesh::Entity elem : *bucket)
        {
            double *elemCentroid = stk::mesh::field_data(centroid, elem);
            for(size_t dim = 0; dim < 3; dim++)
                average[dim] += elemCentroid[dim];
            numElems++;
        }
    }

    double meshCenter = unitTestUtils::get_command_line_option<double>("-dim", "20") / 2.0;
    for(size_t dim = 0; dim < 3; dim++)
    {
        average[dim] /= numElems;
        EXPECT_EQ(meshCenter, average[dim]);
    }
}

TEST_F(NgpFieldPerf, constFieldDataAccess)
{
}

}
