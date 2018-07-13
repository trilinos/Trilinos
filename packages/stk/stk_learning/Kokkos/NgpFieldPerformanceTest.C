#include <gtest/gtest.h>
#include <ngp/Ngp.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/stk_config.h>

namespace ngp_field_perf_test
{

template <typename CoordFieldType>
void calculate_centroid(const ngp::Mesh &ngpMesh, const CoordFieldType &ngpCoords, const stk::mesh::Selector &sel, ngp::Field<double> &ngpCentroid)
{
    ngp::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, sel, KOKKOS_LAMBDA(ngp::Mesh::MeshIndex elem)
    {
        for(size_t count = 0; count < 1000; count++)
        {
            ngp::Mesh::ConnectedNodes nodes = ngpMesh.get_nodes(elem);

            for(unsigned dim = 0; dim < 3; dim++)
                ngpCentroid.get(elem, dim) = 0;

            for(size_t i = 0; i < nodes.size(); i++)
            {
                stk::mesh::FastMeshIndex nodeIndex = ngpMesh.fast_mesh_index(nodes[i]);
                for(unsigned dim = 0; dim < 3; dim++)
                    ngpCentroid.get(elem, dim) += ngpCoords.get(nodeIndex, dim);
            }

            for(unsigned dim = 0; dim < 3; dim++)
                ngpCentroid.get(elem, dim) /= nodes.size();
        }
    });
}

template <typename CoordFieldType>
void calculate_centroid_using_coord_field(const stk::mesh::BulkData &bulk, stk::mesh::FieldBase &centroid)
{
    const stk::mesh::FieldBase& coords = *bulk.mesh_meta_data().coordinate_field();
    CoordFieldType ngpCoords(bulk, coords);
    ngp::Field<double> ngpCentroid(bulk, centroid);
    ngp::Mesh ngpMesh(bulk);

    calculate_centroid(ngpMesh, ngpCoords, bulk.mesh_meta_data().locally_owned_part(), ngpCentroid);

    ngpCentroid.copy_device_to_host(bulk, centroid);
}

std::vector<double> get_centroid_average(stk::mesh::BulkData &bulk, stk::mesh::Field<double, stk::mesh::Cartesian3d> &centroid)
{
    std::vector<double> average = {0, 0, 0};
    size_t numElems = 0;
    for(const stk::mesh::Bucket *bucket : bulk.buckets(stk::topology::ELEM_RANK))
    {
        for(stk::mesh::Entity elem : *bucket)
        {
            double *elemCentroid = stk::mesh::field_data(centroid, elem);
            for(size_t dim = 0; dim < 3; dim++)
                average[dim] += elemCentroid[dim];
            numElems++;
        }
    }

    for(size_t dim = 0; dim < 3; dim++)
        average[dim] /= numElems;

    return average;
}

class NgpFieldPerf : public stk::unit_test_util::MeshFixture
{
protected:
    void declare_centroid_field()
    {
        centroid = &get_meta().declare_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::ELEM_RANK, "centroid");
        stk::mesh::put_field(*centroid, get_meta().universal_part(), 3);
    }
    void generate_mesh()
    {
        std::string meshSpec = stk::unit_test_util::get_mesh_spec("-dim");
        setup_mesh(meshSpec, stk::mesh::BulkData::NO_AUTO_AURA);
    }
    void verify_averaged_centroids_are_center_of_mesh()
    {
        std::vector<double> average = get_centroid_average(get_bulk(), *centroid);
        double meshCenter = stk::unit_test_util::get_command_line_option<double>("-dim", 20) / 2.0;
        for(size_t dim = 0; dim < 3; dim++)
            EXPECT_EQ(meshCenter, average[dim]);
    }
    template <typename CoordFieldType>
    double time_field_data_access()
    {
        double startTime = stk::wall_time();
        calculate_centroid_using_coord_field<CoordFieldType>(get_bulk(), *centroid);
        double duration = stk::wall_time() - startTime;
        verify_averaged_centroids_are_center_of_mesh();
        return duration;
    }
    void warm_up_gpu()
    {
        time_field_data_access<ngp::Field<double>>();
    }

    stk::mesh::Field<double, stk::mesh::Cartesian3d> *centroid;
};


TEST_F(NgpFieldPerf, constFieldDataAccessIsFasterThanFieldDataAccess)
{
    declare_centroid_field();
    generate_mesh();
    warm_up_gpu();
    double nonConstTime = time_field_data_access<ngp::Field<double>>();
    double constTime = time_field_data_access<ngp::ConstField<double>>();
    std::cerr << "non-const time: " << nonConstTime << ", const time: " << constTime << std::endl;
#ifdef KOKKOS_ENABLE_CUDA
    EXPECT_LT(constTime, nonConstTime);
#endif
}

}
