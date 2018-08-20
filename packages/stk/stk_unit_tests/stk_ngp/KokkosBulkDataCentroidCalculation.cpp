#include "KokkosBulkDataCentroidCalculation.hpp"

void calculate_centroids_on_host(const stk::mesh::BulkData& bulkData, const CoordFieldType& coordinates, CoordFieldType& centroid, stk::mesh::Selector selector)
{
    const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::ELEM_RANK, selector);
    for(size_t i=0;i<buckets.size();++i)
    {
        const stk::mesh::Bucket& bucket = *buckets[i];
        for(size_t j=0;j<bucket.size();++j)
        {
            stk::mesh::Entity element = bucket[j];

            unsigned num_nodes = bulkData.num_nodes(element);
            const stk::mesh::Entity* nodes = bulkData.begin_nodes(element);
            double xave = 0, yave = 0, zave = 0;
            for(unsigned k=0;k<num_nodes;++k)
            {
                double *node_coords = stk::mesh::field_data(coordinates, nodes[k]);
                xave += node_coords[0];
                yave += node_coords[1];
                zave += node_coords[2];
            }
            xave /= num_nodes;
            yave /= num_nodes;
            zave /= num_nodes;
            double *centroid_values = stk::mesh::field_data(centroid, element);
            centroid_values[0] = xave;
            centroid_values[1] = yave;
            centroid_values[2] = zave;
        }
    }
}

TEST_F(NGP_Kokkos, calculate_centroid_field_on_host)
{
    MyApp app;

    app.start_timer();
    calculate_centroids_on_host(*app.bulk, *app.coords, app.centroid, app.bulk->mesh_meta_data().locally_owned_part());
    app.stop_timer();
    app.report_bandwidth();

    stk::mesh::Entity element1 = app.bulk->get_entity(stk::topology::ELEM_RANK, 1);
    double *centroid_values = stk::mesh::field_data(app.centroid, element1);
    EXPECT_NEAR(0.5, centroid_values[0], 0.000001);
    EXPECT_NEAR(0.5, centroid_values[1], 0.000001);
    EXPECT_NEAR(0.5, centroid_values[2], 0.000001);
}

