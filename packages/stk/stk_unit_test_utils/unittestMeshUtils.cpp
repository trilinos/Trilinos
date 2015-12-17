
#include "unittestMeshUtils.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/fixtures/QuadFixture.hpp>  // for QuadFixture
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/Field.hpp>      // for Field
#include <unistd.h>
#include <stk_mesh/base/GetEntities.hpp>
#include "ioUtils.hpp"

namespace stk
{
namespace unit_test_util
{


void put_mesh_into_part(stk::mesh::BulkData& bulkData, stk::mesh::Part& part)
{
    stk::mesh::EntityVector entitiesToMakeActive;
    std::vector<stk::mesh::PartVector> add_parts;
    std::vector<stk::mesh::PartVector> rm_parts;
    for(stk::topology::rank_t rank=stk::topology::BEGIN_RANK; rank < bulkData.mesh_meta_data().entity_rank_count(); rank++)
    {
        const stk::mesh::BucketVector &buckets = bulkData.get_buckets(rank, bulkData.mesh_meta_data().locally_owned_part());
        for(const stk::mesh::Bucket *bucket : buckets)
        {
            for(stk::mesh::Entity entity : *bucket)
            {
                entitiesToMakeActive.push_back(entity);
                add_parts.push_back(stk::mesh::PartVector(1, &part));
                rm_parts.push_back(stk::mesh::PartVector());
            }
        }
    }
    bulkData.batch_change_entity_parts(entitiesToMakeActive, add_parts, rm_parts);
}


std::string get_name_of_generated_mesh(int xdim, int ydim, int zdim, const std::string &options)
{
    std::ostringstream os;
    os << "generated:" << xdim << "x" << ydim << "x" << zdim << options;
    return os.str();
}


void move_killed_elements_out_of_parts(stk::mesh::BulkData& bulkData,
                                  const stk::mesh::EntityVector& killedElements,
                                  const stk::mesh::PartVector& removeParts)
{
    std::vector<stk::mesh::PartVector> add_parts(killedElements.size());
    std::vector<stk::mesh::PartVector> rm_parts(killedElements.size());

    for (size_t j=0;j<killedElements.size();++j)
    {
        rm_parts[j] = removeParts;
    }

    bulkData.batch_change_entity_parts(killedElements, add_parts, rm_parts);
}

void convert_quad_fixture_to_my_bulk_data_flavor(unsigned numX, unsigned numY, stk::mesh::BulkData* bulkData)
{
    stk::mesh::fixtures::QuadFixture fixture(bulkData->parallel(), numX, numY, false);

    stk::mesh::Field<double, stk::mesh::Cartesian2d> &coordField = fixture.m_meta.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian2d>>(stk::topology::NODE_RANK, "model_coordinates");
    stk::mesh::put_field(coordField, fixture.m_meta.universal_part(), fixture.m_meta.spatial_dimension());
    stk::mesh::Part& block_1 = fixture.m_meta.declare_part_with_topology("block_1", stk::topology::QUADRILATERAL_4_2D);
    stk::io::put_io_part_attribute(block_1);

    fixture.m_meta.commit();
    fixture.generate_mesh();

    std::vector<double> x;
    std::vector<double> y;

    for(unsigned j=0;j<=numY;++j)
    {
        for(unsigned i=0;i<=numX;i++)
        {
            x.push_back(i); // 0 1 2, 0 1 2, 0 1 2, ...
            y.push_back(j); // 0 0 0, 1 1 1
        }
    }

    stk::mesh::EntityVector nodes;
    stk::mesh::get_selected_entities(fixture.m_meta.universal_part(), fixture.m_bulk_data.buckets(stk::topology::NODE_RANK), nodes);
    for(stk::mesh::Entity node : nodes )
    {
        double* coords = stk::mesh::field_data(coordField, node);
        unsigned id = fixture.m_bulk_data.identifier(node);
        coords[0] = x[id-1];
        coords[1] = y[id-1];
    }

    fixture.m_bulk_data.modification_begin();
    stk::mesh::EntityVector elems;
    stk::mesh::get_selected_entities(fixture.m_meta.locally_owned_part(), fixture.m_bulk_data.buckets(stk::topology::ELEM_RANK), elems);
    for(stk::mesh::Entity element : elems)
    {
        fixture.m_bulk_data.change_entity_parts(element, {&block_1});
    }
    fixture.m_bulk_data.modification_end();

    std::ostringstream os;
    const std::string file_temp("testadfasdasdfas.exo");
    stk::unit_test_util::write_mesh_using_stk_io(file_temp, fixture.m_bulk_data, bulkData->parallel());
    stk::unit_test_util::fill_mesh_using_stk_io(file_temp, *bulkData, bulkData->parallel());

    ThrowRequireMsg(fixture.m_bulk_data.parallel_size()<10, "Testing assumption violated.");
    os << file_temp << "." << fixture.m_bulk_data.parallel_size() << "." << fixture.m_bulk_data.parallel_rank();
    unlink(os.str().c_str());
}

} // namespace unit_test_util
} // namespace stk

