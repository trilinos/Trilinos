#include <unistd.h>
#include <gtest/gtest.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/GetEntities.hpp>

extern int gl_argc;
extern char** gl_argv;

namespace
{

TEST(StkIoHowTo, WriteMeshWithEdges)
{
    std::string filename = "output.exo";
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
        stk::mesh::Part* part = &meta.declare_part_with_topology("edgeBlock", stk::topology::LINE_2);
        stk::io::put_io_part_attribute(*part);
        stk::io::fill_mesh("generated:1x1x1", bulk);
        stk::mesh::create_edges(bulk, meta.universal_part(), part);

        unsigned numEdges = stk::mesh::count_selected_entities(meta.universal_part(), bulk.buckets(stk::topology::EDGE_RANK));
        EXPECT_EQ(12u, numEdges);

        stk::io::StkMeshIoBroker ioBroker;
        ioBroker.set_bulk_data(bulk);
        ioBroker.enable_edge_io();
        stk::io::write_mesh(filename, ioBroker, stk::io::WRITE_RESULTS);
    }

    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
        stk::io::fill_mesh(filename, bulk);

        std::vector<size_t> entityCounts;
        stk::mesh::comm_mesh_counts(bulk, entityCounts);
        EXPECT_EQ(12u, entityCounts[stk::topology::EDGE_RANK]);
        EXPECT_EQ(1u, entityCounts[stk::topology::ELEM_RANK]);

        stk::mesh::EntityVector edges;
        stk::mesh::get_entities(bulk, stk::topology::EDGE_RANK, edges);
        stk::mesh::EntityVector elems;
        stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, elems);

        for(stk::mesh::Entity edge : edges) {
            const stk::mesh::Entity* edgeElement = bulk.begin_elements(edge);
            EXPECT_EQ(1u, bulk.num_elements(edge));
            EXPECT_EQ(edgeElement[0], elems[0]);
        }

        EXPECT_EQ(12u, bulk.num_edges(elems[0]));
    }

    unlink(filename.c_str());
}

}
