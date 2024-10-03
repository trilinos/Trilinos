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
#include <stk_mesh/base/MeshBuilder.hpp>

namespace
{

TEST(StkIoHowTo, WriteRestartWithEdges)
{
  std::string filename = "output.rst";
  unsigned numStates = 1;
  int outputTimeStep = 1;
  double outputTime = 0.0;
  {
    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_spatial_dimension(3);
    std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    stk::mesh::Part* part = &meta.declare_part_with_topology("edgeBlock", stk::topology::LINE_2);
    stk::mesh::Field<double>& edgeField = meta.declare_field<double>(stk::topology::EDGE_RANK, "edgeField", numStates);
    stk::mesh::put_field_on_mesh(edgeField, meta.universal_part(), nullptr);
    stk::io::put_edge_block_io_part_attribute(*part);
    stk::io::fill_mesh("generated:1x1x1", *bulk);
    stk::mesh::create_edges(*bulk, meta.universal_part(), part);

    stk::mesh::EntityVector edges;
    stk::mesh::get_entities(*bulk, stk::topology::EDGE_RANK, meta.universal_part(), edges);

    for(auto edge : edges) {
      double* data = reinterpret_cast<double*>(stk::mesh::field_data(edgeField, edge));

      *data = bulk->identifier(edge);
    }

    stk::io::StkMeshIoBroker ioBroker;
    ioBroker.set_bulk_data(bulk);
    stk::io::write_mesh_with_fields(filename, ioBroker, outputTimeStep, outputTime, stk::io::WRITE_RESTART);
  }

  {
    std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    stk::mesh::Field<double>& edgeField = meta.declare_field<double>(stk::topology::EDGE_RANK, "edgeField", numStates);
    stk::mesh::put_field_on_mesh(edgeField, meta.universal_part(), nullptr);

    stk::io::set_field_role(edgeField, Ioss::Field::TRANSIENT);

    stk::io::StkMeshIoBroker stkIo;
    stk::io::fill_mesh_with_fields(filename, stkIo, *bulk, stk::io::READ_RESTART);

    int numSteps = stkIo.get_num_time_steps();
    EXPECT_EQ(1, numSteps);

    stk::mesh::EntityVector edges;
    stk::mesh::get_entities(*bulk, stk::topology::EDGE_RANK, edges);

    for(auto edge : edges) {
      double* data = reinterpret_cast<double*>(stk::mesh::field_data(edgeField, edge));
      double expectedValue = bulk->identifier(edge);

      EXPECT_EQ(expectedValue, *data);
    }
  }

  unlink(filename.c_str());
}

}
