#ifndef meshcreationhelpers_hpp
#define meshcreationhelpers_hpp

#include <stk_unit_test_utils/meshCreationHelpers.hpp>
#include <stk_io/DatabasePurpose.hpp>   // for DatabasePurpose::READ_MESH, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_topology/topology.hpp>    // for topology, etc

namespace stk
{
namespace unit_test_util
{

size_t create_mesh_without_time_steps(const std::string & filename, MPI_Comm communicator, stk::mesh::BulkData & bulk, stk::io::StkMeshIoBroker & stkIo)
{
    stkIo.set_bulk_data(bulk);
    size_t index = stkIo.add_mesh_database("generated:1x1x4", stk::io::READ_MESH);
    stkIo.set_active_mesh(index);
    stkIo.create_input_mesh();

    stk::mesh::MetaData &meta = bulk.mesh_meta_data();
    stk::mesh::Field<double> & field1 = meta.declare_field<stk::mesh::Field<double>>(stk::topology::ELEM_RANK, "field_1", 1);
    stk::mesh::Field<double> & field2 = meta.declare_field<stk::mesh::Field<double>>(stk::topology::ELEM_RANK, "field_2", 1);
    stk::mesh::Field<double> & field3 = meta.declare_field<stk::mesh::Field<double>>(stk::topology::ELEM_RANK, "field_3", 1);
    stk::mesh::put_field(field1, meta.universal_part(), 1.0);
    stk::mesh::put_field(field2, meta.universal_part(), 2.0);
    stk::mesh::put_field(field3, meta.universal_part(), 3.0);
    stkIo.populate_bulk_data();

    size_t results_output_index = stkIo.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    stkIo.add_field(results_output_index, field1);
    stkIo.add_field(results_output_index, field2);
    stkIo.add_field(results_output_index, field3);

    stkIo.write_output_mesh(results_output_index);
    return results_output_index;
}

void create_mesh_with__field_1__field_2__field_3(const std::string & filename, MPI_Comm communicator)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, communicator);
    stk::io::StkMeshIoBroker stkIo(communicator);
    size_t results_output_index = create_mesh_without_time_steps(filename, communicator, bulk, stkIo);
    double time = 0.0;
    stkIo.begin_output_step(results_output_index, time);
    stkIo.write_defined_output_fields(results_output_index);
    stkIo.end_output_step(results_output_index);
}

}
}

#endif
