// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_unit_test_utils/meshCreationHelpers.hpp>
#include <stddef.h>                    // for size_t
#include <stk_io/DatabasePurpose.hpp>  // for DatabasePurpose::READ_MESH, etc
#include <stk_io/StkMeshIoBroker.hpp>  // for StkMeshIoBroker
#include <stk_mesh/base/Field.hpp>     // for Field
#include <stk_mesh/base/MetaData.hpp>  // for MetaData, put_field
#include <stk_topology/topology.hpp>   // for topology, etc
#include <string>                      // for string
#include "mpi.h"                       // for MPI_Comm, ompi_communicator_t
#include "stk_mesh/base/BulkData.hpp"  // for BulkData
#include "stk_unit_test_utils/BuildMesh.hpp"
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk
{
namespace unit_test_util
{

size_t write_mesh_data__field_1__field_2__field_3(const std::string & filename, MPI_Comm communicator, stk::mesh::BulkData & bulk, stk::io::StkMeshIoBroker & stkIo)
{
    stkIo.set_bulk_data(bulk);
    size_t index = stkIo.add_mesh_database("generated:1x1x4", stk::io::READ_MESH);
    stkIo.set_active_mesh(index);
    stkIo.create_input_mesh();

    stk::mesh::MetaData &meta = bulk.mesh_meta_data();
    stk::mesh::Field<double> & field1 = meta.declare_field<double>(stk::topology::ELEM_RANK, "field_1", 1);
    stk::mesh::Field<double> & field2 = meta.declare_field<double>(stk::topology::ELEM_RANK, "field_2", 1);
    stk::mesh::Field<double> & field3 = meta.declare_field<double>(stk::topology::ELEM_RANK, "field_3", 1);

    double fieldValues[] = {1.0, 2.0, 3.0};
    stk::mesh::put_field_on_mesh(field1, meta.universal_part(), fieldValues);
    stk::mesh::put_field_on_mesh(field2, meta.universal_part(), fieldValues + 1);
    stk::mesh::put_field_on_mesh(field3, meta.universal_part(), fieldValues + 2);
    stkIo.populate_bulk_data();

    size_t results_output_index = stkIo.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    stkIo.add_field(results_output_index, field1);
    stkIo.add_field(results_output_index, field2);
    stkIo.add_field(results_output_index, field3);

    stkIo.write_output_mesh(results_output_index);
    return results_output_index;
}

void create_mesh_without_time_steps(const std::string & filename, MPI_Comm communicator)
{
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(communicator, stk::mesh::BulkData::AUTO_AURA);
    stk::io::StkMeshIoBroker stkIoWriter(communicator);
    write_mesh_data__field_1__field_2__field_3(filename, communicator, *bulk, stkIoWriter);
}

void create_mesh_with__field_1__field_2__field_3(const std::string & filename, MPI_Comm communicator)
{
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(communicator, stk::mesh::BulkData::AUTO_AURA);
    stk::io::StkMeshIoBroker stkIo(communicator);
    size_t results_output_index = write_mesh_data__field_1__field_2__field_3(filename, communicator, *bulk, stkIo);
    double time = 0.0;
    stkIo.begin_output_step(results_output_index, time);
    stkIo.write_defined_output_fields(results_output_index);
    stkIo.end_output_step(results_output_index);
}

namespace simple_fields {

size_t write_mesh_data__field_1__field_2__field_3(const std::string & filename, MPI_Comm communicator, stk::mesh::BulkData & bulk, stk::io::StkMeshIoBroker & stkIo)
{
    stkIo.set_bulk_data(bulk);
    size_t index = stkIo.add_mesh_database("generated:1x1x4", stk::io::READ_MESH);
    stkIo.set_active_mesh(index);
    stkIo.create_input_mesh();

    stk::mesh::MetaData &meta = bulk.mesh_meta_data();
    stk::mesh::Field<double> & field1 = meta.declare_field<double>(stk::topology::ELEM_RANK, "field_1", 1);
    stk::mesh::Field<double> & field2 = meta.declare_field<double>(stk::topology::ELEM_RANK, "field_2", 1);
    stk::mesh::Field<double> & field3 = meta.declare_field<double>(stk::topology::ELEM_RANK, "field_3", 1);

    double fieldValues[] = {1.0, 2.0, 3.0};
    stk::mesh::put_field_on_mesh(field1, meta.universal_part(), fieldValues);
    stk::mesh::put_field_on_mesh(field2, meta.universal_part(), fieldValues + 1);
    stk::mesh::put_field_on_mesh(field3, meta.universal_part(), fieldValues + 2);
    stkIo.populate_bulk_data();

    size_t results_output_index = stkIo.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    stkIo.add_field(results_output_index, field1);
    stkIo.add_field(results_output_index, field2);
    stkIo.add_field(results_output_index, field3);

    stkIo.write_output_mesh(results_output_index);
    return results_output_index;
}

void create_mesh_without_time_steps(const std::string & filename, MPI_Comm communicator)
{
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(communicator, stk::mesh::BulkData::AUTO_AURA);
    stk::io::StkMeshIoBroker stkIoWriter(communicator);
    write_mesh_data__field_1__field_2__field_3(filename, communicator, *bulk, stkIoWriter);
}

void create_mesh_with__field_1__field_2__field_3(const std::string & filename, MPI_Comm communicator)
{
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(communicator, stk::mesh::BulkData::AUTO_AURA);
    stk::io::StkMeshIoBroker stkIo(communicator);
    size_t results_output_index = write_mesh_data__field_1__field_2__field_3(filename, communicator, *bulk, stkIo);
    double time = 0.0;
    stkIo.begin_output_step(results_output_index, time);
    stkIo.write_defined_output_fields(results_output_index);
    stkIo.end_output_step(results_output_index);
}

} // namespace simple_fields

}
}
