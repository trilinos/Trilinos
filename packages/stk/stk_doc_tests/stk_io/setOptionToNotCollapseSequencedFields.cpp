#include <gtest/gtest.h>
#include <unistd.h>                     // for unlink
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/DatabasePurpose.hpp>   // for DatabasePurpose::READ_MESH, etc
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_unit_test_utils/MeshFixture.hpp>

namespace
{

void create_mesh_with__field_1__field_2__field_3(const std::string & filename, MPI_Comm communicator)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, communicator);

    stk::io::StkMeshIoBroker stkIo(communicator);
    stkIo.set_bulk_data(bulk);
    size_t index = stkIo.add_mesh_database("generated:1x1x4", stk::io::READ_MESH);
    stkIo.set_active_mesh(index);
    stkIo.create_input_mesh();

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

    double time = 0.0;
    stkIo.begin_output_step(results_output_index, time);
    stkIo.write_defined_output_fields(results_output_index);
    stkIo.end_output_step(results_output_index);
}

class MultipleNumberedFieldsWithSameBaseName : public stk::unit_test_util::MeshFixture
{
protected:
    MultipleNumberedFieldsWithSameBaseName() : stkIo(get_comm())
    {
        allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    }

    void read_mesh(const std::string &inputFilename)
    {
        stkIo.set_bulk_data(get_bulk());
        stkIo.add_mesh_database(inputFilename, stk::io::READ_MESH);
        stkIo.create_input_mesh();
        stkIo.add_all_mesh_fields_as_input_fields();
        stkIo.populate_bulk_data();
    }

    void TearDown()
    {
        unlink(filename.c_str());
    }
protected:
    const std::string filename = "multipleNumberedFieldsWithSameBaseName.exo";
    stk::io::StkMeshIoBroker stkIo;
};

//-BEGIN
TEST_F(MultipleNumberedFieldsWithSameBaseName, whenReading_collapseToSingleStkField)
{
    create_mesh_with__field_1__field_2__field_3(filename, get_comm());
    read_mesh(filename);
    EXPECT_EQ(1u, get_meta().get_fields(stk::topology::ELEM_RANK).size());
}

TEST_F(MultipleNumberedFieldsWithSameBaseName, whenReadingWitlhoutCollapseOption_threeStkFieldsAreRead)
{
    create_mesh_with__field_1__field_2__field_3(filename, get_comm());
    stkIo.set_option_to_not_collapse_sequenced_fields();
    read_mesh(filename);
    EXPECT_EQ(3u, get_meta().get_fields(stk::topology::ELEM_RANK).size());
}
//-END

}
