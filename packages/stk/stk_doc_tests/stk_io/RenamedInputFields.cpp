#include <gtest/gtest.h>
#include <unistd.h>
#include <stk_io/WriteMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include "stk_mesh/base/FieldBase.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include "stk_io/MeshField.hpp"

namespace
{

class InputNodesetDataCalledDispX : public stk::unit_test_util::MeshFixture
{
protected:
  ~InputNodesetDataCalledDispX()
  {
    delete stkIo;
  }

  void read_meta(const std::string &filename, double initVal)
  {
    delete stkIo;
    stkIo = new stk::io::StkMeshIoBroker;
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    stkIo->set_bulk_data(get_bulk());
    stkIo->add_mesh_database(filename, stk::io::READ_MESH);
    stkIo->create_input_mesh();
    create_both_fields(initVal);
  }

  void create_both_fields(double initVal)
  {
    nodalDisp = create_field("dispx", get_meta().universal_part(), initVal);

    stk::mesh::Part *nodesetPart = get_meta().get_part("nodelist_1");
    ASSERT_TRUE(nodesetPart != nullptr);
    nsDisp = create_field("input_dispx", *nodesetPart, initVal);
  }

  stk::mesh::FieldBase *create_field(const std::string &name, stk::mesh::Part &part, double initVal)
  {
    stk::mesh::FieldBase *field = &get_meta().declare_field<double>(stk::topology::NODE_RANK, name);
    stk::mesh::put_field_on_mesh(*field, part, &initVal);
    return field;
  }

  void write_output()
  {
    size_t outputFileIndex = stkIo->create_output_mesh(outputFileName, stk::io::WRITE_RESULTS);
    stkIo->add_field(outputFileIndex, *nodalDisp);
    stkIo->add_field(outputFileIndex, *nsDisp, nameOfNodesetVarInFile);
    stkIo->write_output_mesh(outputFileIndex);
    stkIo->begin_output_step(outputFileIndex, 1.0);
    stkIo->write_defined_output_fields(outputFileIndex);
    stkIo->end_output_step(outputFileIndex);
  }

  void read_and_check_field_data_at_step(int step)
  {
    stkIo->populate_bulk_data();
    stkIo->read_defined_input_fields(step);
    expect_both_fields_data();
  }

  void expect_both_fields_data()
  {
    expect_field_data(nodalDisp, initialValue);
    expect_field_data(nsDisp, 6.38e-4);
  }

  void expect_field_data(stk::mesh::FieldBase *field, double expected)
  {
    stk::mesh::Entity node11 = get_bulk().get_entity(stk::topology::NODE_RANK, 11);
    double *data = static_cast<double *>(stk::mesh::field_data(*field, node11));
    ASSERT_TRUE(data != nullptr);
    const double epsilon = 1e-6;
    EXPECT_NEAR(expected, *data, epsilon);
  }

  const std::string inputFileName = "nodesetData.exo";
  const std::string outputFileName = "output.exo";
  double initialValue = 13.0;
  stk::io::StkMeshIoBroker *stkIo = nullptr;

  const std::string nameOfNodesetVarInFile = "dispx";
  stk::mesh::FieldBase *nodalDisp = nullptr;
  stk::mesh::FieldBase *nsDisp = nullptr;
};

TEST_F(InputNodesetDataCalledDispX, writingDispDataToFile_nodalDataDispXAndNodesetDataDispX)
{
  if(stk::parallel_machine_size(get_comm()) == 1)
  {
    read_meta(inputFileName, initialValue);
    stkIo->add_input_field(stk::io::MeshField(nsDisp, nameOfNodesetVarInFile));
    read_and_check_field_data_at_step(100);

    write_output();

    reset_mesh();

    read_meta(outputFileName, 0.0);
    stkIo->add_input_field(stk::io::MeshField(nodalDisp, nodalDisp->name()));
    stkIo->add_input_field(stk::io::MeshField(nsDisp, nameOfNodesetVarInFile));
    read_and_check_field_data_at_step(1);

    unlink(outputFileName.c_str());
  }
}

}
