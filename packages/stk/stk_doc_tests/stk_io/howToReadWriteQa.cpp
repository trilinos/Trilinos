#include <gtest/gtest.h>
#include <unistd.h>                     // for unlink
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/DatabasePurpose.hpp>   // for DatabasePurpose::READ_MESH, etc
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/meshCreationHelpers.hpp>

namespace
{

class StkIoHowToQaRecords : public stk::unit_test_util::MeshFixture
{
protected:
  void write_qa_information()
  {
    stk::io::StkMeshIoBroker stkIo(get_comm());
    stkIo.set_bulk_data(get_bulk());
    size_t fileId = stkIo.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    stkIo.set_name_and_version_for_qa_record(fileId, codeName, codeVersion);
    stkIo.add_qa_records(fileId, qaRecords);
    stkIo.add_info_records(fileId, infos);
    stkIo.write_output_mesh(fileId);
  }

  void read_qa_information_and_verify()
  {
    stk::io::StkMeshIoBroker stkIo(get_comm());
    size_t fileId = stkIo.add_mesh_database(filename, stk::io::READ_MESH);
    stkIo.set_active_mesh(fileId);
    std::vector<stk::io::QaRecord> qas = stkIo.get_qa_records();
    std::vector<std::string> infosGotten = stkIo.get_info_records();

    verify_qa_information(qas, infosGotten);
  }
  void verify_qa_information(const std::vector<stk::io::QaRecord> &qas, const std::vector<std::string> &infosGotten)
  {
    ASSERT_EQ(2u, qas.size());
    EXPECT_EQ(qaRecords[0].name,    qas[0].name);
    EXPECT_EQ(qaRecords[0].version, qas[0].version);
    EXPECT_EQ(qaRecords[0].date,    qas[0].date);
    EXPECT_EQ(qaRecords[0].time,    qas[0].time);
    EXPECT_EQ(codeName, qas[1].name);
    EXPECT_EQ(codeVersion, qas[1].version);
    EXPECT_NE(std::find(infosGotten.begin(), infosGotten.end(), infos[0]), infosGotten.end());
    EXPECT_NE(std::find(infosGotten.begin(), infosGotten.end(), infos[1]), infosGotten.end());
  }

  const std::string filename = "file.exo";
  const std::string codeName = "Name" ;
  const std::string codeVersion = "Version";
  const std::vector<stk::io::QaRecord> qaRecords = {{"codeName", "codeVersion", "Today", "RightNow"}};
  const std::vector<std::string> infos = {"info0", "info1"};
};

TEST_F(StkIoHowToQaRecords, write_read)
{
  setup_mesh("generated:1x1x4", stk::mesh::BulkData::NO_AUTO_AURA);
  write_qa_information();
  read_qa_information_and_verify();
  unlink(filename.c_str());
}

}
