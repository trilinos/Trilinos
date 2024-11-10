#include <gtest/gtest.h>
#include "mpi.h"
#include <stk_mesh/base/Comm.hpp>
#include <string>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/BuildMesh.hpp>

#include <stk_mesh/base/GetEntities.hpp>
#include "stk_mesh/base/Field.hpp"

#include <stk_util/parallel/Parallel.hpp>

#include <stk_unit_test_utils/ioUtils.hpp>

#include <stk_io/StkMeshIoBroker.hpp>

#include <iostream>
#include <unistd.h>                     // for unlink
#include <limits>                       // for std::numeric_limits<double>::epsilon

#include "stk_balance/internal/LastStepFieldWriter.hpp"

namespace
{
using stk::unit_test_util::build_mesh;

void read_and_write_mesh_with_added_field_data(const std::string& inputFilename, const std::string& outputFilename, const std::string& fieldName, const double initialVal, const double timeForLastStep)
{
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(MPI_COMM_WORLD);
  stk::mesh::BulkData& bulkData = *bulkPtr;
  stk::mesh::MetaData& meta = bulkData.mesh_meta_data();

  stk::mesh::Field<double>& nodalTestData = meta.declare_field<double>(stk::topology::NODE_RANK, fieldName);
  stk::mesh::put_field_on_mesh(nodalTestData, meta.universal_part(), &initialVal);
  stk::io::set_field_role(nodalTestData, Ioss::Field::TRANSIENT);

  stk::balance::internal::AllStepFieldWriterAutoDecomp ioHelper(bulkData, inputFilename);

  stk::mesh::EntityVector nodes;
  stk::mesh::get_entities(bulkData, stk::topology::NODE_RANK, nodes);

  std::ostringstream os;
  for(stk::mesh::Entity node : nodes)
    if(bulkData.bucket(node).owned())
      os << "node " << bulkData.identifier(node) << " exists.\n";

  std::cerr << os.str();

  ioHelper.set_output_time(timeForLastStep);
  ioHelper.write_mesh(outputFilename);
}

void verify_field_data_is_same_for_all_nodes(const stk::mesh::BulkData& bulkData, const std::string& fieldName, const double initialVal)
{
  const stk::mesh::MetaData& meta = bulkData.mesh_meta_data();
  stk::mesh::EntityVector nodes;
  stk::mesh::get_selected_entities(meta.locally_owned_part(), bulkData.buckets(stk::topology::NODE_RANK), nodes);

  stk::mesh::FieldBase* nodalTestData = meta.get_field(stk::topology::NODE_RANK, fieldName);

  for(stk::mesh::Entity node : nodes)
  {
    double *data = static_cast<double*>(stk::mesh::field_data(*nodalTestData, node));
    EXPECT_EQ(initialVal, *data);
  }
}

void read_mesh_and_verify_field_data_added_correctly(const std::string &outputFilename, const std::string& fieldName, const double initialVal, const double timeForStep1)
{
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(MPI_COMM_WORLD);
  stk::mesh::BulkData& bulkData = *bulkPtr;

  stk::balance::internal::LastStepFieldWriter ioHelper(bulkData, outputFilename);
  EXPECT_NEAR(timeForStep1, ioHelper.get_max_time(), std::numeric_limits<double>::epsilon());
  verify_field_data_is_same_for_all_nodes(bulkData, fieldName, initialVal);
}

TEST(Stk_Balance, read_and_write_stk_mesh_non_ioss_with_auto_decomp_and_compare_field_data)
{
  std::string outputFilename = "output.e";
  std::string inputFilename = "generated:2x2x4";
  std::string fieldName = "nodal_test_data";

  double initialVal = 0.5;
  double timeForStep1 = 0.5;

  read_and_write_mesh_with_added_field_data(inputFilename, outputFilename, fieldName, initialVal, timeForStep1);
  read_mesh_and_verify_field_data_added_correctly(outputFilename, fieldName, initialVal, timeForStep1);
}

size_t get_global_num_nodes_serial(const std::string& inputFilename)
{
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(MPI_COMM_SELF);
  stk::mesh::BulkData& bulkData = *bulkPtr;

  stk::balance::internal::LastStepFieldWriter ioHelper(bulkData, inputFilename);

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(bulkData, counts);
  return counts[stk::topology::NODE_RANK];
}

size_t get_global_num_nodes_parallel(const std::string& inputFilename, MPI_Comm comm)
{
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(comm);
  stk::mesh::BulkData& bulkData = *bulkPtr;

  stk::balance::internal::AllStepFieldWriterAutoDecomp ioHelper(bulkData, inputFilename);

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(bulkData, counts);
  return counts[stk::topology::NODE_RANK];
}

TEST(Stk_Balance, checkParallelAndSerialNumNodesConsistency)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  std::string inputFilename = stk::unit_test_util::get_option("-i", "generated:4x4x4");

  size_t goldGlobalNumNodes = get_global_num_nodes_serial(inputFilename);
  size_t numNodesWIthComm =   get_global_num_nodes_parallel(inputFilename, comm);
  EXPECT_EQ(goldGlobalNumNodes, numNodesWIthComm);
}

}
