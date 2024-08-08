// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocks.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocksImpl.hpp>
#include <stk_tools/mesh_tools/DetectHingesImpl.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_util/environment/WallTime.hpp>
#include "stk_unit_test_utils/getOption.h"
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/EntityLess.hpp"
#include "stk_mesh/base/Comm.hpp"
#include "stk_balance/m2n/TransientFieldTransferById.hpp"
#include "stk_balance/m2n/M2NOutputSerializerBulkData.hpp"
#include "stk_balance/m2n/M2NSubdomainWriter.hpp"
#include "Ioss_Field.h"
#include <stk_tools/mesh_clone/MeshClone.hpp>

#include <unistd.h>
#include <string>
#include <algorithm>

namespace {

using stk::unit_test_util::build_mesh;

void expect_and_unlink_file(const std::string& baseName, int numProc, int procId)
{
  std::string fileName = stk::io::construct_filename_for_serial_or_parallel(baseName, numProc, procId);
  std::ifstream file(fileName);
  EXPECT_TRUE(!file.fail());
  unlink(fileName.c_str());
}

stk::mesh::EntityIdVector get_id_vector(const stk::mesh::BulkData& bulk, const stk::mesh::EntityVector& entityVec)
{
  stk::mesh::EntityIdVector entityIdVector;

  for(stk::mesh::Entity entity : entityVec) {
    entityIdVector.push_back(bulk.identifier(entity));
  }
  return entityIdVector;
}

bool are_entities_equivalent(const stk::mesh::BulkData& bulk1, const stk::mesh::BulkData& bulk2, stk::mesh::EntityRank rank)
{
  const stk::mesh::MetaData& meta1 = bulk1.mesh_meta_data();
  const stk::mesh::MetaData& meta2 = bulk2.mesh_meta_data();

  stk::mesh::EntityVector entityVec1, entityVec2;
  stk::mesh::get_selected_entities(meta1.locally_owned_part(), bulk1.buckets(rank), entityVec1);
  stk::mesh::get_selected_entities(meta2.locally_owned_part(), bulk2.buckets(rank), entityVec2);
  stk::mesh::EntityIdVector entityIdVec1 = get_id_vector(bulk1, entityVec1);
  stk::mesh::EntityIdVector entityIdVec2 = get_id_vector(bulk2, entityVec2);

  return entityIdVec1 == entityIdVec2;
}

bool are_bulk_data_equivalent(const stk::mesh::BulkData& bulk1, const stk::mesh::BulkData& bulk2)
{
  const stk::mesh::MetaData& meta1 = bulk1.mesh_meta_data();
  const stk::mesh::MetaData& meta2 = bulk2.mesh_meta_data();

  if(meta1.entity_rank_names() != meta2.entity_rank_names()) {
    return false;
  }

  std::vector<size_t> entityCountVec1, entityCountVec2;
  stk::mesh::count_entities(meta1.locally_owned_part(), bulk1, entityCountVec1);
  stk::mesh::count_entities(meta2.locally_owned_part(), bulk2, entityCountVec2);

  if(entityCountVec1 != entityCountVec2) {
    return false;
  }

  if(!are_entities_equivalent(bulk1, bulk2, stk::topology::NODE_RANK)) {
    return false;
  }

  if(!are_entities_equivalent(bulk1, bulk2, stk::topology::ELEMENT_RANK)) {
    return false;
  }

  return true;
}

stk::io::EntitySharingInfo get_four_hex_mesh_node_sharing_info(const stk::mesh::BulkData& inputBulk, int subdomain)
{
  unsigned numProcs = inputBulk.parallel_size();
  STK_ThrowRequireMsg( (numProcs == 1 || numProcs == 2), "Invalid number of procs: " << numProcs);
  stk::io::EntitySharingInfo nodeSharingInfo;

  if(numProcs == 1) { return nodeSharingInfo; }

  if(inputBulk.parallel_rank() == 0) {
    if(subdomain == 0) {
      nodeSharingInfo = { {5,1}, {6,1}, {7,1}, {8,1} };
    }
    else if(subdomain == 1) {
      nodeSharingInfo = { {5,0}, {6,0}, {7,0}, {8,0}, {9,2}, {10,2}, {11,2}, {12,2} };
    }
  } else {
    if(subdomain == 2) {
      nodeSharingInfo = { {9,1}, {10,1}, {11,1}, {12,1}, {13,3}, {14,3}, {15,3}, {16,3} };
    }
    else if(subdomain == 3) {
      nodeSharingInfo = { {13,2}, {14,2}, {15,2}, {16,2} };
    }
  }
  return nodeSharingInfo;
}

stk::io::EntitySharingInfo get_two_hex_mesh_node_sharing_info(const stk::mesh::BulkData& inputBulk)
{
  unsigned numProcs = inputBulk.parallel_size();
  STK_ThrowRequireMsg(numProcs == 1 || numProcs == 2, "Invalid number of procs: " << numProcs);
  stk::io::EntitySharingInfo nodeSharingInfo;

  if(numProcs == 2) {
    if(inputBulk.parallel_rank() == 0) {
      nodeSharingInfo = { {5,1}, {6,1}, {7,1}, {8,1} };
    } else {
      nodeSharingInfo = { {5,0}, {6,0}, {7,0}, {8,0} };
    }
  }
  return nodeSharingInfo;
}

void test_static_mesh_output(const stk::mesh::BulkData& inputBulk, const std::string& fileName)
{
  std::shared_ptr<stk::mesh::BulkData> outputBulk = build_mesh(3, MPI_COMM_WORLD);

  stk::io::fill_mesh(fileName, *outputBulk);

  EXPECT_TRUE(are_bulk_data_equivalent(inputBulk, *outputBulk));
}

void clean_up(int targetNumDomains, const std::string& fileName) {
  stk::parallel_machine_barrier(MPI_COMM_WORLD);
  if(stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    for(int i = 0; i < targetNumDomains; i++) {
      expect_and_unlink_file(fileName, targetNumDomains, i);
    }
  }
}

TEST(TransientFieldTransfer, writeStaticMesh)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) { return; }

  stk::balance::m2n::OutputSerializerBulkData bulk(3, MPI_COMM_WORLD);

  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, "generated:1x1x2", bulk);

  stk::io::EntitySharingInfo nodeSharingInfo = get_two_hex_mesh_node_sharing_info(bulk);

  std::string fileName = "test_m2n_output.g";
  int targetNumDomains = bulk.parallel_size();
  int globalNumNodes = 12;
  int globalNumElems = 2;
  int subdomain = bulk.parallel_rank();

  stk::balance::m2n::TransientFieldTransferById m2nIo(ioBroker, targetNumDomains);
  m2nIo.setup_subdomain(bulk, fileName, subdomain, nodeSharingInfo, globalNumNodes, globalNumElems);
  m2nIo.get_subdomain_writer(subdomain).write_mesh();

  test_static_mesh_output(bulk, fileName);
  clean_up(targetNumDomains, fileName);
}

double get_field_data(const stk::mesh::BulkData& bulk, stk::mesh::Entity entity, double time)
{
  stk::mesh::EntityId id = bulk.identifier(entity);
  return 10*time+id;
} 

void set_field_data(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase& field, double time)
{
  stk::mesh::EntityRank rank = field.entity_rank();
  for(stk::mesh::Bucket* bucket : bulk.buckets(rank)) {
    for(stk::mesh::Entity entity : *bucket) {
      double* data = (double*)stk::mesh::field_data(field, entity);
      *data = get_field_data(bulk, entity, time);
    }
  }
}

void test_field_data(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase& field, double time)
{
  stk::mesh::EntityRank rank = field.entity_rank();
  for(stk::mesh::Bucket* bucket : bulk.buckets(rank)) {
    if(!bucket->owned()) { continue; }

    for(stk::mesh::Entity entity : *bucket) {
      double* data = (double*)stk::mesh::field_data(field, entity);
      double expectedValue = get_field_data(bulk, entity, time);
      EXPECT_DOUBLE_EQ(expectedValue, *data) << "Proc: " << bulk.parallel_rank() << " entity: " << bulk.entity_key(entity) << " time: " << time << " data = " << *data << std::endl;
    }
  }
}

void create_n_hex_mesh_with_transient_field(int numElems, const std::string& fileName, unsigned timeSteps, const std::string& fieldName, stk::mesh::EntityRank rank)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::mesh::FieldBase& field = meta.declare_field<double>(rank, fieldName);

  stk::io::set_field_role(field, Ioss::Field::TRANSIENT);

  std::ostringstream os;
  os << "generated:1x1x" << numElems;

  stk::mesh::put_field_on_mesh(field, meta.universal_part(), 1, nullptr);
  stk::io::fill_mesh(os.str(), *bulk);

  stk::io::StkMeshIoBroker ioBroker;
  ioBroker.set_bulk_data(*bulk);
  unsigned dbIndex = ioBroker.create_output_mesh(fileName, stk::io::WRITE_RESULTS);
  ioBroker.add_field(dbIndex, field);
  ioBroker.write_output_mesh(dbIndex);
  for(unsigned i = 0; i <= timeSteps; i++) {
    double time = i;
    set_field_data(*bulk, field, time);
    ioBroker.begin_output_step(dbIndex, time);
    ioBroker.write_defined_output_fields(dbIndex);
    ioBroker.end_output_step(dbIndex);
  }
}

void create_n_hex_mesh_with_transient_field_and_global_data(int numElems, const std::string& fileName, unsigned timeSteps, const std::string& fieldName, const std::string& globalVarName, stk::mesh::EntityRank rank)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::mesh::FieldBase& field = meta.declare_field<double>(rank, fieldName);

  stk::io::set_field_role(field, Ioss::Field::TRANSIENT);

  std::ostringstream os;
  os << "generated:1x1x" << numElems;

  stk::mesh::put_field_on_mesh(field, meta.universal_part(), 1, nullptr);
  stk::io::fill_mesh(os.str(), *bulk);

  stk::io::StkMeshIoBroker ioBroker;
  ioBroker.set_bulk_data(*bulk);
  unsigned dbIndex = ioBroker.create_output_mesh(fileName, stk::io::WRITE_RESULTS);
  ioBroker.add_field(dbIndex, field);
  ioBroker.add_global(dbIndex, globalVarName, 1, Ioss::Field::DOUBLE);
  ioBroker.write_output_mesh(dbIndex);
  for(unsigned i = 0; i <= timeSteps; i++) {
    double time = i;
    double globalVar = time;
    set_field_data(*bulk, field, time);
    ioBroker.begin_output_step(dbIndex, time);
    ioBroker.write_defined_output_fields(dbIndex);
    ioBroker.write_global(dbIndex, globalVarName, globalVar);
    ioBroker.end_output_step(dbIndex);
  }
}

void test_transient_mesh_output(stk::io::StkMeshIoBroker& inputBroker, const std::string& fileName, const std::string& fieldName)
{
  std::shared_ptr<stk::mesh::BulkData> outputBulk = build_mesh(3, MPI_COMM_WORLD);
  stk::io::StkMeshIoBroker outputBroker;

  stk::io::fill_mesh_preexisting(outputBroker, fileName, *outputBulk);

  EXPECT_TRUE(are_bulk_data_equivalent(inputBroker.bulk_data(), *outputBulk));

  std::vector<double> inputTimeSteps = inputBroker.get_time_steps();
  std::vector<double> outputTimeSteps = outputBroker.get_time_steps();

  EXPECT_TRUE(inputTimeSteps == outputTimeSteps);
}

TEST(TransientFieldTransfer, writeTransientMesh)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) { return; }

  std::string inputFileName = "inputFile.e";
  std::string outputFileName = "outputFile.e";
  std::string fieldName = "field";
  unsigned numSteps = 1;
  int targetNumDomains = stk::parallel_machine_size(MPI_COMM_WORLD);
  int globalNumNodes = 12;
  int globalNumElems = 2;

  create_n_hex_mesh_with_transient_field(globalNumElems, inputFileName, numSteps, fieldName, stk::topology::NODE_RANK);

  stk::balance::m2n::OutputSerializerBulkData inputBulk(3, MPI_COMM_WORLD);
  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, inputFileName, inputBulk);
  int subdomain = inputBulk.parallel_rank();

  stk::io::EntitySharingInfo nodeSharingInfo = get_two_hex_mesh_node_sharing_info(inputBulk);

  stk::balance::m2n::TransientFieldTransferById m2nIo(ioBroker, targetNumDomains);
  m2nIo.setup_subdomain(inputBulk, outputFileName, subdomain, nodeSharingInfo, globalNumNodes, globalNumElems);
  m2nIo.get_subdomain_writer(subdomain).write_mesh();

  for(unsigned i = 0; i <= numSteps; i++) {
    double outputTime = (double)i;
    m2nIo.get_subdomain_writer(subdomain).write_transient_data(outputTime);
  }

  test_transient_mesh_output(ioBroker, outputFileName, fieldName);
  clean_up(targetNumDomains, inputFileName);
  clean_up(targetNumDomains, outputFileName);
}

TEST(TransientFieldTransfer, invalidSubdomainIndex)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3, MPI_COMM_WORLD);
  stk::io::StkMeshIoBroker ioBroker;
  ioBroker.set_bulk_data(*bulk);

  stk::balance::m2n::TransientFieldTransferById m2nIo(ioBroker, 2);
  EXPECT_THROW(m2nIo.get_subdomain_writer(0).write_mesh(), std::logic_error);
}

TEST(TransientFieldTransfer, setupSubdomainError)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;

  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, "generated:1x1x1", bulk);

  stk::io::EntitySharingInfo nodeSharingInfo;

  std::string fileName = "test_m2n_output.g";

  stk::balance::m2n::TransientFieldTransferById m2nIo(ioBroker, 1);
  EXPECT_THROW(m2nIo.get_subdomain_writer(bulk.parallel_rank()).write_mesh(), std::logic_error) << " setup_subdomain() was not called";
}

void load_time_step(stk::io::StkMeshIoBroker& ioBroker, double time)
{
  ioBroker.read_defined_input_fields(time);
}

void test_load_transient_data(stk::io::StkMeshIoBroker& ioBroker, const std::string& fieldName, stk::mesh::EntityRank fieldRank,
                              unsigned expectedNumSteps)
{
  const stk::mesh::BulkData& inputBulk = ioBroker.bulk_data();
  const stk::mesh::MetaData& inputMeta = inputBulk.mesh_meta_data();
  unsigned inputNumSteps = ioBroker.get_num_time_steps();
  EXPECT_EQ(inputNumSteps, expectedNumSteps+1);

  stk::mesh::FieldBase* field = inputMeta.get_field(fieldRank, fieldName);
  EXPECT_TRUE(field != nullptr);

  std::vector<double> inputTimeSteps = ioBroker.get_time_steps();
  for(unsigned i = 0; i <= expectedNumSteps; i++) {
    double time = inputTimeSteps[i];
    ioBroker.read_defined_input_fields(time);
    test_field_data(inputBulk, *field, time);
  }
}

TEST(TransientFieldTransfer, loadTransientData)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) { return; }

  std::string inputFileName = "inputFile.e";
  std::string fieldName = "field";
  unsigned numSteps = 1;
  stk::mesh::EntityRank fieldRank = stk::topology::NODE_RANK;

  create_n_hex_mesh_with_transient_field(2, inputFileName, numSteps, fieldName, fieldRank);

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::BulkData& inputBulk = *bulkPtr;
  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, inputFileName, inputBulk);

  test_load_transient_data(ioBroker, fieldName, fieldRank, numSteps);
  clean_up(inputBulk.parallel_size(), inputFileName);
}

TEST(TransientFieldTransfer, transferTransientData)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) { return; }

  std::string inputFileName = "inputFile.e";
  std::string fieldName = "field";
  unsigned numSteps = 1;
  std::string outputFileName = "outputFile.e";
  stk::mesh::EntityRank fieldRank = stk::topology::NODE_RANK;
  int globalNumNodes = 12;
  int globalNumElems = 2;

  create_n_hex_mesh_with_transient_field(globalNumElems, inputFileName, numSteps, fieldName, fieldRank);

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::BulkData& inputBulk = *bulkPtr;
  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, inputFileName, inputBulk);

  test_load_transient_data(ioBroker, fieldName, fieldRank, numSteps);

  stk::balance::m2n::OutputSerializerBulkData outputBulk(MPI_COMM_WORLD);
  int subdomain = inputBulk.parallel_rank();

  stk::tools::copy_mesh(inputBulk, inputBulk.mesh_meta_data().universal_part(), outputBulk);

  stk::io::EntitySharingInfo nodeSharingInfo = get_two_hex_mesh_node_sharing_info(inputBulk);
  stk::balance::m2n::TransientFieldTransferById m2nIo(ioBroker, inputBulk.parallel_size());
  m2nIo.setup_subdomain(outputBulk, outputFileName, subdomain, nodeSharingInfo, globalNumNodes, globalNumElems);
  m2nIo.get_subdomain_writer(subdomain).write_mesh();

  std::vector<double> inputTimeSteps = ioBroker.get_time_steps();
  stk::mesh::FieldBase* outputField = outputBulk.mesh_meta_data().get_field(fieldRank, fieldName);
  EXPECT_TRUE(outputField != nullptr);

  for(double time : inputTimeSteps) {
    load_time_step(ioBroker, time);
    m2nIo.transfer_transient_data(subdomain);
    test_field_data(outputBulk, *outputField, time);
  }
  clean_up(inputBulk.parallel_size(), inputFileName);
  clean_up(outputBulk.parallel_size(), outputFileName);
}

void test_transient_data_from_file(const std::string& filename, stk::mesh::EntityRank fieldRank, const std::string& fieldName, int numSubdomains)
{
  if(stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    for(int i = 0; i < numSubdomains; i++) {
      std::string subdomainName = stk::io::construct_filename_for_serial_or_parallel(filename, numSubdomains, i);

      std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_SELF);
      stk::mesh::BulkData& inputBulk = *bulkPtr;
      stk::io::StkMeshIoBroker ioBroker;
      stk::io::fill_mesh_preexisting(ioBroker, subdomainName, inputBulk);

      std::vector<double> inputTimeSteps = ioBroker.get_time_steps();
      stk::mesh::FieldBase* field = inputBulk.mesh_meta_data().get_field(fieldRank, fieldName);
      EXPECT_TRUE(field != nullptr);

      for(double time : inputTimeSteps) {
        load_time_step(ioBroker, time);
        test_field_data(inputBulk, *field, time);
      }
    }
  }
}

void test_transient_data_from_file(const std::string& filename, stk::mesh::EntityRank fieldRank, const std::string& fieldName)
{
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::BulkData& inputBulk = *bulkPtr;
  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, filename, inputBulk);

  std::vector<double> inputTimeSteps = ioBroker.get_time_steps();
  stk::mesh::FieldBase* field = inputBulk.mesh_meta_data().get_field(fieldRank, fieldName);
  EXPECT_TRUE(field != nullptr);

  for(double time : inputTimeSteps) {
    load_time_step(ioBroker, time);
    test_field_data(inputBulk, *field, time);
  }
}

void test_transient_and_global_data_from_file(const std::string& filename, stk::mesh::EntityRank fieldRank, const std::string& fieldName, const std::string& globalVarName)
{
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::BulkData& inputBulk = *bulkPtr;
  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, filename, inputBulk);

  std::vector<double> inputTimeSteps = ioBroker.get_time_steps();
  stk::mesh::FieldBase* field = inputBulk.mesh_meta_data().get_field(fieldRank, fieldName);
  EXPECT_TRUE(field != nullptr);
  EXPECT_TRUE(ioBroker.has_input_global(globalVarName));
  double globalVar;

  for(double time : inputTimeSteps) {
    load_time_step(ioBroker, time);
    test_field_data(inputBulk, *field, time);
    EXPECT_TRUE(ioBroker.get_global(globalVarName, globalVar));
    EXPECT_DOUBLE_EQ(time, globalVar);
  }
}

TEST(TransientFieldTransfer, testTransferAndWrite)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) { return; }

  std::string inputFileName = "inputFile.e";
  std::string fieldName = "field";
  unsigned numSteps = 2;
  std::string outputFileName = "outputFile.e";
  stk::mesh::EntityRank fieldRank = stk::topology::NODE_RANK;
  int globalNumNodes = 12;
  int globalNumElems = 2;

  create_n_hex_mesh_with_transient_field(globalNumElems, inputFileName, numSteps, fieldName, fieldRank);

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::BulkData& inputBulk = *bulkPtr;
  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, inputFileName, inputBulk);

  test_load_transient_data(ioBroker, fieldName, fieldRank, numSteps);

  stk::balance::m2n::OutputSerializerBulkData outputBulk(MPI_COMM_WORLD);
  int subdomain = inputBulk.parallel_rank();

  stk::tools::copy_mesh(inputBulk, inputBulk.mesh_meta_data().universal_part(), outputBulk);

  stk::io::EntitySharingInfo nodeSharingInfo = get_two_hex_mesh_node_sharing_info(inputBulk);
  stk::balance::m2n::TransientFieldTransferById m2nIo(ioBroker, inputBulk.parallel_size());

  m2nIo.setup_subdomain(outputBulk, outputFileName, subdomain, nodeSharingInfo, globalNumNodes, globalNumElems);
  m2nIo.transfer_and_write_transient_data(subdomain);

  test_transient_data_from_file(outputFileName, fieldRank, fieldName);
  clean_up(inputBulk.parallel_size(), inputFileName);
  clean_up(outputBulk.parallel_size(), outputFileName);
}

TEST(TransientFieldTransfer, testTransferAndWriteWithGlobalVar)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) { return; }

  std::string inputFileName = "inputFile.e";
  std::string fieldName = "field";
  std::string globalVarName = "global";
  unsigned numSteps = 2;
  std::string outputFileName = "outputFile.e";
  stk::mesh::EntityRank fieldRank = stk::topology::NODE_RANK;
  int globalNumNodes = 12;
  int globalNumElems = 2;

  create_n_hex_mesh_with_transient_field_and_global_data(globalNumElems, inputFileName, numSteps, fieldName, globalVarName, fieldRank);

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::BulkData& inputBulk = *bulkPtr;
  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, inputFileName, inputBulk);

  test_load_transient_data(ioBroker, fieldName, fieldRank, numSteps);

  stk::balance::m2n::OutputSerializerBulkData outputBulk(MPI_COMM_WORLD);
  int subdomain = inputBulk.parallel_rank();

  stk::tools::copy_mesh(inputBulk, inputBulk.mesh_meta_data().universal_part(), outputBulk);

  stk::io::EntitySharingInfo nodeSharingInfo = get_two_hex_mesh_node_sharing_info(inputBulk);
  stk::balance::m2n::TransientFieldTransferById m2nIo(ioBroker, inputBulk.parallel_size());

  m2nIo.setup_subdomain(outputBulk, outputFileName, subdomain, nodeSharingInfo, globalNumNodes, globalNumElems);
  m2nIo.transfer_and_write_transient_data(subdomain);

  test_transient_and_global_data_from_file(outputFileName, fieldRank, fieldName, globalVarName);
  clean_up(inputBulk.parallel_size(), inputFileName);
  clean_up(outputBulk.parallel_size(), outputFileName);
}

TEST(TransientFieldTransfer, testMockBalance2x4)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) { return; }

  std::string inputFileName = "inputFile.e";
  std::string fieldName = "field";
  unsigned numSteps = 5;
  std::string outputFileName = "outputFile.e";
  stk::mesh::EntityRank fieldRank = stk::topology::NODE_RANK;

  create_n_hex_mesh_with_transient_field(4, inputFileName, numSteps, fieldName, fieldRank);

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::BulkData& inputBulk = *bulkPtr;
  stk::mesh::MetaData& inputMeta = inputBulk.mesh_meta_data();
  stk::mesh::PartVector partVector(4);
  partVector[0] = &inputMeta.declare_part("element_1_part", stk::topology::ELEMENT_RANK);
  partVector[1] = &inputMeta.declare_part("element_2_part", stk::topology::ELEMENT_RANK);
  partVector[2] = &inputMeta.declare_part("element_3_part", stk::topology::ELEMENT_RANK);
  partVector[3] = &inputMeta.declare_part("element_4_part", stk::topology::ELEMENT_RANK);

  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, inputFileName, inputBulk);

  stk::mesh::EntityIdVector entityIdVec;
  if(inputBulk.parallel_size() == 1) {
    entityIdVec = {1, 2, 3, 4};
  } else {
    if(inputBulk.parallel_rank() == 0) {
      entityIdVec = {1, 2};
    } else {
      entityIdVec = {3, 4};
    }
  }

  inputBulk.modification_begin();

  for(stk::mesh::EntityId id : entityIdVec) {
    stk::mesh::Entity entity = inputBulk.get_entity(stk::topology::ELEMENT_RANK, id);
    EXPECT_TRUE(inputBulk.is_valid(entity));
    int subdomain = id-1;
    inputBulk.change_entity_parts(entity, stk::mesh::PartVector{ partVector[subdomain] }, stk::mesh::PartVector{});
  }

  inputBulk.modification_end();

  test_load_transient_data(ioBroker, fieldName, fieldRank, numSteps);
  int targetNumDomains = 4;
  int globalNumNodes = 20;
  int globalNumElems = 4;

  for(stk::mesh::EntityId id : entityIdVec) {
    int subdomain = id-1;

    stk::balance::m2n::OutputSerializerBulkData outputBulk(MPI_COMM_SELF);
    stk::tools::copy_mesh(inputBulk, *partVector[subdomain], outputBulk);

    stk::io::EntitySharingInfo nodeSharingInfo = get_four_hex_mesh_node_sharing_info(inputBulk, subdomain);
    stk::balance::m2n::TransientFieldTransferById m2nIo(ioBroker, targetNumDomains);

    m2nIo.setup_subdomain(outputBulk, outputFileName, subdomain, nodeSharingInfo, globalNumNodes, globalNumElems);
    m2nIo.transfer_and_write_transient_data(subdomain);
  }

  test_transient_data_from_file(outputFileName, fieldRank, fieldName, targetNumDomains);
  clean_up(inputBulk.parallel_size(), inputFileName);
  clean_up(targetNumDomains, outputFileName);
}

}
