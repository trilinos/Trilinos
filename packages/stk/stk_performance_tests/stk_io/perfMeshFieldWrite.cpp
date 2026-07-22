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

#include <gtest/gtest.h>
#include <stk_util/stk_config.h>
#include <stk_unit_test_utils/getOption.h>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_unit_test_utils/timer.hpp>

#include <unistd.h>

namespace stk_perf_io_mesh_field_write
{

void create_nodal_vector_field(stk::mesh::MetaData& meshMetaData,
                               const std::string &fieldName,
                               const std::string &partName)
{
  double initial_value[3] = {0, 0, 0};
  stk::mesh::Field<double> &field = meshMetaData.declare_field<double>(stk::topology::NODE_RANK, fieldName);
  stk::io::set_field_role(field, Ioss::Field::TRANSIENT);
  stk::mesh::Part* part = meshMetaData.get_part(partName);
  ASSERT_TRUE(part != nullptr)<<partName<<" not found";
  stk::mesh::put_field_on_mesh(field, *part, 3, initial_value);
}

void create_nodal_vector_fields(stk::mesh::MetaData& meshMetaData,
                                const std::string& fieldBaseName,
                                unsigned numFields,
                                const std::string& partName)
{
  for(unsigned i=0; i<numFields; ++i) {
    std::string fieldName = fieldBaseName+std::to_string(i);
    create_nodal_vector_field(meshMetaData, fieldName, partName);
  }
}

size_t set_field_data(stk::mesh::BulkData& mesh, const stk::mesh::Selector& select, int outputStep)
{
  const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
  const stk::mesh::FieldVector& fields = meta.get_fields(stk::topology::NODE_RANK);
  const stk::mesh::BucketVector& nodeBuckets = mesh.get_buckets(stk::topology::NODE_RANK, select);

  for (unsigned fieldIndex = 0; fieldIndex < fields.size(); ++fieldIndex) {
    const stk::mesh::FieldBase& field = *fields[fieldIndex];
    auto fieldData = field.data<double, stk::mesh::ReadWrite>();
    for (stk::mesh::Bucket* bucket : nodeBuckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entityIdx : bucket->entities()) {
        for (stk::mesh::ScalarIdx scalar : bucketValues.scalars()) {
          bucketValues(entityIdx, scalar) = static_cast<double>(fieldIndex + 10*outputStep);
        }
      }
    }
  }

  size_t numNodes = 0;
  for (stk::mesh::Bucket* bucket : nodeBuckets) {
    numNodes += bucket->size();
  }

  return numNodes;
}

void setup_mesh(const std::string& meshSpec, stk::mesh::BulkData& bulkData, unsigned numFields)
{
  stk::ParallelMachine communicator = bulkData.parallel();
  stk::io::StkMeshIoBroker exodusFileReader(communicator);

  exodusFileReader.set_bulk_data(bulkData);
  exodusFileReader.add_mesh_database(meshSpec, stk::io::READ_MESH);

  int myProc = bulkData.parallel_rank();

  if (myProc==0) {
    std::cerr << "Starting To Read Mesh: " << std::endl;
  }

  exodusFileReader.create_input_mesh();
  std::shared_ptr<stk::mesh::MetaData> metaData = exodusFileReader.meta_data_ptr();

  create_nodal_vector_fields(*metaData, "nodalField", numFields, "block_1");

  metaData->commit();

  bool delay_field_data_allocation = true;
  exodusFileReader.populate_mesh(delay_field_data_allocation);
  exodusFileReader.populate_field_data();
}

void add_output_fields(stk::io::StkMeshIoBroker& outputBroker, size_t outputFileIndex, stk::mesh::BulkData& bulkData)
{
  const stk::mesh::FieldVector fields = bulkData.mesh_meta_data().get_fields();
  for(stk::mesh::FieldBase* field : fields)
  {
    const Ioss::Field::RoleType* fieldRole = stk::io::get_field_role(*field);
    if(*fieldRole == Ioss::Field::TRANSIENT)
      outputBroker.add_field(outputFileIndex, *field);
  }
}

void do_output(stk::io::StkMeshIoBroker& outputBroker, size_t outputFileIndex,
               stk::mesh::BulkData& bulkData, unsigned outputStep,
               bool DO_FLUSH = false)
{
  set_field_data(bulkData, bulkData.mesh_meta_data().universal_part(), outputStep);

  double outputTime = static_cast<double>(outputStep);
  outputBroker.begin_output_step(outputFileIndex, outputTime);
  outputBroker.write_defined_output_fields(outputFileIndex);
  outputBroker.end_output_step(outputFileIndex);

  if(DO_FLUSH) {
    outputBroker.flush_output();
  }
}

TEST(StkIo, meshFieldWrite_hex_noAura)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int nProc = stk::parallel_machine_size(communicator);
  int myProc = stk::parallel_machine_rank(communicator);

  if (nProc > 16) { GTEST_SKIP(); }

  const unsigned NUM_RUNS = 1;
  const unsigned NUM_ITERS = 500;
  const unsigned NUM_FIELDS = 15;

  int ELEMS_PER_DIM = stk::unit_test_util::get_command_line_option("--ne", 80);
  bool DO_FLUSH = stk::unit_test_util::has_option("--flush");

  std::string elems = std::to_string(ELEMS_PER_DIM);
  std::string meshSpec = "generated:"+elems+"x"+elems+"x"+elems;
  std::string outputFileName = "meshFieldWrite.e";
  std::string perProcOutputFileName = stk::io::construct_filename_for_serial_or_parallel(outputFileName, nProc, myProc);

  if (myProc == 0) {
    std::cout << "Using mesh-spec: " << meshSpec << std::endl;
  }

  std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(communicator)
                                                              .set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA)
                                                              .create();
  setup_mesh(meshSpec, *bulkPtr, NUM_FIELDS);

  stk::unit_test_util::BatchTimer batchTimer(communicator);

  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    unlink(perProcOutputFileName.c_str());

    stk::parallel_machine_barrier(communicator);

    if (myProc == 0) {
      std::cout << "Starting run: " << j << std::endl;
    }

    stk::io::StkMeshIoBroker outputBroker(communicator);
    outputBroker.set_bulk_data(*bulkPtr);
    size_t outputFileIndex = outputBroker.create_output_mesh(outputFileName, stk::io::WRITE_RESULTS);

    batchTimer.start_batch_timer();

    add_output_fields(outputBroker, outputFileIndex, *bulkPtr);

    outputBroker.write_output_mesh(outputFileIndex);

    for(unsigned i=0; i<NUM_ITERS; ++i) {
      do_output(outputBroker, outputFileIndex, *bulkPtr, i,DO_FLUSH);
    }

    batchTimer.stop_batch_timer();

    unlink(perProcOutputFileName.c_str());
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

}
