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
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/ExodusTranslator.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_unit_test_utils/timer.hpp>
#include <stk_performance_tests/stk_mesh/multi_block.hpp>

#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/NgpFieldParallel.hpp>
#include <vector>
#include <string>

#include <stk_unit_test_utils/getOption.h>

namespace
{

double initial_value1[3] = {-1, 2, -0.3};

void createNodalVectorField(stk::mesh::MetaData& meshMetaData,
                            const std::string &fieldName,
                            const std::string &partName)
{
    stk::mesh::Field<double> &field1 = meshMetaData.declare_field<double>(stk::topology::NODE_RANK, fieldName);
    stk::mesh::Part* part = meshMetaData.get_part(partName);
    ASSERT_TRUE(part != nullptr)<<partName<<" not found";
    stk::mesh::put_field_on_mesh(field1, *part, 3, initial_value1);
}

void createNodalVectorFields(stk::mesh::MetaData& meshMetaData,
                             const std::string& fieldBaseName,
                             unsigned numFields,
                             const std::string& partName)
{
  for(unsigned i=0; i<numFields; ++i) {
    std::string fieldName = fieldBaseName+std::to_string(i);
    createNodalVectorField(meshMetaData, fieldName, partName);
  }
}

size_t setFieldData(stk::mesh::BulkData& mesh, stk::mesh::Selector select)
{
    const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
    const stk::mesh::FieldVector& fields = meta.get_fields(stk::topology::NODE_RANK);
    const stk::mesh::BucketVector& node_buckets = mesh.get_buckets(stk::topology::NODE_RANK, select);
    size_t num_nodes = 0;
    for(size_t i=0; i<node_buckets.size(); ++i) {
        stk::mesh::Bucket& bucket = *node_buckets[i];
        num_nodes += bucket.size();
        for(size_t f=0; f<fields.size(); ++f) {
            double* data = static_cast<double*>(stk::mesh::field_data(*fields[f], bucket));
            unsigned flen = stk::mesh::field_scalars_per_entity(*fields[f], bucket);
            for(size_t n=0; n<flen*bucket.size(); ++n) {
                data[n] = static_cast<double>(f);
            }
        }
    }

    return num_nodes;
}

void createMetaAndBulkData(stk::io::StkMeshIoBroker &exodusFileReader,
                           const std::string& genMeshSpec,
                           unsigned numBlocks,
                           unsigned numFields)
{
    std::string exodusFileName = stk::unit_test_util::get_option("-i", "NO_FILE_SPECIFIED");
    if (exodusFileName == "NO_FILE_SPECIFIED") {
      exodusFileName = genMeshSpec;
    }

    exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);

    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int myProc = stk::parallel_machine_rank(communicator);

    if (myProc==0) {
      std::cerr << "Starting To Read Mesh: " << exodusFileName << std::endl;
    }

    exodusFileReader.create_input_mesh();
    std::shared_ptr<stk::mesh::MetaData> stkMeshMetaData = exodusFileReader.meta_data_ptr();
    if (numBlocks==1) {
      createNodalVectorFields(*stkMeshMetaData, "nodalField", numFields, "{universal}");
    }
    else {
      stk::performance_tests::setup_multiple_blocks(*stkMeshMetaData, numBlocks);
      stk::mesh::PartVector blockParts;
      stk::mesh::fill_element_block_parts(*stkMeshMetaData, stk::topology::HEX_8, blockParts);
      for(const stk::mesh::Part* part : blockParts) {
        std::string fieldBaseName = part->name()+"nodalField";
        createNodalVectorFields(*stkMeshMetaData, fieldBaseName, numFields, part->name());
      }
    }
    stkMeshMetaData->commit();

    std::shared_ptr<stk::mesh::BulkData> arg_bulk_data = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_aura_option(stk::mesh::BulkData::AUTO_AURA).create(stkMeshMetaData);
    exodusFileReader.set_bulk_data(arg_bulk_data);
    stk::mesh::BulkData& stkMeshBulkData = *arg_bulk_data;

    bool delay_field_data_allocation = true;
    exodusFileReader.populate_mesh(delay_field_data_allocation);
    exodusFileReader.populate_field_data();

    size_t num_ghostings = stkMeshBulkData.ghostings().size();
    std::vector<size_t> num_ghost_nodes(num_ghostings,0);
    for (size_t ghost_i = 0 ; ghost_i < num_ghostings ; ++ghost_i) {
        num_ghost_nodes[ghost_i] = setFieldData(stkMeshBulkData, stkMeshBulkData.ghosting_part(*stkMeshBulkData.ghostings()[ghost_i]));
    }

    for(int i=0; i<stkMeshBulkData.parallel_size(); ++i) {
        stk::parallel_machine_barrier(MPI_COMM_WORLD);
        if (stkMeshBulkData.parallel_rank() == i) {
            for (size_t ghost_i = 0 ; ghost_i<num_ghostings ; ++ghost_i) {
                std::cerr <<"proc "<<stkMeshBulkData.parallel_rank()<<", Number of " << stkMeshBulkData.ghostings()[ghost_i]->name() << " nodes: " << num_ghost_nodes[ghost_i] << std::endl;
            }
        }

        stk::parallel_machine_barrier(MPI_COMM_WORLD);
        std::cerr << std::flush;
        stk::parallel_machine_barrier(MPI_COMM_WORLD);
    }
    if (myProc==0) {
      std::cerr << "Finished Reading Mesh" << std::endl;
    }

    stk::mesh::Selector allEntities = stkMeshMetaData->universal_part();
    std::vector<size_t> entityCounts;
    stk::mesh::count_entities(allEntities, stkMeshBulkData, entityCounts);
    size_t numElements = entityCounts[stk::topology::ELEMENT_RANK];
    size_t numNodes = entityCounts[stk::topology::NODE_RANK];

    if (myProc==0) {
      std::cerr << "Number of elements: " << numElements << std::endl;
      std::cerr << "Number of nodes: " << numNodes << std::endl;
    }
}

void createNgpFields(stk::mesh::BulkData& bulk)
{
  auto stkFields = bulk.mesh_meta_data().get_fields();
  for(auto stkField : stkFields) {
    stk::mesh::get_updated_ngp_field<double>(*stkField);
  }
}

template<typename T>
void set_modify_on_device(std::vector<stk::mesh::NgpField<T>*> ngpFields)
{
  for(auto ngpField : ngpFields) {
    ngpField->modify_on_device();
  }
}

void addPartToGhosting(stk::mesh::BulkData & bulk, const std::string & partName, stk::mesh::Ghosting& ghost)
{
    const int numProcs = bulk.parallel_size();
    const int myProc = bulk.parallel_rank();
    stk::mesh::MetaData &meta = bulk.mesh_meta_data();
    std::vector<stk::mesh::EntityProc> entities_to_ghost;
    stk::mesh::Selector surface_selector = *meta.get_part(partName) & meta.locally_owned_part();
    stk::mesh::EntityVector surface_nodes;
    stk::mesh::get_selected_entities(surface_selector,bulk.buckets(stk::topology::NODE_RANK), surface_nodes);
    const int next_proc = (myProc+1)%numProcs;
    for (size_t i=0 ; i<surface_nodes.size() ; ++i) {
        entities_to_ghost.push_back(stk::mesh::EntityProc(surface_nodes[i],next_proc));
    }
    bulk.change_ghosting( ghost, entities_to_ghost );
}

void test_communicate_field_data_all_ghosting(stk::ParallelMachine communicator, int num_iters)
{
    stk::unit_test_util::BatchTimer batchTimer(communicator);
    batchTimer.initialize_batch_timer();

    stk::io::StkMeshIoBroker exodusFileReader(communicator);

    std::string genMeshSpec = "generated:60x60x48|sideset:xXyY";
    const unsigned numBlocks = 1;
    const unsigned numFields = 8;
    createMetaAndBulkData(exodusFileReader,genMeshSpec, numBlocks, numFields);
    stk::mesh::BulkData &mesh = exodusFileReader.bulk_data();

    mesh.modification_begin();
    stk::mesh::Ghosting & ghosting_1 = mesh.create_ghosting( "CUSTOM_1" );
    stk::mesh::Ghosting & ghosting_2 = mesh.create_ghosting( "CUSTOM_2" );
    stk::mesh::Ghosting & ghosting_3 = mesh.create_ghosting( "CUSTOM_3" );
    stk::mesh::Ghosting & ghosting_4 = mesh.create_ghosting( "CUSTOM_4" );

    addPartToGhosting(mesh, "surface_1", ghosting_1);
    addPartToGhosting(mesh, "surface_2", ghosting_2);
    addPartToGhosting(mesh, "surface_3", ghosting_3);
    addPartToGhosting(mesh, "surface_4", ghosting_4);
    mesh.modification_end();

    size_t num_ghostings = mesh.ghostings().size();
    std::ostringstream oss;
    for (size_t ghost_i = 0 ; ghost_i < num_ghostings ; ++ghost_i) {
        size_t num_nodes = stk::mesh::count_selected_entities(mesh.ghosting_part(*mesh.ghostings()[ghost_i]),mesh.buckets(stk::topology::NODE_RANK));
        oss <<"proc "<<mesh.parallel_rank()<<", Number of " << mesh.ghostings()[ghost_i]->name() << " nodes: " << num_nodes << std::endl;
    }

    std::cerr << oss.str() << std::endl;

    const int my_proc = mesh.parallel_rank();
    if (my_proc == 0) {
        std::cerr << "Calling communicate_field_data " << num_iters << " times"<<std::endl;
    }

    const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
    const stk::mesh::FieldVector& fields = meta.get_fields();
    std::vector<const stk::mesh::FieldBase*> const_fields(fields.size());
    for(size_t i=0; i<fields.size(); ++i) {
        const_fields[i] = fields[i];
    }

    const unsigned NUM_RUNS = 5;
  
    for (unsigned j = 0; j < NUM_RUNS; j++) {
      stk::parallel_machine_barrier(mesh.parallel());
      batchTimer.start_batch_timer();

      for(int iter=0; iter<num_iters; ++iter) {
          for (size_t ghost_i = 0 ; ghost_i<mesh.ghostings().size() ; ++ghost_i) {
              stk::mesh::communicate_field_data(*mesh.ghostings()[ghost_i], const_fields);
          }
      }

      for(int iter=0; iter<num_iters; ++iter) {
          stk::mesh::communicate_field_data(mesh, const_fields);
      }

      batchTimer.stop_batch_timer();
    }
    batchTimer.print_batch_timing(num_iters);
}

void test_communicate_field_data_ghosting(MPI_Comm communicator,
                                          unsigned numBlocks, unsigned numFields, int num_iters)
{
    stk::unit_test_util::BatchTimer batchTimer(communicator);
    batchTimer.initialize_batch_timer();

    stk::io::StkMeshIoBroker exodusFileReader(communicator);
    std::string genMeshSpec = "generated:100x100x48|sideset:xXyY";
    createMetaAndBulkData(exodusFileReader,genMeshSpec, numBlocks, numFields);

    stk::mesh::BulkData& mesh = exodusFileReader.bulk_data();
    const int my_proc = mesh.parallel_rank();
    if (my_proc == 0) {
        std::cerr << "Calling communicate_field_data " << num_iters << " times"<<std::endl;
    }

    const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
    const stk::mesh::FieldVector& fields = meta.get_fields();
    std::vector<const stk::mesh::FieldBase*> const_fields(fields.size());
    for(size_t i=0; i<fields.size(); ++i) {
        const_fields[i] = fields[i];
    }

    const unsigned NUM_RUNS = 5;
  
    for (unsigned j = 0; j < NUM_RUNS; j++) {
      stk::parallel_machine_barrier(mesh.parallel());
      batchTimer.start_batch_timer();

      for(int iter=0; iter<num_iters; ++iter) {
        stk::mesh::communicate_field_data(mesh, const_fields);
      }
      batchTimer.stop_batch_timer();
    }
    batchTimer.print_batch_timing(num_iters);
}

void test_communicate_field_data_ngp_ghosting(int num_iters, bool syncToHostEveryIter)
{
  const int meshDim = stk::unit_test_util::get_command_line_option("-s", 100);
  std::string meshDimStr = std::to_string(meshDim);
  std::string meshSpec = "generated:" + meshDimStr + "x" + meshDimStr + "x" + meshDimStr;

  stk::ParallelMachine comm = MPI_COMM_WORLD;
  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  stk::io::StkMeshIoBroker exodusFileReader(comm);

  const unsigned numBlocks = 1;
  const unsigned numFields = 8;
  createMetaAndBulkData(exodusFileReader, meshSpec, numBlocks, numFields);

  stk::mesh::BulkData &mesh = exodusFileReader.bulk_data();
  const int my_proc = mesh.parallel_rank();

  if (my_proc == 0) {
      std::cerr << "Calling communicate_field_data " << num_iters << " times"<<std::endl;
  }

  createNgpFields(mesh);
  const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
  const stk::mesh::FieldVector& fields = meta.get_fields();
  std::vector<stk::mesh::NgpField<double>*> ngpFields;
  for(auto field : fields) {
    ngpFields.push_back(&stk::mesh::get_updated_ngp_field<double>(*field));
  }

  const unsigned NUM_RUNS = 5;
  
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    stk::parallel_machine_barrier(mesh.parallel());
    batchTimer.start_batch_timer();

    for(int iter=0; iter<num_iters; ++iter) {
      if(syncToHostEveryIter) {
        set_modify_on_device<double>(ngpFields);
      }

      stk::mesh::communicate_field_data(mesh, ngpFields, true);
    }
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(num_iters);
}

TEST(CommunicateFieldData, copy_to_all)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) < 2) { GTEST_SKIP(); }

  test_communicate_field_data_all_ghosting(comm, 300);
}


TEST(CommunicateFieldData, Ghosting)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs < 2) {
    return;
  }

  const unsigned numBlocks = 1;
  const unsigned numFields = 8;
  test_communicate_field_data_ghosting(communicator, numBlocks, numFields, 1000);
}

TEST(CommunicateFieldData, Ghosting_MultiBlock)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs < 2) {
    return;
  }

  const unsigned numBlocks = 1;
  const unsigned numFields = 8;
  test_communicate_field_data_ghosting(communicator, numBlocks, numFields, 300);
}

TEST(CommunicateFieldData, NgpGhosting)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) < 2) { GTEST_SKIP(); }

  int iter = stk::unit_test_util::get_command_line_option("-t", 1000);
  bool syncToHostEveryIter = stk::unit_test_util::get_command_line_option("-h", false);
  test_communicate_field_data_ngp_ghosting(iter, syncToHostEveryIter);
}

}
