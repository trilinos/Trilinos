#include <gtest/gtest.h>                // for TEST
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for BulkData
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>


namespace
{

void test_that_ids_are_unique(const stk::mesh::BulkData &bulkData, stk::mesh::EntityRank rank, std::vector<stk::mesh::EntityId>& requestedIds)
{
  std::vector<stk::mesh::EntityId> ids_in_use;
  stk::mesh::for_each_entity_run(bulkData, rank, bulkData.mesh_meta_data().locally_owned_part(),
    [&](const stk::mesh::BulkData& bulk, stk::mesh::Entity entity) {
      ids_in_use.push_back(bulkData.identifier(entity));
    });

  // ========

  std::sort(requestedIds.begin(), requestedIds.end());
  std::vector<stk::mesh::EntityId>::iterator iter1 = std::unique(requestedIds.begin(), requestedIds.end());
  STK_ThrowRequireMsg(iter1 == requestedIds.end(), "Ids not unique. " << __FILE__ << __LINE__);

  stk::CommSparse comm(bulkData.parallel());

  for(int phase = 0; phase < 2; ++phase)
  {
    for(int i = 0; i < bulkData.parallel_size(); ++i)
    {
      if(i != bulkData.parallel_rank())
      {
        for(size_t j = 0; j < requestedIds.size(); ++j)
        {
          bool is_id_unique = std::binary_search(ids_in_use.begin(), ids_in_use.end(), requestedIds[j]);
          STK_ThrowRequireMsg(is_id_unique == false, "ID="<<requestedIds[j]<<" already in use in the mesh." <<  __FILE__<< __LINE__);
          comm.send_buffer(i).pack<uint64_t>(requestedIds[j]);
        }
      }
    }

    if(phase == 0)
    {
      comm.allocate_buffers();
    }
    else
    {
      comm.communicate();
    }
  }

  for(int i = 0; i < bulkData.parallel_size(); ++i) {
    if(i != bulkData.parallel_rank()) {
      while(comm.recv_buffer(i).remaining()) {
        uint64_t id;
        comm.recv_buffer(i).unpack(id);
        bool is_other_procs_id_on_this_proc = std::binary_search(requestedIds.begin(), requestedIds.end(), id);
        STK_ThrowRequireMsg(is_other_procs_id_on_this_proc == false, "Id requested on proc " << i << " also requested on proc " << bulkData.parallel_rank()<< ". " << __FILE__<< __LINE__);
        bool is_id_already_in_use = std::binary_search(ids_in_use.begin(), ids_in_use.end(), id);
        STK_ThrowRequireMsg(is_id_already_in_use == false, "Id requested on proc " << i << " already in use on proc " << bulkData.parallel_rank() << ". " << __FILE__ << __LINE__);
      }
    }
  }
}

//BEGIN_TEST_1
TEST(StkMeshHowTo, use_generate_new_ids)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  int num_procs = stk::parallel_machine_size(communicator);
  std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(communicator).create();

  const std::string generatedMeshSpecification = "generated:1x1x" + std::to_string(num_procs);
  stk::io::fill_mesh(generatedMeshSpecification, *bulkPtr);

  // Given a mesh, request 10 unique node ids

  std::vector<stk::mesh::EntityId> requestedIds;
  unsigned numRequested = 10;

  bulkPtr->generate_new_ids(stk::topology::NODE_RANK, numRequested, requestedIds);

  test_that_ids_are_unique(*bulkPtr, stk::topology::NODE_RANK, requestedIds);
}
//END_TEST_1

}

