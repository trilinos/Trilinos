#include "gtest/gtest.h"
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_unit_test_utils/timer.hpp>

namespace
{

void build_mesh_in_batches(stk::mesh::BulkData& mesh,
                           size_t numElemsPerProc,
                           size_t batchSize)
{
  unsigned numOwnedNodes = (numElemsPerProc+1)*2;
  unsigned numSharedNodes = (numElemsPerProc+1)*2;

  std::vector<stk::mesh::EntityId> ownedNodes(numOwnedNodes);
  std::vector<stk::mesh::EntityId> sharedNodes(numSharedNodes);

  int localProc = mesh.parallel_rank();
  int numProcs = mesh.parallel_size();
  stk::mesh::EntityId firstLocalOwnedNode = 1 + numOwnedNodes*localProc;
  stk::mesh::EntityId firstLocalSharedNode = 1 + numOwnedNodes*numProcs;

  for(size_t inode=0; inode<numOwnedNodes; ++inode) {
    ownedNodes[inode] = firstLocalOwnedNode + inode;
    sharedNodes[inode] = firstLocalSharedNode + inode;
  }

  size_t ielem = 0;
  stk::mesh::EntityId firstLocalElemId = 1 + numElemsPerProc*localProc;
  constexpr size_t numNodesPerElem = 8;
  constexpr size_t numNodesPerSide = 4;
  stk::mesh::EntityIdVector nodeIds(numNodesPerElem);

  stk::mesh::Part& hexPart = mesh.mesh_meta_data().get_topology_root_part(stk::topology::HEX_8);
  int otherProc = 1 - localProc;
  while(ielem < numElemsPerProc) {
    mesh.modification_begin();
    stk::mesh::EntityProcVec sharedNodesAndProcs;

    for(size_t i=0; i<batchSize; ++i) {
      stk::mesh::EntityId elemId = firstLocalElemId + ielem;
      unsigned startNodeIdx = ielem*2;
      for(size_t n=0; n<numNodesPerSide; ++n) {
        nodeIds[n] = ownedNodes[startNodeIdx+n];
      }
      for(size_t n=0; n<numNodesPerSide; ++n) {
        nodeIds[numNodesPerSide+n] = sharedNodes[startNodeIdx+n];
      }

      stk::mesh::declare_element(mesh, hexPart, elemId, nodeIds);
      for(size_t n=0; n<numNodesPerSide; ++n) {
        stk::mesh::Entity sharedNode = mesh.get_entity(stk::topology::NODE_RANK, sharedNodes[startNodeIdx+n]);
        sharedNodesAndProcs.push_back(stk::mesh::EntityProc(sharedNode, otherProc));
      }
      ++ielem;
    }
    mesh.add_node_sharing(sharedNodesAndProcs);

    mesh.modification_end();
  }
}

TEST(AddNodeSharing, Timing)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }

  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();

  constexpr unsigned numRuns = 5;
#ifndef NDEBUG
  constexpr size_t numElemsPerProc = 2000;
#else
  constexpr size_t numElemsPerProc = 20000;
#endif
  constexpr size_t batchSize = 100;

  for(unsigned run=0; run<numRuns; ++run) {
    batchTimer.start_batch_timer();

    std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
    build_mesh_in_batches(*bulk, numElemsPerProc, batchSize);

    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(numRuns);
}

}

