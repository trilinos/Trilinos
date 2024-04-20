#include "test_stk_tools.hpp"

#include <iostream>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_tools/mesh_clone/MeshClone.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocks.hpp>
#include <stk_tools/block_extractor/ExtractBlocks.hpp>
#include <stk_tools/transfer_utils/TransientFieldTransferById.hpp>

namespace test_stk_lib {

void test_stk_tools()
{
  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    std::cout << "The stk_tools test says 'Hello, World!'" << std::endl;
  }
}

}

