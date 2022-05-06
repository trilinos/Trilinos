#include "test_stk_io.hpp"

#include <stk_topology/topology.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <stk_io/FillMesh.hpp>

namespace test_stk_lib {

void test_stk_io(stk::ParallelMachine comm, const std::string& meshSource, bool useAutoDecomp)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).create();
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  if (bulk->parallel_rank() == 0) {
    std::cout << "test_stk_io: meshSource="<<meshSource<<std::endl;
  }

  if (useAutoDecomp) {
    stk::io::fill_mesh_with_auto_decomp(meshSource, *bulk);
  }
  else {
    stk::io::fill_mesh(meshSource, *bulk);
  }

  std::vector<size_t> meshCounts(meta.entity_rank_count());
  stk::mesh::comm_mesh_counts(*bulk, meshCounts);

  if (bulk->parallel_rank() == 0) {
    std::cout<<"   global number of elements = "<<meshCounts[stk::topology::ELEM_RANK]<<std::endl;
  }
}

}

