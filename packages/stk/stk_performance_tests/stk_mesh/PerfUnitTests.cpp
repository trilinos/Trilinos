#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/perf_util.hpp>
#include <stk_util/util/AllocatorMemoryUsage.hpp>
#include <stk_util/util/memory_util.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>
#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>

#include <stk_mesh/base/BoundaryAnalysis.hpp>
#include <stk_mesh/base/SkinMesh.hpp>

#include <sstream>
#include <valgrind/callgrind.h>

namespace {

using namespace stk::mesh;

void print_memory_sum_all_procs(stk::ParallelMachine pm)
{
  const size_t p_rank = stk::parallel_machine_rank(pm);
  const size_t my_peak = stk::allocator_memory_usage<void>::peak_memory();
  size_t peak_sum = 0;
  int err = MPI_Reduce((void*)&my_peak, (void*)&peak_sum, 1 /*size*/, MPI_LONG_LONG, MPI_SUM, 0 /*root*/, pm);
  ThrowRequire(err == MPI_SUCCESS);

  if (p_rank == 0) {
    std::cout << "\nSTKPERF peak memory sum: " << peak_sum << std::endl;
  }
}

void check_valgrind_version()
{
  STKUNIT_ASSERT_EQ(__VALGRIND_MAJOR__, 3);
  STKUNIT_ASSERT_EQ(__VALGRIND_MINOR__, 8);
}

void print_debug_skip(stk::ParallelMachine pm)
{
#ifndef NDEBUG
  // We're in debug; need to tell test script not to validate cycle count
  const size_t p_rank = stk::parallel_machine_rank(pm);
  if (p_rank == 0) {
    std::cout << "\nSTKPERF SKIP VALIDATION" << std::endl;
  }
#endif
}

STKUNIT_UNIT_TEST( stk_mesh_perf_unit_test, induced_part )
{
  check_valgrind_version();

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const size_t p_size = stk::parallel_machine_size(pm);

  // serial only
  ThrowRequire(p_size == 1);

  const int x_dim = 500;
  const int y_dim = 3;
  const int z_dim = 3;

  // Set up dead-simple mesh, 2 quads sharing 2 nodes
  stk::mesh::fixtures::HexFixture mesh(pm, x_dim, y_dim, z_dim);
  MetaData & meta = mesh.m_meta;
  BulkData & bulk = mesh.m_bulk_data;

  const int num_inducable_parts = 4;
  const int num_non_inducable_parts = 10;

  PartVector node_parts;
  PartVector elem_parts;
  for (int i = 0; i < num_inducable_parts + num_non_inducable_parts; ++i) {
    std::ostringstream elem_part_name, node_part_name;
    std::string inducable_str = i < num_inducable_parts ? "Inducable" : "Non-inducable";
    elem_part_name << inducable_str << " Elem Part " << i;
    node_part_name << inducable_str << " Node Part " << i;
    if (i < num_inducable_parts) {
      elem_parts.push_back(&meta.declare_part(elem_part_name.str(), stk::topology::ELEMENT_RANK));
      node_parts.push_back(&meta.declare_part(node_part_name.str(), stk::topology::NODE_RANK));
    }
    else {
      elem_parts.push_back(&meta.declare_part(elem_part_name.str()));
      node_parts.push_back(&meta.declare_part(node_part_name.str()));
    }
  }

  mesh.add_elem_parts(elem_parts.begin(), elem_parts.size());
  mesh.add_node_parts(node_parts.begin(), node_parts.size());

  meta.commit();

  mesh.generate_mesh();

  // Grab handles to the elements and nodes we'll be touching
  EntityVector elems_to_chg;
  EntityVector nodes_to_chg;
  elems_to_chg.reserve(x_dim);
  nodes_to_chg.reserve(x_dim);

  for (int x = 0; x < x_dim; ++x) {
    Entity elem = mesh.elem(x, 1 /*y*/, 1 /*z*/);
    elems_to_chg.push_back(elem);
    nodes_to_chg.push_back(bulk.begin_nodes(elem)[0]);
  }

  const int num_iterations = 20;

  bulk.modification_begin();

  CALLGRIND_START_INSTRUMENTATION;
  CALLGRIND_TOGGLE_COLLECT;
  // Remove and reattach the shared nodes from/to element1 over and over.
  // We don't want the node to change parts because we don't want to measure
  // entity movement costs. We just want to time the induced-part logic.
  for (int i = 0; i < num_iterations; ++i) {
    for (int x = 0; x < x_dim; ++x) {
      bulk.destroy_relation(elems_to_chg[x], nodes_to_chg[x], 0 /*rel id*/);
      bulk.declare_relation(elems_to_chg[x], nodes_to_chg[x], 0 /*rel id*/);
    }
  }
  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_STOP_INSTRUMENTATION;

  bulk.modification_end();

  print_memory_sum_all_procs(pm);

  print_debug_skip(pm);
}

}
