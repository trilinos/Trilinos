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
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/BoundaryAnalysis.hpp>

#include <sstream>
#include <valgrind/callgrind.h>

namespace {

using namespace stk::mesh;
using stk::mesh::fixtures::HexFixture;


#define PERFORMANCE_TEST_PREAMBLE(expected_np)          \
  CALLGRIND_START_INSTRUMENTATION;                      \
                                                        \
  check_valgrind_version();                             \
                                                        \
  stk::ParallelMachine pm = MPI_COMM_WORLD;             \
                                                        \
  const int p_size = stk::parallel_machine_size(pm);    \
                                                        \
  ThrowRequire(p_size == expected_np)


#define PERFORMANCE_TEST_POSTAMBLE()            \
  print_memory_sum_all_procs(pm);               \
                                                \
  print_debug_skip(pm)


void print_memory_sum_all_procs(stk::ParallelMachine pm)
{
  const size_t p_rank = stk::parallel_machine_rank(pm);
  size_t my_peak = stk::allocator_memory_usage<void>::peak_memory();
  size_t peak_sum = 0;
  int err = MPI_Reduce(&my_peak, &peak_sum, 1 /*size*/, MPI_LONG_LONG, MPI_SUM, 0 /*root*/, pm);
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
  PERFORMANCE_TEST_PREAMBLE(1 /*num procs*/);

  const int x_dim = 500;
  const int y_dim = 3;
  const int z_dim = 3;

  HexFixture mesh(pm, x_dim, y_dim, z_dim);
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

  const int num_iterations = 40;

  bulk.modification_begin();

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

  PERFORMANCE_TEST_POSTAMBLE();
}

void mesh_create_hex_with_edges_test(stk::ParallelMachine pm)
{
  const int p_size = stk::parallel_machine_size(pm);

  const int x_dim = 10 * p_size;
  const int y_dim = 10;
  const int z_dim = 10;

  // Set up mesh
  HexFixture mesh(pm, x_dim, y_dim, z_dim);
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

  CALLGRIND_TOGGLE_COLLECT;

  meta.commit();

  mesh.generate_mesh();

  bulk.modification_begin();

  stk::mesh::create_edges(bulk);

  bulk.modification_end();

  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_STOP_INSTRUMENTATION;
}

STKUNIT_UNIT_TEST( stk_mesh_perf_unit_test, mesh_create_hex_with_edges_serial )
{
  PERFORMANCE_TEST_PREAMBLE(1 /*num procs*/);

  mesh_create_hex_with_edges_test(pm);

  PERFORMANCE_TEST_POSTAMBLE();
}

STKUNIT_UNIT_TEST( stk_mesh_perf_unit_test, mesh_create_hex_with_edges_parallel )
{
  PERFORMANCE_TEST_PREAMBLE(8 /*num procs*/);

  mesh_create_hex_with_edges_test(pm);

  PERFORMANCE_TEST_POSTAMBLE();
}

// caller will need to delete returned ptr
const bool TIME_CHANGE_PARTS = true;
const bool INDUCE_ELEMENT_PARTS = true;
const bool ALLOCATE_FIELDS = true;

typedef Field<double> SimpleField;

template <bool TimeChangeParts, bool Induce, bool AllocateFields>
HexFixture* create_hex_with_complex_parts(stk::ParallelMachine pm, int x_dim, int y_dim, int z_dim, int dim_span,
                                          std::vector<PartVector>& element_parts,
                                          std::vector<std::vector<SimpleField*> >* fields = NULL)
{
  ThrowRequire(x_dim % dim_span == 0);
  ThrowRequire(y_dim % dim_span == 0);
  ThrowRequire(z_dim % dim_span == 0);

  // Set up mesh
  HexFixture* mesh = new HexFixture(pm, x_dim, y_dim, z_dim);
  MetaData & meta = mesh->m_meta;
  BulkData & bulk = mesh->m_bulk_data;

  element_parts.clear();
  element_parts.resize(3);
  if (AllocateFields) {
    fields->clear();
    fields->resize(3);
  }

  int dims[3] = {x_dim, y_dim, z_dim};
  char dim_names[3] = {'x', 'y', 'z'};

  for (int dim = 0; dim < 3; ++dim) {
    int dim_size = dims[dim];
    char dim_name = dim_names[dim];
    for (int i = 0; i < dim_size; ++i) {
      if (i % dim_span == 0) {
        std::ostringstream oss;
        oss << dim_name << "_" << i / dim_span;
        if (Induce) {
          element_parts[dim].push_back(&meta.declare_part(std::string("part_") + oss.str(), stk::topology::ELEMENT_RANK));
        }
        else {
          element_parts[dim].push_back(&meta.declare_part(std::string("part_") + oss.str()));
        }

        if (AllocateFields) {
          double init_val = 0.0;
          (*fields)[dim].push_back(&meta.declare_field<SimpleField>(std::string("field_") + oss.str()));
          put_field(*((*fields)[dim].back()), stk::topology::ELEMENT_RANK, *element_parts[dim].back(), &init_val);
        }
      }
    }
  }

  ThrowRequire(element_parts[0].size() == static_cast<size_t>(x_dim));
  ThrowRequire(element_parts[1].size() == static_cast<size_t>(y_dim));
  ThrowRequire(element_parts[2].size() == static_cast<size_t>(z_dim));

  if (AllocateFields) {
    ThrowRequire((*fields)[0].size() == static_cast<size_t>(x_dim));
    ThrowRequire((*fields)[1].size() == static_cast<size_t>(y_dim));
    ThrowRequire((*fields)[2].size() == static_cast<size_t>(z_dim));
  }

  meta.commit();

  mesh->generate_mesh();

  bulk.modification_begin();

  if (TimeChangeParts) {
    CALLGRIND_TOGGLE_COLLECT;
  }

  for (int x = 0; x < x_dim; ++x) {
    for (int y = 0; y < y_dim; ++y) {
      for (int z = 0; z < z_dim; ++z) {
        PartVector add(3);
        add[0] = element_parts[0][x / dim_span];
        add[1] = element_parts[1][y / dim_span];
        add[2] = element_parts[2][z / dim_span];

        Entity elem = mesh->elem(x, y, z);
        bulk.change_entity_parts(elem, add);
      }
    }
  }

  if (TimeChangeParts) {
    CALLGRIND_TOGGLE_COLLECT;
    CALLGRIND_STOP_INSTRUMENTATION;
  }

  bulk.modification_end();

  return mesh;
}

STKUNIT_UNIT_TEST( stk_mesh_perf_unit_test, frag_mesh_selector )
{
  PERFORMANCE_TEST_PREAMBLE(1 /*num procs*/);

  const int x_dim = 10;
  const int y_dim = 10;
  const int z_dim = 10;
  const int dim_span = 1; // max fragmentation

  std::vector<PartVector> element_parts;
  HexFixture* mesh =
    create_hex_with_complex_parts< !TIME_CHANGE_PARTS, !INDUCE_ELEMENT_PARTS, !ALLOCATE_FIELDS >(pm, x_dim, y_dim, z_dim, dim_span, element_parts);
  BulkData & bulk = mesh->m_bulk_data;

  const int num_iterations = 5;

  bulk.modification_begin();

  CALLGRIND_TOGGLE_COLLECT;

  BucketVector output_buckets;

  for (int i = 0; i < num_iterations; ++i) {
    for (int x = 0; x < x_dim; ++x) {
      for (int y = 0; y < y_dim; ++y) {
        int size_sum = 0;
        for (int z = 0; z < z_dim; ++z) {
          Selector sel = *element_parts[0][x];
          sel &= *element_parts[1][y];
          sel &= *element_parts[2][z];

          get_buckets(sel, bulk.buckets(stk::topology::ELEMENT_RANK), output_buckets);
          for (int j = 0, je = output_buckets.size(); j < je; ++j) {
            size_sum += output_buckets[j]->size();
          }
        }
        STKUNIT_EXPECT_EQ(size_sum, z_dim);
      }
    }
  }

  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_STOP_INSTRUMENTATION;

  bulk.modification_end();

  PERFORMANCE_TEST_POSTAMBLE();

  delete mesh;
}

STKUNIT_UNIT_TEST( stk_mesh_perf_unit_test, frag_mesh_selector_fast )
{
  PERFORMANCE_TEST_PREAMBLE(1 /*num procs*/);

  const int x_dim = 10;
  const int y_dim = 10;
  const int z_dim = 10;
  const int dim_span = 1; // max fragmentation

  std::vector<PartVector> element_parts;
  HexFixture* mesh =
    create_hex_with_complex_parts< !TIME_CHANGE_PARTS, !INDUCE_ELEMENT_PARTS, !ALLOCATE_FIELDS >(pm, x_dim, y_dim, z_dim, dim_span, element_parts);
  BulkData & bulk = mesh->m_bulk_data;

  const int num_iterations = 5;

  CALLGRIND_TOGGLE_COLLECT;

  for (int i = 0; i < num_iterations; ++i) {
    for (int x = 0; x < x_dim; ++x) {
      for (int y = 0; y < y_dim; ++y) {
        int size_sum = 0;
        for (int z = 0; z < z_dim; ++z) {
          Selector sel = *element_parts[0][x];
          sel &= *element_parts[1][y];
          sel &= *element_parts[2][z];

          BucketVector const& output_buckets = bulk.get_buckets(stk::topology::ELEMENT_RANK, sel);
          for (int j = 0, je = output_buckets.size(); j < je; ++j) {
            size_sum += output_buckets[j]->size();
          }
        }
        STKUNIT_EXPECT_EQ(size_sum, z_dim);
      }
    }
  }

  STKUNIT_EXPECT_EQ(1000u, bulk.buckets(stk::topology::ELEMENT_RANK).size());
  STKUNIT_EXPECT_EQ(3u,    bulk.buckets(stk::topology::NODE_RANK).size());

  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_STOP_INSTRUMENTATION;

  PERFORMANCE_TEST_POSTAMBLE();

  delete mesh;
}

STKUNIT_UNIT_TEST( stk_mesh_perf_unit_test, frag_mesh_single_part_selector )
{
  PERFORMANCE_TEST_PREAMBLE(1 /*num procs*/);

  const int x_dim = 10;
  const int y_dim = 10;
  const int z_dim = 10;
  const int dim_span = 1; // max fragmentation

  std::vector<PartVector> element_parts;
  HexFixture* mesh =
    create_hex_with_complex_parts< !TIME_CHANGE_PARTS, !INDUCE_ELEMENT_PARTS, !ALLOCATE_FIELDS >(pm, x_dim, y_dim, z_dim, dim_span, element_parts);
  BulkData & bulk = mesh->m_bulk_data;

  const int num_iterations = 5;

  bulk.modification_begin();

  CALLGRIND_TOGGLE_COLLECT;

  BucketVector output_buckets;

  for (int i = 0; i < num_iterations; ++i) {
    for (int x = 0; x < x_dim; ++x) {
      for (int y = 0; y < y_dim; ++y) {
        int size_sum = 0;
        for (int z = 0; z < z_dim; ++z) {
          Selector sel = *element_parts[2][z];

          get_buckets(sel, bulk.buckets(stk::topology::ELEMENT_RANK), output_buckets);
          for (int j = 0, je = output_buckets.size(); j < je; ++j) {
            size_sum += output_buckets[j]->size();
          }
        }
        STKUNIT_EXPECT_EQ(size_sum, x_dim * y_dim * z_dim);
      }
    }
  }

  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_STOP_INSTRUMENTATION;

  bulk.modification_end();

  PERFORMANCE_TEST_POSTAMBLE();

  delete mesh;
}

STKUNIT_UNIT_TEST( stk_mesh_perf_unit_test, frag_mesh_memory)
{
  PERFORMANCE_TEST_PREAMBLE(1 /*num procs*/);

  const int x_dim = 10;
  const int y_dim = 10;
  const int z_dim = 10;
  const int dim_span = 1; // max fragmentation

  // Set up mesh
  // This mesh uses incredible amounts of memory with part induction turned on
  std::vector<PartVector> element_parts;
  HexFixture* mesh =
    create_hex_with_complex_parts< TIME_CHANGE_PARTS, !INDUCE_ELEMENT_PARTS, !ALLOCATE_FIELDS >(pm, x_dim, y_dim, z_dim, dim_span, element_parts);

  PERFORMANCE_TEST_POSTAMBLE();

  delete mesh;
}

STKUNIT_UNIT_TEST( stk_mesh_perf_unit_test, field_access)
{
  PERFORMANCE_TEST_PREAMBLE(1 /*num procs*/);

  const int x_dim = 10;
  const int y_dim = 10;
  const int z_dim = 10;
  const int dim_span = 1; // max fragmentation

  // Set up mesh
  // This mesh uses incredible amounts of memory with part induction turned on
  std::vector<PartVector> element_parts;
  std::vector<std::vector<SimpleField*> > fields;
  HexFixture* mesh =
    create_hex_with_complex_parts< !TIME_CHANGE_PARTS, !INDUCE_ELEMENT_PARTS, ALLOCATE_FIELDS>(pm, x_dim, y_dim, z_dim, dim_span, element_parts, &fields);
  BulkData & bulk = mesh->m_bulk_data;

  BucketVector all_elem_buckets = bulk.buckets(stk::topology::ELEMENT_RANK);

  std::vector<std::vector<BucketVector> > bucket_map;
  bucket_map.resize(x_dim);

  for (int x = 0; x < x_dim; ++x) {
    bucket_map[x].resize(y_dim);
    for (int y = 0; y < y_dim; ++y) {
      bucket_map[x][y].resize(z_dim);
      for (int z = 0; z < z_dim; ++z) {
        Selector sel = *element_parts[0][x];
        sel &= *element_parts[1][y];
        sel &= *element_parts[2][z];

        BucketVector selected_buckets;
        get_buckets(sel, all_elem_buckets, selected_buckets);
        ThrowRequire(selected_buckets.size() == 1u);
        bucket_map[x][y][z] = selected_buckets[0];
      }
    }
  }

  CALLGRIND_TOGGLE_COLLECT;

  const int num_iterations = 5000;
  int dummy = 0;
  for (int i = 0; i < num_iterations; ++i) {
    for (int x = 0; x < x_dim; ++x) {
      SimpleField& x_field = *fields[0][x];
      for (int y = 0; y < y_dim; ++y) {
        SimpleField& y_field = *fields[1][y];
        for (int z = 0; z < z_dim; ++z) {
          SimpleField& z_field = *fields[2][z];
          Bucket& bucket = *bucket_map[x][y][z];

          const double* x_field_data = bulk.field_data(x_field, bucket);
          const double* y_field_data = bulk.field_data(y_field, bucket);
          const double* z_field_data = bulk.field_data(z_field, bucket);

          if (x_field_data != NULL && y_field_data != NULL && z_field_data != NULL) {
            ++dummy;
          }
        }
      }
    }
  }

  STKUNIT_EXPECT_EQ(dummy, x_dim * y_dim * z_dim * num_iterations);

  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_STOP_INSTRUMENTATION;


  PERFORMANCE_TEST_POSTAMBLE();

  delete mesh;
}

}
