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
HexFixture* create_hex_with_complex_parts(stk::ParallelMachine pm, int x_dim, int y_dim, int z_dim, int chunk_span,
                                          std::vector<PartVector>& element_parts,
                                          std::vector<std::vector<std::vector<SimpleField*> > >* fields = NULL,
                                          int num_fields_per_chunk = 1)
{
  ThrowRequire(x_dim % chunk_span == 0);
  ThrowRequire(y_dim % chunk_span == 0);
  ThrowRequire(z_dim % chunk_span == 0);

  // Set up mesh
  HexFixture* mesh = new HexFixture(pm, x_dim, y_dim, z_dim);
  MetaData & meta = mesh->m_meta;
  BulkData & bulk = mesh->m_bulk_data;

  element_parts.clear();
  element_parts.resize(3);
  if (AllocateFields) {
    ThrowRequire(fields != NULL);
    ThrowRequire(num_fields_per_chunk > 0);
    fields->clear();
    fields->resize(3);
  }

  int dims[3] = {x_dim, y_dim, z_dim};
  char dim_names[3] = {'x', 'y', 'z'};
  const EntityRank field_entity_rank = Induce ? stk::topology::NODE_RANK : stk::topology::ELEMENT_RANK;

  for (int dim = 0; dim < 3; ++dim) {
    const int dim_size = dims[dim];
    const char dim_name = dim_names[dim];
    const int dim_chunks = dim_size / chunk_span;

    if (AllocateFields) {
      (*fields)[dim].resize(dim_chunks);
    }

    for (int i = 0; i < dim_chunks; ++i) {
      std::ostringstream oss;
      oss << dim_name << "_chunk_" << i;
      if (Induce) {
        element_parts[dim].push_back(&meta.declare_part(std::string("part_") + oss.str(), stk::topology::ELEMENT_RANK));
      }
      else {
        element_parts[dim].push_back(&meta.declare_part(std::string("part_") + oss.str()));
      }

      if (AllocateFields) {
        double init_val = 0.0;
        for (int f = 0; f < num_fields_per_chunk; ++f) {
          std::ostringstream append;
          append << f << "_";
          (*fields)[dim][i].push_back(&meta.declare_field<SimpleField>(static_cast<stk::topology::rank_t>(field_entity_rank), std::string("field_") + append.str() + oss.str()));
          put_field(*((*fields)[dim][i].back()), *element_parts[dim].back(), &init_val);
        }
      }
    }
  }

  ThrowRequire(element_parts[0].size() == static_cast<size_t>(x_dim / chunk_span));
  ThrowRequire(element_parts[1].size() == static_cast<size_t>(y_dim / chunk_span));
  ThrowRequire(element_parts[2].size() == static_cast<size_t>(z_dim / chunk_span));

  if (AllocateFields) {
    ThrowRequire((*fields)[0].size() == static_cast<size_t>(x_dim / chunk_span));
    ThrowRequire((*fields)[1].size() == static_cast<size_t>(y_dim / chunk_span));
    ThrowRequire((*fields)[2].size() == static_cast<size_t>(z_dim / chunk_span));
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
        add[0] = element_parts[0][x / chunk_span];
        add[1] = element_parts[1][y / chunk_span];
        add[2] = element_parts[2][z / chunk_span];

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
  const int chunk_span = 1; // max fragmentation

  std::vector<PartVector> element_parts; // element_parts[dim][chunk]
  HexFixture* mesh =
    create_hex_with_complex_parts< !TIME_CHANGE_PARTS, !INDUCE_ELEMENT_PARTS, !ALLOCATE_FIELDS >(pm, x_dim, y_dim, z_dim, chunk_span, element_parts);
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

          bulk.get_buckets(stk::topology::ELEMENT_RANK, sel, output_buckets);
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
  const int chunk_span = 1; // max fragmentation

  std::vector<PartVector> element_parts; // element_parts[dim][chunk]
  HexFixture* mesh =
    create_hex_with_complex_parts< !TIME_CHANGE_PARTS, !INDUCE_ELEMENT_PARTS, !ALLOCATE_FIELDS >(pm, x_dim, y_dim, z_dim, chunk_span, element_parts);
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
  const int chunk_span = 1; // max fragmentation

  std::vector<PartVector> element_parts; // element_parts[dim][chunk]
  HexFixture* mesh =
    create_hex_with_complex_parts< !TIME_CHANGE_PARTS, !INDUCE_ELEMENT_PARTS, !ALLOCATE_FIELDS >(pm, x_dim, y_dim, z_dim, chunk_span, element_parts);
  BulkData & bulk = mesh->m_bulk_data;

  const int num_iterations = 5;

  bulk.modification_begin();

  CALLGRIND_TOGGLE_COLLECT;

  for (int i = 0; i < num_iterations; ++i) {
    for (int x = 0; x < x_dim; ++x) {
      for (int y = 0; y < y_dim; ++y) {
        int size_sum = 0;
        for (int z = 0; z < z_dim; ++z) {
          Selector sel = *element_parts[2][z];

          BucketVector const& output_buckets = bulk.get_buckets(stk::topology::ELEMENT_RANK, sel);
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
  const int chunk_span = 1; // max fragmentation

  // Set up mesh
  // This mesh uses incredible amounts of memory with part induction turned on
  std::vector<PartVector> element_parts; // element_parts[dim][chunk]
  HexFixture* mesh =
    create_hex_with_complex_parts< TIME_CHANGE_PARTS, !INDUCE_ELEMENT_PARTS, !ALLOCATE_FIELDS >(pm, x_dim, y_dim, z_dim, chunk_span, element_parts);

  PERFORMANCE_TEST_POSTAMBLE();

  delete mesh;
}

void fill_buckets_for_field_tests(BulkData const& bulk,
                                  std::vector<PartVector> const& element_parts,
                                  std::vector<std::vector<std::vector<BucketVector> > >& bucket_map,
                                  int x_chunks, int y_chunks, int z_chunks, int chunk_span)
{
  bucket_map.clear();

  bucket_map.resize(x_chunks);

  const size_t nodes_per_chunk = (chunk_span + 1) * (chunk_span + 1) * (chunk_span + 1);

  for (int x = 0; x < x_chunks; ++x) {
    bucket_map[x].resize(y_chunks);
    for (int y = 0; y < y_chunks; ++y) {
      bucket_map[x][y].resize(z_chunks);
      for (int z = 0; z < z_chunks; ++z) {
        Selector sel = *element_parts[0][x];
        sel &= *element_parts[1][y];
        sel &= *element_parts[2][z];

        BucketVector& chunk_buckets = bucket_map[x][y][z];
        bulk.get_buckets(stk::topology::NODE_RANK, sel, chunk_buckets);

        size_t chunk_size = 0;
        for (int b = 0, be = chunk_buckets.size(); b < be; ++b) {
          chunk_size += chunk_buckets[b]->size();
        }
        STKUNIT_EXPECT_EQ(nodes_per_chunk, chunk_size);
      }
    }
  }
}

STKUNIT_UNIT_TEST( stk_mesh_perf_unit_test, field_access)
{
  PERFORMANCE_TEST_PREAMBLE(1 /*num procs*/);

  const int x_dim = 20;
  const int y_dim = 20;
  const int z_dim = 20;
  const int chunk_span = 5; // low fragmentation
  const int num_fields_per_chunk = 3;

  // Set up mesh
  std::vector<PartVector> element_parts; // element_parts[dim][chunk]
  std::vector<std::vector<std::vector<SimpleField*> > > fields; // fields[dim][chunk][field_offset]
  HexFixture* mesh =
    create_hex_with_complex_parts< !TIME_CHANGE_PARTS, INDUCE_ELEMENT_PARTS, ALLOCATE_FIELDS>(pm, x_dim, y_dim, z_dim, chunk_span, element_parts, &fields, num_fields_per_chunk);
  BulkData & bulk = mesh->m_bulk_data;

  const int x_chunks = x_dim / chunk_span;
  const int y_chunks = y_dim / chunk_span;
  const int z_chunks = z_dim / chunk_span;

  const size_t nodes_per_chunk = (chunk_span + 1) * (chunk_span + 1) * (chunk_span + 1);

  std::vector<std::vector<std::vector<BucketVector> > > bucket_map; // bucket_map[x_chunk][y_chunk][z_chunk]
  fill_buckets_for_field_tests(bulk, element_parts, bucket_map, x_chunks, y_chunks, z_chunks, chunk_span);

  CALLGRIND_TOGGLE_COLLECT;

  const int num_iterations = 100;
  size_t dummy = 0;
  for (int i = 0; i < num_iterations; ++i) {
    for (int x = 0; x < x_chunks; ++x) {
      std::vector<SimpleField*> const& x_fields = fields[0][x];
      for (int y = 0; y < y_chunks; ++y) {
        std::vector<SimpleField*> const& y_fields = fields[1][y];
        for (int z = 0; z < z_chunks; ++z) {
          std::vector<SimpleField*> const& z_fields = fields[2][z];
          BucketVector const& chunk_buckets = bucket_map[x][y][z];
          for (int b = 0, be = chunk_buckets.size(); b < be; ++b) {
            Bucket& bucket = *bucket_map[x][y][z][b];

            // We are intentionally not using the much-faster style of accessing
            // fields by bucket instead of by entity.

            for (size_t e = 0, ee = bucket.size(); e < ee; ++e) {
              Entity node = bucket[e];

              stk::mesh::MeshIndex mi = bulk.mesh_index(node);
              for (int f = 0; f < num_fields_per_chunk; ++f) {
                const double* x_field_data = bulk.field_data(*x_fields[f], *mi.bucket, mi.bucket_ordinal);
                const double* y_field_data = bulk.field_data(*y_fields[f], *mi.bucket, mi.bucket_ordinal);
                const double* z_field_data = bulk.field_data(*z_fields[f], *mi.bucket, mi.bucket_ordinal);

                if (x_field_data != NULL && y_field_data != NULL && z_field_data != NULL) {
                  ++dummy;
		}
              }
            }
          }
        }
      }
    }
  }

  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_STOP_INSTRUMENTATION;

  STKUNIT_EXPECT_EQ(dummy, nodes_per_chunk * x_chunks * y_chunks * z_chunks * num_iterations * num_fields_per_chunk);

  PERFORMANCE_TEST_POSTAMBLE();

  delete mesh;
}

STKUNIT_UNIT_TEST( stk_mesh_perf_unit_test, field_access_sm_style)
{
  PERFORMANCE_TEST_PREAMBLE(1 /*num procs*/);

  const int x_dim = 20;
  const int y_dim = 20;
  const int z_dim = 20;
  const int chunk_span = 5; // low fragmentation
  const int num_fields_per_chunk = 3;

  // Set up mesh
  std::vector<PartVector> element_parts; // element_parts[dim][chunk]
  std::vector<std::vector<std::vector<SimpleField*> > > fields; // fields[dim][chunk][field_offset_within_chunk]
  HexFixture* mesh =
    create_hex_with_complex_parts< !TIME_CHANGE_PARTS, INDUCE_ELEMENT_PARTS, ALLOCATE_FIELDS>(pm, x_dim, y_dim, z_dim, chunk_span, element_parts, &fields, num_fields_per_chunk);
  BulkData & bulk = mesh->m_bulk_data;

  const int x_chunks = x_dim / chunk_span;
  const int y_chunks = y_dim / chunk_span;
  const int z_chunks = z_dim / chunk_span;

  const size_t nodes_per_chunk = (chunk_span + 1) * (chunk_span + 1) * (chunk_span + 1);

  std::vector<std::vector<std::vector<BucketVector> > > bucket_map; // bucket_map[x_chunk][y_chunk][z_chunk]
  fill_buckets_for_field_tests(bulk, element_parts, bucket_map, x_chunks, y_chunks, z_chunks, chunk_span);



  CALLGRIND_TOGGLE_COLLECT;

  //
  // Model caching of the bucket index and ordinals
  //
  std::vector<std::vector<unsigned int  > >bucket_ind; 
  std::vector<std::vector<unsigned short> >bucket_ord; 
  bucket_ind.resize(345);
  bucket_ord.resize(345);
  for (int x = 0; x < x_chunks; ++x) {
    for (int y = 0; y < y_chunks; ++y) {
      for (int z = 0; z < z_chunks; ++z) {
        BucketVector const& chunk_buckets = bucket_map[x][y][z];
        for (int b = 0, be = chunk_buckets.size(); b < be; ++b) {
          Bucket& bucket = *bucket_map[x][y][z][b];
          int bucket_id = bucket.bucket_id();
          if(bucket_ind[bucket_id].size() > 0) continue;
          for (size_t e = 0, ee = bucket.size(); e < ee; ++e) {
            bucket_ind[bucket_id].push_back(bucket_id);
            bucket_ord[bucket_id].push_back(e);
	  }
	}
      }
    }
  }
  
  FieldMetaDataVector const* default_val = NULL;
  std::vector<FieldMetaDataVector const*> x_field_meta(num_fields_per_chunk, default_val);
  std::vector<FieldMetaDataVector const*> y_field_meta(num_fields_per_chunk, default_val);
  std::vector<FieldMetaDataVector const*> z_field_meta(num_fields_per_chunk, default_val);
  
  const int num_iterations = 100;
  size_t dummy = 0;
   for (int i = 0; i < num_iterations; ++i) {
    for (int x = 0; x < x_chunks; ++x) {
      std::vector<SimpleField*> const& x_fields = fields[0][x];
      for (int y = 0; y < y_chunks; ++y) {
        std::vector<SimpleField*> const& y_fields = fields[1][y];
        for (int z = 0; z < z_chunks; ++z) {
          std::vector<SimpleField*> const& z_fields = fields[2][z];

          for (int f = 0; f < num_fields_per_chunk; ++f) {
            x_field_meta[f] =  &x_fields[f]->get_meta_data_for_field();
            y_field_meta[f] =  &y_fields[f]->get_meta_data_for_field();
            z_field_meta[f] =  &z_fields[f]->get_meta_data_for_field();
          }

          BucketVector const& chunk_buckets = bucket_map[x][y][z];
          for (int b = 0, be = chunk_buckets.size(); b < be; ++b) {
            Bucket& bucket = *bucket_map[x][y][z][b];
            int bid = bucket.bucket_id();

            std::vector<unsigned int  >& bucket_indZ = bucket_ind[bid]; 
            std::vector<unsigned short>& bucket_ordZ = bucket_ord[bid]; 


            for (size_t e = 0, ee = bucket_indZ.size(); e < ee; ++e) {
              unsigned int bucket_id = bucket_indZ[e];
              unsigned short bucket_ord = bucket_ordZ[e];

              // We are intentionally not using the much-faster style of accessing
              // fields by bucket instead of by entity.

	      for (int f = 0; f < num_fields_per_chunk; ++f) {

                const double* x_field_data = reinterpret_cast<const double*>((*x_field_meta[f])[bucket_id].m_data) + bucket_ord;
                const double* y_field_data = reinterpret_cast<const double*>((*y_field_meta[f])[bucket_id].m_data) + bucket_ord;
                const double* z_field_data = reinterpret_cast<const double*>((*z_field_meta[f])[bucket_id].m_data) + bucket_ord;


                if (x_field_data != NULL && y_field_data != NULL && z_field_data != NULL) {
                  ++dummy;
                }
              }
            }
          }
        }
      }
    }
  }

  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_STOP_INSTRUMENTATION;

  STKUNIT_EXPECT_EQ(dummy, nodes_per_chunk * x_chunks * y_chunks * z_chunks * num_iterations * num_fields_per_chunk);

  PERFORMANCE_TEST_POSTAMBLE();

  delete mesh;
}

}
