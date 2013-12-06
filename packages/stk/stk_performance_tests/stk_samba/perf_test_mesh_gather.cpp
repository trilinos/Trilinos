#if 0 // re-enable with 4.7.2

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <samba_fixtures/box_fixture.hpp>

#include <stk_performance_test_includes/calculate_centroid.hpp>
#include <stk_performance_test_includes/cpu_time.hpp>
#include <stk_performance_test_includes/memory_size.hpp>

#include <iostream>
#include <sys/resource.h>

namespace {

typedef samba::fixtures::box_fixture::const_coordinate_field_type coordinate_field;

template <typename Index>
struct gather_elem
{
  coordinate_field coordinates;
  double avg_centroid[3];

  gather_elem(coordinate_field arg_coordinates)
    : coordinates(arg_coordinates)
    , avg_centroid()
  {
    clear_centroid();
  }

  void clear_centroid()
  {
    avg_centroid[0] = 0;
    avg_centroid[1] = 0;
    avg_centroid[2] = 0;
  }

  void operator()( samba::partition_proxy partition)
  {
    //const int num_nodes = static_cast<int>(samba::num_nodes(partition.topology()));
    const int num_nodes = 8;

    double elem_centroid[3] = {0};
    std::vector<double> elem_node_coords(num_nodes*3,0);

    samba::partition_offset offset, end_offset;
    for(offset=0,end_offset=partition.size(); offset<end_offset; ++offset) {
      const Index * nodes = partition.begin_nodes<Index>(offset);

      for (int i=0; i<num_nodes; ++i, ++nodes) {
        Index node = *nodes;
        elem_node_coords[i*3+0] = coordinates[node][0];
        elem_node_coords[i*3+1] = coordinates[node][1];
        elem_node_coords[i*3+2] = coordinates[node][2];
      }
      stk::performance_tests::calculate_centroid_3d(num_nodes,&elem_node_coords[0],elem_centroid);

      avg_centroid[0] += elem_centroid[0]; elem_centroid[0] = 0;
      avg_centroid[1] += elem_centroid[1]; elem_centroid[1] = 0;
      avg_centroid[2] += elem_centroid[2]; elem_centroid[2] = 0;
    }
  }
};

template <typename Gather>
void run_test(samba::mesh mesh, coordinate_field coordinates, int iterations, const double expected, const double tolerance)
{
  const samba::set_expression hex_partitions = samba::entity_rank::element() & samba::entity_topology::hex_8();

  Gather gather(coordinates);

  const size_t num_elements = mesh.num_entities(samba::entity_rank::element());

  samba::partition_id partition;
  const samba::partition_id partition_end = samba::partition_id::create(mesh.num_partitions());

  for(int iteration=0; iteration<iterations; ++iteration)
  {

    for (partition=0; partition < partition_end; ++partition)
    {
      if ( hex_partitions.contains(mesh[partition]) ) {
        gather(mesh[partition]);
      }
    }

    for (int d=0; d<3; ++d) {
      EXPECT_LT(std::abs(gather.avg_centroid[d]/num_elements - expected), tolerance);
    }

    gather.clear_centroid();
  }
}

} // namespace

STKUNIT_UNIT_TEST(mesh_gather, mesh_gather)
{
// 01/30/12 tscoffe:  Note:  The Callgrind performance test gold standard is based on a 10x10x10 mesh with 100 iterations.
#if defined(STK_SAMBA_CALLGRIND)
  const uint32_t NX = 10;
  const uint32_t NY = 10;
  const uint32_t NZ = 10;
#elif defined(STK_SAMBA_LARGE_PERF_TEST)
  const uint32_t NX = 200;
  const uint32_t NY = 200;
  const uint32_t NZ = 200;
#else
  const uint32_t NX = 20;
  const uint32_t NY = 20;
  const uint32_t NZ = 20;
#endif

  samba::detail::report_memory_size(std::cout,"Initial");

  double start_time = cpu_time();

  samba::fixtures::box_fixture fixture(NX,NY,NZ);

  double end_time = cpu_time() - start_time;

  const double tolerance = 1.e-6;
  const double expected = NX/2.0;

  const int num_iterations = 100;

  samba::mesh mesh = fixture.mesh();
  coordinate_field coordinates = fixture.coordinates();

  std::cout << std::endl;
  std::cout << "Samba partition test:\n";
  std::cout << "\tNum Nodes: " << fixture.num_nodes() << std::endl;
  std::cout << "\tNum Elements: " << fixture.num_elements() << std::endl;

  std::cout << std::endl;
  std::cout << "Time to create mesh: " << end_time << std::endl;
  samba::detail::report_memory_size(std::cout,"Mesh Create");

  //compress mesh
  {
    start_time = cpu_time();

    mesh.end_modification();

    end_time = cpu_time() - start_time;

    std::cout << std::endl;
    std::cout << "Time to compress mesh: " << end_time << std::endl;
    samba::detail::report_memory_size(std::cout,"Mesh Compress");
  }

  // gather entity_key
  {
    start_time = cpu_time();
    run_test< gather_elem<samba::entity_key> >(mesh, coordinates, num_iterations, expected, tolerance);
    end_time = cpu_time() - start_time;

    std::cout << "\nTime to compute centroids (entity_key) ("<<num_iterations<<" iterations): " << end_time << std::endl;
  }

  // gather partition_index
  {
    start_time = cpu_time();
    run_test< gather_elem<samba::partition_index> >(mesh, coordinates, num_iterations, expected, tolerance);
    end_time = cpu_time() - start_time;

    std::cout << "\nTime to compute centroids (partition_index) ("<<num_iterations<<" iterations): " << end_time << std::endl;
  }

  // gather node_index
  {
    start_time = cpu_time();
    run_test< gather_elem<samba::node_index> >(mesh, coordinates, num_iterations, expected, tolerance);
    end_time = cpu_time() - start_time;

    std::cout << "\nTime to compute centroids (node_index) ("<<num_iterations<<" iterations): " << end_time << std::endl;
  }
}

#endif
