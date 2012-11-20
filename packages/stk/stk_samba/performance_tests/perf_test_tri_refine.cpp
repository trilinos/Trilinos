#include <gtest/gtest.h>

#include <samba_fixtures/tri_fixture.hpp>
#include <samba/mesh.hpp>

#include <iostream>
#include <iterator>
#include <algorithm>
#include <sstream>
#include <vector>
#include <cstdio>

#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/equal.hpp>

#include <performance_tests/cpu_time.hpp>
#include <performance_tests/memory_size.hpp>

#include <sys/resource.h>

namespace {

// 01/30/12 tscoffe:  Note:  The Callgrind performance test gold standard is based on a 10x10x10 mesh with 100 iterations.
#if defined(STK_SAMBA_CALLGRIND)
  static const uint32_t NX = 10;
  static const uint32_t NY = 10;
  //static const uint32_t NZ = 10;
#elif defined(STK_SAMBA_LARGE_PERF_TEST)
  static const uint32_t NX = 2000;
  static const uint32_t NY = 2000;
  //static const uint32_t NZ = 200;
#else
  static const uint32_t NX = 200;
  static const uint32_t NY = 200;
  //static const uint32_t NZ = 20;
#endif

#define DO_MEMORY_ACCOUNTING 1

typedef size_t MemorySizeType;

struct MemoryInfo
{
  MemorySizeType m_rss_current;
  MemorySizeType m_rss_high_water_mark;
  static const MemorySizeType MB = 1024*1024;

  MemoryInfo() { get_memory_usage(); }

  void get_memory_usage() 
  {
    m_rss_current = 0;
    m_rss_high_water_mark = 0;
    stk::get_memory_usage(m_rss_current, m_rss_high_water_mark);
  }
  void set_state() { get_memory_usage(); }
  void get_increment() {
    MemoryInfo old_state = *this;
    get_memory_usage();
    m_rss_current -= old_state.m_rss_current;
  }

};

inline double MegaByte(MemorySizeType x) { return  ((double)x/1024.0/1024.0); }

std::ostream& operator<<(std::ostream& os, const MemoryInfo& mem)
{
  char buf[1024];
  sprintf(buf, "\n%20s %20s\n%20g %20g\n", "current_rss [MB]", "max_rss [MB]",
          MegaByte(mem.m_rss_current), MegaByte(mem.m_rss_high_water_mark) );
  os << buf;
  return os;
}

/// first cut at a refinement algorithm - just adds 4 new tris and 12 new nodes for each original tri,
/// as a very simple model of what happens in a real refinement algorithm
/// Also, we don't yet deal with coordinates, fields, etc. that need to be interpolated at new nodes;
/// and, we don't deal with trying to ensure new nodes are shared on edges, etc.
static void refine_partition( samba::partition_proxy partition, samba::mesh mesh)
{
  const size_t num_nodes = samba::num_nodes(partition.topology());
  //std::cout << "num_nodes= " << num_nodes << std::endl;
  if (num_nodes != 3) throw std::runtime_error("refine_partition 1");

  std::vector<samba::entity_key> to_remove;

  samba::partition_offset offset, end_offset;
  for(offset=0,end_offset=partition.size(); offset<end_offset; ++offset) {
    //samba::entity_connectivity_iterator nodes = partition.begin_extended_nodes(offset);

    // for (size_t i=0; i<num_nodes; ++i, ++nodes) {
    //   samba::entity_key node = nodes->key();
    //   (void)node;
    // }

    // FIXME, for now we just create new nodes at midpoints and don't try to re-use, giving a nonconformal mesh
    const size_t num_new_tris = 4;
    samba::entity_key_interval new_nodes = mesh.add_entities( samba::entity_topology::node(), num_new_tris*num_nodes );
    samba::entity_key_interval new_tris = mesh.add_entities( samba::entity_topology::triangle_3(), num_new_tris );
    size_t i_node=0;
    for (size_t i=0; i<num_new_tris; ++i) {
      samba::connectivity_ordinal ordinal = {0};
      for (size_t j=0; j < num_nodes; ++j) {
        mesh.add_connectivity(new_tris[i],new_nodes[i_node],ordinal);
        mesh.add_connectivity(new_nodes[i_node],new_tris[i],ordinal);
        ++ordinal;
        ++i_node;
      }
    }

    samba::entity_key  k = partition[offset].key();
    to_remove.push_back(k);
  }

  mesh.remove_entities(to_remove.begin(), to_remove.end());
}

static void refine_tri_mesh(samba::mesh mesh)
{
  mesh.begin_modification();

  samba::partition_id partition,partition_end;
  for (partition=0,partition_end=mesh.num_partitions(); partition < partition_end; ++partition) {
    // ??
    if (samba::contains(samba::entity_topology::triangle_3() - samba::entity_rank::node(), mesh[partition])) {
      refine_partition(mesh[partition], mesh);
    }
  }

  mesh.end_modification();
}

static void count_memory(samba::mesh mesh)
{
  double time = -cpu_time();
  // a more refined test of memory per node, element, with some connectivity
  const size_t num_new_tris = NX*NY;
  const size_t num_nodes_per_tri = 3;
  const size_t num_new_nodes = num_new_tris*num_nodes_per_tri;
  MemoryInfo mem_delta_node;

  //samba::detail::report_memory_size(std::cout,"Exact: before create nodes");

  mesh.begin_modification();
  samba::entity_key_interval new_nodes = mesh.add_entities( samba::entity_topology::node(), num_new_nodes );
  mesh.end_modification();

  //samba::detail::report_memory_size(std::cout,"Exact: after create nodes");

  mem_delta_node.get_increment();
  double mem_per_node = double(mem_delta_node.m_rss_current)/double(num_new_nodes);
  std::cout << "\nsamba count_memory mem_per_node = " << mem_per_node << "\n" << std::endl;

  MemoryInfo mem_delta_elem_0, mem_delta_elem_1;

  //samba::detail::report_memory_size(std::cout,"Exact: before create elems");

  mesh.begin_modification();
  samba::entity_key_interval new_tris = mesh.add_entities( samba::entity_topology::triangle_3(), num_new_tris );
  mesh.end_modification();

  mem_delta_elem_0.get_increment();
  //samba::detail::report_memory_size(std::cout,"Exact: after create elems");

  //samba::detail::report_memory_size(std::cout,"Exact: before add conn");

  mesh.begin_modification();
  size_t i_node=0;
  for (size_t i=0; i<num_new_tris; ++i) {
    samba::connectivity_ordinal ordinal = {0};
    for (size_t j=0; j < num_nodes_per_tri; ++j) {
      mesh.add_connectivity(new_tris[i],new_nodes[i_node],ordinal);
      mesh.add_connectivity(new_nodes[i_node],new_tris[i],ordinal);
      ++ordinal;
      ++i_node;
    }
  }
  mesh.end_modification();
  //samba::detail::report_memory_size(std::cout,"Exact: after add conn");

  mem_delta_elem_1.get_increment();
  double mem_per_elem_0 = double(mem_delta_elem_0.m_rss_current)/double(num_new_tris);
  double mem_per_elem_1 = double(mem_delta_elem_1.m_rss_current)/double(num_new_tris);
  time += cpu_time();
  std::cout << "\nsamba count_memory mem_per_elem (no connectivity) = " << mem_per_elem_0 << " with connectivity= " << mem_per_elem_1 << " cpu = " << time << std::endl;


}

} // namespace

TEST(samba, tri_refine)
{

  samba::mesh mesh( samba::connectivity_map::default_map_2d());

  double mesh_create_time = -cpu_time();

  MemoryInfo mem_delta, mem_absolute;
  std::cout << "Memory: initial= " << mem_absolute << std::endl;

  samba::fixtures::tri_fixture tf(NX,NY);
  tf.samba_mesh_create(mesh);

  mesh_create_time += cpu_time();

  double mesh_refine_time = -cpu_time();

  unsigned num_nodes_0 = mesh.num_nodes();
  unsigned num_elems_0 = mesh.num_elements();
  mem_delta.get_increment();
  mem_absolute.set_state();
  double mem_per_node = double(mem_delta.m_rss_current)/double(mesh.num_nodes());
  double mem_per_element = double(mem_delta.m_rss_current)/double(mesh.num_elements());
  std::cout << "Create mesh: Time= " << mesh_refine_time << " Memory delta= " << mem_delta << "\n   absolute= " << mem_absolute 
            << "\n approx mem_per_node = " << mem_per_node << " [B] approx mem_per_element= " << mem_per_element << " [B]"
            << "\n" << std::endl;

  std::cout << "orig mesh.num_nodes()=  " << mesh.num_nodes() << std::endl;
  std::cout << "orig mesh.num_elements()=  " << mesh.num_elements() << std::endl;
  EXPECT_TRUE(mesh.num_nodes() == (NX+1)*(NY+1));
  EXPECT_TRUE(mesh.num_elements() == NX*NY*2);
  mem_delta.set_state();

  refine_tri_mesh(mesh);

  mesh_refine_time += cpu_time();
  unsigned num_nodes_delta = mesh.num_nodes()-num_nodes_0;
  unsigned num_elems_delta = mesh.num_elements() - num_elems_0;
  std::cout << " num_nodes_delta= " << num_nodes_delta << " num_elems_delta= " << num_elems_delta << std::endl;

  mem_delta.get_increment();
  mem_absolute.set_state();
  mem_per_node = double(mem_delta.m_rss_current)/double(num_nodes_delta);
  mem_per_element = double(mem_delta.m_rss_current)/double(num_elems_delta);
  std::cout << "Refine mesh: Time= " << mesh_refine_time << " Memory delta= " << mem_delta << "\n   absolute= " << mem_absolute 
            << "\n approx mem_per_node = " << mem_per_node << " [B] approx mem_per_element= " << mem_per_element << " [B]"
            << "\n" << std::endl;

  std::cout << "ref mesh.num_nodes()=  " << mesh.num_nodes() << std::endl;
  std::cout << "ref mesh.num_elements()=  " << mesh.num_elements() << std::endl;

  EXPECT_TRUE(mesh.num_nodes() == (NX+1)*(NY+1)+(NX*NY*2)*12);
  EXPECT_TRUE(mesh.num_elements() == 4*NX*NY*2);

  count_memory(mesh);
}

TEST(samba, tri_count_memory)
{
  samba::mesh mesh( samba::connectivity_map::default_map_2d());

  double mesh_create_time = -cpu_time();

  MemoryInfo mem_delta, mem_absolute;
  std::cout << "Memory: initial= " << mem_absolute << std::endl;

  samba::fixtures::tri_fixture tf(NX,NY);
  tf.samba_mesh_create(mesh);

  mesh_create_time += cpu_time();

  count_memory(mesh);
}
