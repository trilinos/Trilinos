#include <gtest/gtest.h>

#include <samba_fixtures/tri_fixture.hpp>
#include <samba/mesh.hpp>

#include <iostream>
#include <iterator>
#include <algorithm>
#include <sstream>
#include <vector>

#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/equal.hpp>

typedef uint64_t MemorySizeType;

#define DO_MEMORY_ACCOUNTING 0

#if DO_MEMORY_ACCOUNTING
struct MemoryInfo
{
  MemorySizeType malloc_used;
  MemorySizeType malloc_footprint;
  MemorySizeType malloc_max_footprint;
  static const MemorySizeType MB = 1024*1024;

  void get_memory_usage() 
  {
#if !defined(SIERRA_PTMALLOC3_ALLOCATOR) && !defined(SIERRA_PTMALLOC2_ALLOCATOR)
          std::cout << "WARNING: ptmalloc2|3 not compiled in so malloc_used info unavailable.  Recompile with e.g. 'bake allocator=ptmalloc2 (or 3)'.  Printing zeros..." << std::endl;
#else

#endif
    malloc_used = malloc_used();
    malloc_footprint = malloc_footprint();
    malloc_max_footprint = malloc_max_footprint();
  }
  void set_state() { get_memory_usage(); }
  void get_increment() {
    MemoryInfo old_state = *this;
    get_memory_usage();
    malloc_used -= old_state.malloc_used;
    malloc_footprint -= old_state.malloc_footprint;
  }

};

std::ostream& operator<<(std::ostream& os, const MemoryInfo& mem)
{
  os << "malloc_used= " << mem.malloc_used/mem.MB << " [MB] malloc_footprint= " << mem.malloc_footprint/mem.MB << " [MB] malloc_max_footprint= " << malloc_max_footprint/mem.MB << " [MB]";
  return os;
}
#endif

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

TEST(samba, tri_refine)
{
  samba::mesh mesh( samba::connectivity_map::default_map_2d());

  samba::fixtures::tri_fixture tf;
  tf.samba_mesh_create(mesh);

  std::cout << "orig mesh.num_nodes()=  " << mesh.num_nodes() << std::endl;
  std::cout << "orig mesh.num_elements()=  " << mesh.num_elements() << std::endl;
  EXPECT_TRUE(mesh.num_nodes() == 121);
  EXPECT_TRUE(mesh.num_elements() == 200);
  refine_tri_mesh(mesh);
  std::cout << "ref mesh.num_nodes()=  " << mesh.num_nodes() << std::endl;
  std::cout << "ref mesh.num_elements()=  " << mesh.num_elements() << std::endl;
  EXPECT_TRUE(mesh.num_nodes() == 2521);
  EXPECT_TRUE(mesh.num_elements() == 800);

}
