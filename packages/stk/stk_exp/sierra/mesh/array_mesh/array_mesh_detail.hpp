#ifndef STK_SIERRA_MESH_ARRAYMESH_DETAIL_HPP
#define STK_SIERRA_MESH_ARRAYMESH_DETAIL_HPP

#include <string>
#include <vector>

namespace sierra {
namespace mesh {
namespace detail {

struct BlockIndex {
  int value;
  operator int() const { return value; }
};

struct SidesetIndex {
  int value;
  operator int() const { return value; }
};

struct NodesetIndex {
  int value;
  operator int() const { return value; }
};

struct Locator {
  BlockIndex block;
  unsigned offset_into_block;
};

struct Block {
  int id;
  std::string name;
  int topology;
  int rank;
  int num_elems;
  int nodes_per_elem;
  std::vector<int> elem_nums;
  std::vector<int> connectivity;
};

struct Sideset {
  int id;
  std::string name;
  std::vector<int> elem_nums;
  std::vector<int> elem_local_sides;
};

struct Nodeset {
	int id;
	std::string name;
    std::vector<int> node_nums;
};

}//namespace detail
}//namespace mesh
}//namespace sierra

#endif

