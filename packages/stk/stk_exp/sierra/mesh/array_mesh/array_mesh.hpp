/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_SIERRA_MESH_ARRAYMESH_HPP
#define STK_SIERRA_MESH_ARRAYMESH_HPP

#include <stk_topology/topology.hpp>
#include <sierra/mesh/array_mesh/array_utils.hpp>

#include <sierra/mesh/array_mesh/array_mesh_detail.hpp>

#include <vector>
#include <iterator>
#include <stdexcept>

namespace sierra {
namespace mesh {

/** 'array_mesh' stores everything in arrays (actually std::vectors...)
  *  Very similar to an in-memory version of the exodus format.
  *
  * array_mesh essentially stores the following kinds of data "collections":
  *
  *   1. blocks (usually element-blocks but other blocks are possible)
  *   2. sidesets
  *   3. nodesets
  *
  *   See the exodus manual for descriptions of these.
  *
  *
  * Note: some array_mesh methods reference node/element 'numbers' and 'ids'.
  *       These follow the exodus convention that a node-number is an index into
  *       the local array of nodes, while a node-id is a globally unique identifier.
  *
  *  array_mesh is currently purely serial. No knowledge of sharing etc.
  *
  * The best documentation for this class is the gather-centroid program, located
  * here: stk/performance_tests/mesh/array_mesh_hex_gather.cpp
  * Also see the I/O example: stk/examples/io/io_array_mesh.cpp
  */
class array_mesh {
 public:

  static const int Node       = stk::topology::NODE_RANK;
  static const int Element    = stk::topology::ELEMENT_RANK;
  static const int Num_ranks  = stk::topology::NUM_RANKS;

  typedef std::vector<int>::const_iterator             ConstIntIterator;
  typedef std::pair<ConstIntIterator,ConstIntIterator> ConstIntRange;

  typedef detail::BlockIndex                       BlockIndex;
  typedef std::vector<BlockIndex>::const_iterator  BlockIterator;
  typedef std::pair<BlockIterator,BlockIterator>   BlockRange;

  typedef detail::SidesetIndex                       SidesetIndex;
  typedef std::vector<SidesetIndex>::const_iterator  SidesetIterator;
  typedef std::pair<SidesetIterator,SidesetIterator> SidesetRange;

  typedef detail::NodesetIndex                       NodesetIndex;
  typedef std::vector<NodesetIndex>::const_iterator  NodesetIterator;
  typedef std::pair<NodesetIterator,NodesetIterator> NodesetRange;

  /** Constructor.
     Creation of upward connectivity is optional. Currently upward connections are
     only created from nodes to elements.
  */
  array_mesh(bool create_upward_connectivity = false)
   : m_blocks(),
     m_block_indices(),
     m_sidesets(),
     m_sideset_indices(),
     m_nodesets(),
     m_nodeset_indices(),
     m_global_ids(Num_ranks),
     m_create_upward_connectivities(create_upward_connectivity),
     m_upward_connectivities(Num_ranks),
     m_empty()
  {
  }

  ~array_mesh(){}

  /** reserve makes later node/element additions more efficient, similar to std::vector::reserve.
  */
  void reserve(size_t num_elems, size_t num_nodes);

  size_t get_num_nodes() const { return m_global_ids[Node].size(); }
  size_t get_num_elements() const { return m_global_ids[Element].size(); }

  /** add global-ids for nodes. This method pushes-back onto the internal node-map.
    So don't add repeated nodes...
  */
  template<class IdIterator>
  void add_node_ids(IdIterator id_begin, IdIterator id_end);

  template<class IdIterator>
  void add_element_ids(IdIterator id_begin, IdIterator id_end);

  int get_node_id(int node_num) const;
  int get_element_id(int elem_num) const;

  /** a ConstIntRange is just a pair of vector<int>::const_iterator
  */
  ConstIntRange get_node_ids() const
  { return std::make_pair(m_global_ids[stk::topology::NODE_RANK].begin(),m_global_ids[stk::topology::NODE_RANK].end()); }

  ConstIntRange get_element_ids() const
  { return std::make_pair(m_global_ids[Element].begin(),m_global_ids[Element].end()); }


  ///////////////// block methods /////////////////////

  /** Declare the existence and attributes of a block. (Usually an element-block, but edge/face blocks
    * could fit into this scheme also.
    * TODO: should we guard against creating blocks with already-existing names/ids?
  */
  BlockIndex add_block( int rank,
                        int ID,
                        int num_elems,
                        stk::topology T,
                        const std::string& name = std::string(""));

  /** Add connectivity to an already-declared block.
  */
  template<class NodeIterator>
  void add_connectivity(BlockIndex block_index, int elem_num, NodeIterator nodes_begin, NodeIterator nodes_end);

  BlockIndex get_block(const std::string& name) const;
  BlockIndex get_block(int ID) const;
  BlockIndex get_element_block(int ID) const;

  /** a BlockRange is a pair of vector<BlockIndex>::const_iterator
  */
  BlockRange get_blocks() const
  { return std::make_pair(m_block_indices.begin(), m_block_indices.end()); }

  size_t get_num_node_blocks() const
  {
    // Assumes 1 node block in exodus
    return 1;
  }

  size_t get_num_element_blocks() const
  { return get_num_blocks() - get_num_node_blocks(); }

  size_t get_num_blocks() const
  { return m_block_indices.size(); }

  /** getters for block attributes.
    TODO: Document these getters.
  */
  const std::string& get_name(BlockIndex block_index) const { return m_blocks[block_index].name; }
  int get_id(BlockIndex block_index) const { return m_blocks[block_index].id; }
  stk::topology get_topology(BlockIndex block_index) const { return m_blocks[block_index].topology; }
  int get_rank(BlockIndex block_index) const { return m_blocks[block_index].rank; }
  int get_num_elems(BlockIndex block_index) const { return m_blocks[block_index].num_elems; }
  int get_num_nodes_per_elem(BlockIndex block_index) const { return m_blocks[block_index].nodes_per_elem; }
  const std::vector<int>& get_block_elem_nums(BlockIndex block_index) const { return m_blocks[block_index].elem_nums; }
  const std::vector<int>& get_block_connectivity(BlockIndex block_index) const { return m_blocks[block_index].connectivity; }

  /////////////// end of block methods ////////////////////


  /////////////// sideset methods //////////////////////////

  /** create new sideset.
   * *TODO Should we guard against creating sidesets with already-existing names/ids?
   */
  SidesetIndex add_sideset(int ID, const std::string& name = std::string(""));
  void add_side(SidesetIndex sideset_index, int elem_num, int elem_local_side);

  SidesetIndex get_sideset(const std::string& name) const;
  SidesetIndex get_sideset(int ID) const;

  size_t get_num_sidesets() const
  { return m_sidesets.size(); }

  SidesetRange get_sidesets() const;
  const std::string& get_name(SidesetIndex sideset_index) const { return m_sidesets[sideset_index].name; }
  int get_id(SidesetIndex sideset_index) const {return m_sidesets[sideset_index].id; }
  const std::vector<int>& get_sideset_elem_nums(SidesetIndex sideset_index) const;
  const std::vector<int>& get_sideset_elem_local_sides(SidesetIndex sideset_index) const;

  unsigned get_num_side_nodes(int elem_num, int elem_local_side) const;
  void get_side_nodes(int elem_num, int elem_local_side, std::vector<int>& node_nums) const;
  void get_side_nodes(SidesetIndex sideset_index, int side, std::vector<int>& node_nums) const;

  ////////////// end of sideset methods /////////////////

  /////////////// nodeset methods ////////////////////////

  /** create new nodeset
   */
  NodesetIndex add_nodeset(int ID, const std::string& name);
  void add_node(NodesetIndex, int node_num);
  template<typename NodeIterator>
  void add_nodes(NodesetIndex, NodeIterator nodes_begin, NodeIterator nodes_end);

  NodesetRange get_nodesets() const;
  NodesetIndex get_nodeset(const std::string& name) const;
  NodesetIndex get_nodeset(int ID) const;
  size_t get_size(NodesetIndex nodeset_index) const;
  ConstIntRange get_nodes(NodesetIndex nodeset_index) const;
  const std::string& get_name(NodesetIndex nodeset_index) const;
  int get_id(NodesetIndex nodeset_index) const;

  //////////////end of nodeset methods ///////////////////

  /** getter for the nodes connected to a given element.
    It is more efficient to iterate element-node connectivity by iterating a block's connectivity.
  */
  ConstIntRange get_connected_nodes(int elem_num) const;

  stk::topology get_element_topology(int elem_num) const;

  /** Upward connectivity: node-to-element.
    Doesn't exist unless the optional constructor argument was used to enable upward connectivity.
  */
  ConstIntRange get_connected_elements(int node_num) const;

 private:
  typedef std::vector<int> id_vector;
  typedef std::vector<std::vector<int> > upward_connectivity;

  std::vector<detail::Block> m_blocks;
  std::vector<BlockIndex> m_block_indices;
  std::vector<detail::Sideset> m_sidesets;
  std::vector<SidesetIndex> m_sideset_indices;
  std::vector<detail::Nodeset> m_nodesets;
  std::vector<NodesetIndex> m_nodeset_indices;

  std::vector<id_vector> m_global_ids;
  std::vector<detail::Locator> m_elem_locations;

  bool m_create_upward_connectivities;
  std::vector<upward_connectivity> m_upward_connectivities;

  std::vector<int> m_empty;
};

inline
void array_mesh::reserve(size_t num_elems, size_t  num_nodes)
{
  m_global_ids[Node].reserve(num_nodes);

  m_global_ids[Element].reserve(num_elems);
  m_elem_locations.reserve(num_elems);

  if(m_create_upward_connectivities) {
    m_upward_connectivities[Node].reserve(num_nodes);
  }
}

inline
array_mesh::BlockIndex array_mesh::add_block(int rank,
                                             int ID,
                                             int num_elems,
                                             stk::topology T,
                                             const std::string& name)
{

  detail::Block blk;
  blk.id = ID;
  blk.name = name;
  blk.topology = T;
  blk.rank = rank;
  blk.num_elems = num_elems;
  int nodes_per_elem = T.num_nodes();
  blk.nodes_per_elem = nodes_per_elem;

  BlockIndex block_index;
  block_index.value = m_blocks.size();
  m_blocks.push_back(blk);

  m_blocks[block_index].elem_nums.reserve(num_elems);
  m_blocks[block_index].connectivity.reserve(num_elems*nodes_per_elem);

  m_block_indices.push_back(block_index);

  return block_index;
}

template<class NodeIterator>
inline
void array_mesh::add_connectivity(BlockIndex block_index, int elem_num,
                                  NodeIterator nodes_begin, NodeIterator nodes_end)
{
  std::vector<int>& conn = m_blocks[block_index].connectivity;

  detail::Locator loc;
  loc.block = block_index;
  loc.offset_into_block = m_blocks[block_index].elem_nums.size();
  stk::topology topo = m_blocks[block_index].topology;
  if ( topo != stk::topology::NODE )
  {
    m_elem_locations.push_back(loc);
  }

  m_blocks[block_index].elem_nums.push_back(elem_num);

  for(NodeIterator nd_iter=nodes_begin; nd_iter!=nodes_end; ++nd_iter) {
    conn.push_back(*nd_iter);

    if (m_create_upward_connectivities) {
      size_t this_node_num = *nd_iter;
      if (m_upward_connectivities[Node].size() <= this_node_num) {
        m_upward_connectivities[Node].resize(this_node_num+1);
      }

      std::vector<int>& node_elems = m_upward_connectivities[Node][this_node_num];
      insert_sorted(node_elems, elem_num);
    }
  }
}

template<class IdIterator>
inline
void array_mesh::add_node_ids(IdIterator id_begin, IdIterator id_end)
{
  for(IdIterator id_iter=id_begin; id_iter != id_end; ++id_iter) {
    m_global_ids[Node].push_back(*id_iter);
  }
}

template<class IdIterator>
inline
void array_mesh::add_element_ids(IdIterator id_begin, IdIterator id_end)
{
  for(IdIterator id_iter=id_begin; id_iter != id_end; ++id_iter) {
    m_global_ids[Element].push_back(*id_iter);
  }
}

inline
int array_mesh::get_node_id(int node_idx) const
{ return m_global_ids[Node][node_idx]; }

inline
int array_mesh::get_element_id(int elem_idx) const
{ return m_global_ids[Element][elem_idx]; }

inline
array_mesh::ConstIntRange array_mesh::get_connected_nodes(int elem_num) const
{
  const detail::Locator& loc = m_elem_locations[elem_num];
  const std::vector<int>& conn = m_blocks[loc.block].connectivity;
  const int nnodes = m_blocks[loc.block].nodes_per_elem;
  const int offset = loc.offset_into_block*nnodes;

  return std::make_pair(conn.begin()+offset, conn.begin()+offset+nnodes);
}

inline
stk::topology array_mesh::get_element_topology(int elem_num) const
{
  const detail::Locator& loc = m_elem_locations[elem_num];
  return m_blocks[loc.block].topology;
}

inline
array_mesh::ConstIntRange array_mesh::get_connected_elements(int node_num) const
{
  if (m_create_upward_connectivities) {
    const std::vector<int>& conn = m_upward_connectivities[Node][node_num];
    return std::make_pair(conn.begin(), conn.end());
  }

  return std::make_pair(m_empty.end(), m_empty.end());
}

inline
array_mesh::BlockIndex array_mesh::get_block(const std::string& name) const
{
	BlockIndex bidx = {-1};
	for(size_t i=0; i<m_blocks.size(); ++i) {
		if (m_blocks[i].name == name) {
			bidx.value = i;
			return bidx;
		}
	}
	return bidx;
}

inline
array_mesh::BlockIndex array_mesh::get_block(int ID) const
{
	BlockIndex bidx = {-1};
	for(size_t i=0; i<m_blocks.size(); ++i) {
		if (m_blocks[i].id == ID) {
			bidx.value = i;
			return bidx;
		}
	}
	return bidx;
}

inline
array_mesh::BlockIndex array_mesh::get_element_block(int ID) const
{
    BlockIndex bidx = {-1};
    for(size_t i = 0; i < m_blocks.size(); ++i)
    {
        if(m_blocks[i].id == ID && m_blocks[i].rank == Element)
        {
            bidx.value = i;
            return bidx;
        }
    }
    return bidx;
}

inline
array_mesh::SidesetIndex array_mesh::add_sideset(int ID, const std::string& name)
{
  int ss_idx = m_sidesets.size();
  SidesetIndex sideset_index = { ss_idx };
  detail::Sideset sideset;
  sideset.id = ID;
  sideset.name = name;
  m_sidesets.push_back(sideset);
  m_sideset_indices.push_back(sideset_index);
  return sideset_index;
}

inline
void array_mesh::add_side(SidesetIndex side_index, int elem_num, int elem_local_side)
{
  detail::Sideset& sideset = m_sidesets[side_index];
  sideset.elem_nums.push_back(elem_num);
  sideset.elem_local_sides.push_back(elem_local_side);
}

inline
array_mesh::SidesetIndex array_mesh::get_sideset(const std::string& name) const
{
	SidesetIndex sidx = {-1};
	for(size_t i=0; i<m_sidesets.size(); ++i) {
		if (m_sidesets[i].name == name) {
			sidx.value = i;
			return sidx;
		}
	}
	return sidx;
}

inline
array_mesh::SidesetIndex array_mesh::get_sideset(int ID) const
{
	SidesetIndex sidx = {-1};
	for(size_t i=0; i<m_sidesets.size(); ++i) {
		if (m_sidesets[i].id == ID) {
			sidx.value = i;
			return sidx;
		}
	}
	return sidx;
}

inline
array_mesh::SidesetRange array_mesh::get_sidesets() const
{
  return std::make_pair( m_sideset_indices.begin(), m_sideset_indices.end() );
}

inline
const std::vector<int>& array_mesh::get_sideset_elem_nums(SidesetIndex sideset_index) const
{
  return m_sidesets[sideset_index].elem_nums;
}

inline
const std::vector<int>& array_mesh::get_sideset_elem_local_sides(SidesetIndex sideset_index) const
{
  return m_sidesets[sideset_index].elem_local_sides;
}

inline
unsigned array_mesh::get_num_side_nodes(int elem_num, int elem_local_side) const
{
  stk::topology topo = get_element_topology(elem_num);
  stk::topology side_topo = topo.side_topology(elem_local_side);

  return side_topo.num_nodes();
}

inline
void array_mesh::get_side_nodes(SidesetIndex sideset_index, int side, std::vector<int>& node_nums) const
{
  int elem_num = m_sidesets[sideset_index].elem_nums[side];
  stk::topology elem_topo = get_element_topology(elem_num);
  int elem_local_side = m_sidesets[sideset_index].elem_local_sides[side];
  get_side_nodes(elem_num, elem_local_side, node_nums);
}

inline
void array_mesh::get_side_nodes(int elem_num, int elem_local_side, std::vector<int>& node_nums) const
{
  node_nums.clear();
  ConstIntRange nodes = get_connected_nodes(elem_num);
  ConstIntIterator node_iterator = nodes.first;

  stk::topology elem_topo = get_element_topology(elem_num);
  std::vector<int> side_node_vector;
  elem_topo.face_nodes(node_iterator, elem_local_side, std::back_inserter(side_node_vector));
  for(unsigned i=0; i<get_num_side_nodes(elem_num, elem_local_side); ++i) {
    node_nums.push_back(side_node_vector[i]);
  }
}

inline
array_mesh::NodesetIndex array_mesh::add_nodeset(int ID, const std::string& name)
{
  int ns_idx = m_nodesets.size();
	NodesetIndex nodeset_index = {ns_idx};
	detail::Nodeset nodeset;
	m_nodesets.push_back(nodeset);
	m_nodeset_indices.push_back(nodeset_index);
	m_nodesets[nodeset_index].id = ID;
	m_nodesets[nodeset_index].name = name;
	return nodeset_index;
}

inline
 void array_mesh::add_node(NodesetIndex nodeset_index, int node_num)
{
	m_nodesets[nodeset_index].node_nums.push_back(node_num);
}

 template<typename NodeIterator>
inline
void array_mesh::add_nodes(NodesetIndex nodeset_index, NodeIterator nodes_begin, NodeIterator nodes_end)
 {
	 for(NodeIterator nodes_iter=nodes_begin; nodes_iter!=nodes_end; ++nodes_iter) {
		 m_nodesets[nodeset_index].node_nums.push_back(*nodes_iter);
	 }
 }

inline
array_mesh::NodesetRange array_mesh::get_nodesets() const
{
	return std::make_pair(m_nodeset_indices.begin(), m_nodeset_indices.end());
}

inline
array_mesh::NodesetIndex array_mesh::get_nodeset(const std::string& name) const
{
	NodesetIndex nodeset_index = {-1};
	for(size_t i=0; i<m_nodesets.size(); ++i) {
		if (m_nodesets[i].name == name) {
			nodeset_index.value = i;
			break;
		}
	}
	return nodeset_index;
}

inline
array_mesh::NodesetIndex array_mesh::get_nodeset(int ID) const
{
	NodesetIndex nodeset_index = {-1};
	for(size_t i=0; i<m_nodesets.size(); ++i) {
		if (m_nodesets[i].id == ID) {
			nodeset_index.value = i;
			break;
		}
	}
	return nodeset_index;
}

inline
size_t array_mesh::get_size(NodesetIndex nodeset_index) const
{
	if (nodeset_index.value < 0) throw std::runtime_error("array_mesh::get_size(nodeset_index) ERROR, bad nodeset_index");
	return m_nodesets[nodeset_index].node_nums.size();
}

inline
array_mesh::ConstIntRange array_mesh::get_nodes(NodesetIndex nodeset_index) const
{
	if (nodeset_index.value < 0) return std::make_pair(m_empty.end(), m_empty.end());
	return std::make_pair(m_nodesets[nodeset_index].node_nums.begin(),
							m_nodesets[nodeset_index].node_nums.end());
}

inline
const std::string& array_mesh::get_name(NodesetIndex nodeset_index) const
{
	if (nodeset_index.value < 0) throw std::runtime_error("array_mesh::get_name ERROR, bad nodeset_index");
	return m_nodesets[nodeset_index].name;
}

inline
int array_mesh::get_id(NodesetIndex nodeset_index) const
{
	if (nodeset_index.value < 0) throw std::runtime_error("array_mesh::get_name ERROR, bad nodeset_index");
	return m_nodesets[nodeset_index].id;
}

} // mesh
} // sierra

#endif

