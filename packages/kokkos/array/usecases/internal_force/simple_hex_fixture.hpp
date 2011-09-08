#ifndef HEX_FIXTURE_HPP
#define HEX_FIXTURE_HPP

#include <vector>

namespace fixture {

/** Simplest possible flat-array mesh that can at least
  represent multiple element-blocks.
 This "mesh" is purely serial, and currently only represents element-node
  connectivity. Extra connectivity-tables could be added for other relations...
*/
struct simple_mesh {
  ~simple_mesh(){}

  void allocate(const std::vector<std::pair<int,int> >& elem_blks, const int num_nodes)
  {
    elem_blocks = elem_blks;
    int total_num_elems = 0;
    int total_elem_node_connectivity_length = 0;
    for(size_t i=0; i<elem_blks.size(); ++i) {
      int num_elems = elem_blks[i].first;
      int nodes_per_elem = elem_blks[i].second;
      total_num_elems += num_elems;
      total_elem_node_connectivity_length += num_elems*nodes_per_elem;
    }
    elem_node_connectivity.resize(total_elem_node_connectivity_length);

    node_elem_offset.resize(num_nodes);
    node_elem_ids.reserve(num_nodes*8);

  }

  /** The 'elem_blocks' vector of pairs holds, for each element-block,
     the number of elements in the block and the number of nodes per element.
  */
  std::vector<std::pair<int,int> > elem_blocks;

  /** elem_node_connectivity will have length:
         sum(elem_blocks[i].first*elem_blocks[i].second).
  */
  std::vector<int> elem_node_connectivity;

  std::vector<int> node_elem_offset;
  std::vector<int> node_elem_ids;
};

/**
 * A 3-dimensional X*Y*Z hex fixture.
 *
 * A coordinate field will be created for all nodes.
 */
class simple_hex_fixture
{
 public:
  /**
    The parameters nx, ny, nz specify the number of elements in each
    spatial dimension. Thus, the number of nodes in each dimension will
    be 1 greater (nx+1, ny+1, nz+1).
   */
  simple_hex_fixture(unsigned nx, unsigned ny, unsigned nz);

  const unsigned                m_nx;
  const unsigned                m_ny;
  const unsigned                m_nz;
  const unsigned                m_num_nodes;
  const unsigned                m_num_elements;
  simple_mesh                   m_mesh;
  std::vector<double>           m_coord_field;

  /**
   * Thinking in terms of a 3D grid of nodes, get the index of the node in
   * the (x, y, z) position.
   */
  int node_index( unsigned x , unsigned y , unsigned z ) const  {
    return x + ( m_nx + 1 ) * ( y + ( m_ny + 1 ) * z );
  }

  /**
   * Thinking in terms of a 3D grid of elements, get the id of the
   * element in the (x, y, z) position.
   */
  int elem_index( unsigned x , unsigned y , unsigned z ) const  {
    return m_num_nodes + x + m_nx * ( y + m_ny * z );
  }

  /**
   * Thinking in terms of a 3D grid of nodes, compute the (x, y, z) position
   * of a node given its index.
   */
  void node_x_y_z( int entity_key, unsigned &x , unsigned &y , unsigned &z ) const;

  /**
   * Thinking in terms of a 3D grid of elements, compute the (x, y, z) position
   * of an element given it's id.
   */
  void elem_x_y_z( int entity_key, unsigned &x , unsigned &y , unsigned &z ) const;

  /**
   * Create the mesh (into m_mesh).
   */
  void generate_mesh();

  void populate_field_data();

 private:

  simple_hex_fixture();
  simple_hex_fixture( const simple_hex_fixture &);
  simple_hex_fixture & operator = (const simple_hex_fixture &);
};

} // fixture

#endif

