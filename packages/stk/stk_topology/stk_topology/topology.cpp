#include <stk_topology/topology.hpp>
#include <ostream>
#include <sstream>
#include <iomanip>



namespace stk {

const char * topology::rank_names[] =
{
    "NODE_RANK"
  , "EDGE_RANK"
  , "FACE_RANK"
  , "ELEMENT_RANK"
  , "CONSTRAINT_RANK"
};

const char * topology::topology_names[] =
{
    "NODE"
  , "LINE_2"
  , "LINE_3"
  , "TRIANGLE_3"
  , "TRIANGLE_4"
  , "TRIANGLE_6"
  , "QUADRILATERAL_4"
  , "QUADRILATERAL_8"
  , "QUADRILATERAL_9"
  , "PARTICLE"
  , "LINE_2_1D"
  , "LINE_3_1D"
  , "BEAM_2"
  , "BEAM_3"
  , "SHELL_LINE_2"
  , "SHELL_LINE_3"
  , "TRIANGLE_3_2D"
  , "TRIANGLE_4_2D"
  , "TRIANGLE_6_2D"
  , "QUADRILATERAL_4_2D"
  , "QUADRILATERAL_8_2D"
  , "QUADRILATERAL_9_2D"
  , "SHELL_TRIANGLE_3"
  , "SHELL_TRIANGLE_4"
  , "SHELL_TRIANGLE_6"
  , "SHELL_QUADRILATERAL_4"
  , "SHELL_QUADRILATERAL_8"
  , "SHELL_QUADRILATERAL_9"
  , "TETRAHEDRON_4"
  , "TETRAHEDRON_8"
  , "TETRAHEDRON_10"
  , "TETRAHEDRON_11"
  , "PYRAMID_5"
  , "PYRAMID_13"
  , "PYRAMID_14"
  , "WEDGE_6"
  , "WEDGE_15"
  , "WEDGE_18"
  , "HEXAHEDRON_8"
  , "HEXAHEDRON_20"
  , "HEXAHEDRON_27"
  , "INVALID_TOPOLOGY"
};


std::ostream & operator<<(std::ostream &out, topology::rank_t r)
{
  if ( r < topology::END_RANK )
    return out << topology::rank_names[r];
  else if ( r < topology::INVALID_RANK)
    return out << "RANK_" << static_cast<int>(r);
  return out << "INVALID_RANK";
}

std::ostream & operator<<(std::ostream &out, topology t)
{
  if ( t <= topology::END_TOPOLOGY )
    return out << t.name();
  if ( t.is_superelement() )
    return out << "SUPERELEMENT_TOPOLOGY_" << (static_cast<int>(t) - topology::SUPERELEMENT_START);
  return out << "ARBITRARY_NODE_TOPOLOGY_" << (static_cast<int>(t) - topology::END_TOPOLOGY);
}

const char *topology::superelement_name() const
{
  static char basename[128];
  std::ostringstream oss;
  oss << *this;
  std::strncpy(basename, oss.str().c_str(), 128);
  basename[127] = '\0';
  return basename;
}
void verbose_print_topology(std::ostream &out, topology t)
{
  int shiftwidth = 34;

  int node_ordinals[27] = {0};

  out << std::boolalpha;

  out << t << std::endl;;
  out << std::setw(shiftwidth) << "is valid: " << t.is_valid() << std::endl;
  out << std::setw(shiftwidth) << "base: " << t.base() << std::endl;
  out << std::setw(shiftwidth) << "is shell: " << t.is_shell() << std::endl;
  out << std::setw(shiftwidth) << "rank: " << t.rank() << std::endl;
  out << std::setw(shiftwidth) << "side rank: " << t.side_rank() << std::endl;
  out << std::setw(shiftwidth) << "dimension: " << t.dimension() << std::endl;
  out << std::setw(shiftwidth) << "num nodes: " << t.num_nodes() << std::endl;
  out << std::setw(shiftwidth) << "num vertices: " << t.num_vertices() << std::endl;

  out << std::setw(shiftwidth) << "(1d, 2d, 3d): ";
  for (int i=1; i<4; ++i)
    out << t.defined_on_spatial_dimension(i) << ", ";
  out << "\b\b  " << std::endl;

  out << std::setw(shiftwidth) << "num edges: " << t.num_edges() << std::endl;
  if (t.num_edges() > 0) {
    const int num_edge_nodes = t.edge_topology().num_nodes();
    out << std::setw(shiftwidth) << t.edge_topology() << std::endl;
    for (int i=0, e=t.num_edges(); i<e; ++i) {
      out << std::setw(shiftwidth) << " " << i << ": (";
      t.edge_node_ordinals(i,node_ordinals);
      for (int j=0, ne = num_edge_nodes; j < ne; ++j) {
        out << node_ordinals[j] << ", ";
      }
      out << "\b\b)  " << std::endl;
    }
  }

  out << std::setw(shiftwidth) << "num faces: " << t.num_faces() << std::endl;
  if (t.num_faces() > 0) {
    for (int i=0, e=t.num_faces(); i<e; ++i) {
      out << std::setw(shiftwidth) << t.face_topology(i) << " " << i << ": (";
      t.face_node_ordinals(i,node_ordinals);
      for (int j=0, ne = t.face_topology(i).num_nodes(); j < ne; ++j) {
        out << node_ordinals[j] << ", ";
      }
      out << "\b\b)  " << std::endl;
    }
  }

  out << std::setw(shiftwidth) << "num permutations: " << t.num_permutations() << std::endl;
  out << std::setw(shiftwidth) << "num positive permutations: " << t.num_positive_permutations() << std::endl;
  if (t.num_permutations() > 0) {
    for (int i=0, e=t.num_positive_permutations(); i<e; ++i) {
      out << std::setw(shiftwidth) << i << ": (";
      t.permutation_node_ordinals(i,node_ordinals);
      for (int j=0, ne = t.num_nodes(); j < ne; ++j) {
        out << node_ordinals[j] << ", ";
      }
      out << "\b\b)  " << std::endl;
    }
    out << std::setw(shiftwidth) << "num negative permutations: " << t.num_permutations() - t.num_positive_permutations() << std::endl;
    if (t.num_positive_permutations() < t.num_permutations()) {
      for (int i=t.num_positive_permutations(), e=t.num_permutations(); i<e; ++i) {
        out << std::setw(shiftwidth) << i << ": (";
        t.permutation_node_ordinals(i,node_ordinals);
        for (int j=0, ne = t.num_nodes(); j < ne; ++j) {
          out << node_ordinals[j] << ", ";
        }
        out << "\b\b)  " << std::endl;
      }
    }
  }

  out << std::endl;
  out << std::noboolalpha;

};

} //namespace stk


