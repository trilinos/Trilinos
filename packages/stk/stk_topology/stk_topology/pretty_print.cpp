
#include <stk_topology/pretty_print.hpp>
#include <ostream>
#include <iomanip>

namespace stk {

std::ostream & operator<<(std::ostream &out, topology::rank_t r)
{
  return out << topology::rank_names[r];
}

std::ostream & operator<<(std::ostream &out, topology t)
{
  return out << t.name();
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
    for (int i=0, e=t.num_edges(); i<e; ++i) {
      out << std::setw(shiftwidth) << t.edge_topology(i) << " " << i << ": (";
      t.edge_node_ordinals(i,node_ordinals);
      for (int j=0, ne = t.edge_topology(i).num_nodes(); j < ne; ++j) {
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
  if (t.num_permutations() > 0) {
    for (int i=0, e=t.num_permutations(); i<e; ++i) {
      out << std::setw(shiftwidth) << i << ": (";
      t.permutation_node_ordinals(i,node_ordinals);
      for (int j=0, ne = t.num_nodes(); j < ne; ++j) {
        out << node_ordinals[j] << ", ";
      }
      out << "\b\b)  " << std::endl;
    }
  }

  out << std::endl;
  out << std::noboolalpha;

};

} //namespace stk
