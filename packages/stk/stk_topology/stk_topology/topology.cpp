// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stk_topology/topology.hpp>
#include <ostream>
#include <sstream>
#include <iomanip>



namespace stk {

std::string topology::name() const
{
  switch (m_value)
  {
  case INVALID_TOPOLOGY:         return "INVALID_TOPOLOGY";
  case NODE:             return "NODE";
  case LINE_2:           return "LINE_2";
  case LINE_3:           return "LINE_3";
  case TRI_3:            return "TRIANGLE_3";
  case TRI_4:            return "TRIANGLE_4";
  case TRI_6:            return "TRIANGLE_6";
  case QUAD_4:           return "QUADRILATERAL_4";
  case QUAD_8:           return "QUADRILATERAL_8";
  case QUAD_9:           return "QUADRILATERAL_9";
  case PARTICLE:         return "PARTICLE";
  case LINE_2_1D:        return "LINE_2_1D";
  case LINE_3_1D:        return "LINE_3_1D";
  case BEAM_2:           return "BEAM_2";
  case BEAM_3:           return "BEAM_3";
  case SHELL_LINE_2:     return "SHELL_LINE_2";
  case SHELL_LINE_3:     return "SHELL_LINE_3";
  case TRI_3_2D:         return "TRIANGLE_3_2D";
  case TRI_4_2D:         return "TRIANGLE_4_2D";
  case TRI_6_2D:         return "TRIANGLE_6_2D";
  case QUAD_4_2D:        return "QUADRILATERAL_4_2D";
  case QUAD_8_2D:        return "QUADRILATERAL_8_2D";
  case QUAD_9_2D:        return "QUADRILATERAL_9_2D";
  case SHELL_TRI_3:      return "SHELL_TRIANGLE_3";
  case SHELL_TRI_4:      return "SHELL_TRIANGLE_4";
  case SHELL_TRI_6:      return "SHELL_TRIANGLE_6";
  case SHELL_QUAD_4:     return "SHELL_QUADRILATERAL_4";
  case SHELL_QUAD_8:     return "SHELL_QUADRILATERAL_8";
  case SHELL_QUAD_9:     return "SHELL_QUADRILATERAL_9";
  case TET_4:            return "TETRAHEDRON_4";
  case TET_8:            return "TETRAHEDRON_8";
  case TET_10:           return "TETRAHEDRON_10";
  case TET_11:           return "TETRAHEDRON_11";
  case PYRAMID_5:        return "PYRAMID_5";
  case PYRAMID_13:       return "PYRAMID_13";
  case PYRAMID_14:       return "PYRAMID_14";
  case WEDGE_6:          return "WEDGE_6";
  case WEDGE_15:         return "WEDGE_15";
  case WEDGE_18:         return "WEDGE_18";
  case HEX_8:            return "HEXAHEDRON_8";
  case HEX_20:           return "HEXAHEDRON_20";
  case HEX_27:           return "HEXAHEDRON_27";
  default: break;
  }

// jvo: Superelements require additional calculations? What is value calculated from m_value and SUPERELEMENT_START?
  std::ostringstream oss;
  if ( is_superelement() )
    oss << "SUPERELEMENT_TOPOLOGY_" << (static_cast<unsigned>(m_value) - topology::SUPERELEMENT_START);
  else
    oss << "UNKNOWN_TOPOLOGY_" << (static_cast<unsigned>(m_value));

  return oss.str();
}

// jvo: Overloading << operator?
std::ostream & operator<<(std::ostream &out, topology::rank_t r)
{
  switch (r)
  {
  case topology::NODE_RANK:    out << "NODE_RANK"; break;
  case topology::EDGE_RANK:    out << "EDGE_RANK"; break;
  case topology::FACE_RANK:    out << "FACE_RANK"; break;
  case topology::ELEMENT_RANK: out << "ELEMENT_RANK"; break;
  case topology::INVALID_RANK: out << "INVALID_RANK"; break;
  default:                     out << "RANK_" << static_cast<unsigned>(r); break;
  }
  return out;
}

std::ostream & operator<<(std::ostream &out, topology t)
{
  return out << t.name();
}


void verbose_print_topology(std::ostream &out, topology t)
{
  unsigned shiftwidth = 34;
  // jvo: Is 27 the max number of nodes for any topology? What about superelements?
  unsigned node_ordinals[27] = {0};

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
  for (unsigned i=1; i<4; ++i)
    out << t.defined_on_spatial_dimension(i) << ", ";
  out << "\b\b  " << std::endl;

  out << std::setw(shiftwidth) << "num edges: " << t.num_edges() << std::endl;
  if (t.num_edges() > 0) {
    const unsigned num_edge_nodes = t.edge_topology().num_nodes();
    out << std::setw(shiftwidth) << t.edge_topology() << std::endl;
    for (unsigned i=0, e=t.num_edges(); i<e; ++i) {
      out << std::setw(shiftwidth) << " " << i << ": (";
      t.edge_node_ordinals(i,node_ordinals);
      for (unsigned j=0, ne = num_edge_nodes; j < ne; ++j) {
        out << node_ordinals[j] << ", ";
      }
      out << "\b\b)  " << std::endl;
    }
  }

  out << std::setw(shiftwidth) << "num faces: " << t.num_faces() << std::endl;
  if (t.num_faces() > 0) {
    for (unsigned i=0, e=t.num_faces(); i<e; ++i) {
      out << std::setw(shiftwidth) << t.face_topology(i) << " " << i << ": (";
      t.face_node_ordinals(i,node_ordinals);
      for (unsigned j=0, ne = t.face_topology(i).num_nodes(); j < ne; ++j) {
        out << node_ordinals[j] << ", ";
      }
      out << "\b\b)  " << std::endl;
    }
  }

  // jvo: is positive permutation according to right-hand-rule?
  out << std::setw(shiftwidth) << "num permutations: " << t.num_permutations() << std::endl;
  out << std::setw(shiftwidth) << "num positive permutations: " << t.num_positive_permutations() << std::endl;
  if (t.num_permutations() > 0) {
    for (unsigned i=0, e=t.num_positive_permutations(); i<e; ++i) {
      out << std::setw(shiftwidth) << i << ": (";
      t.permutation_node_ordinals(i,node_ordinals);
      for (unsigned j=0, ne = t.num_nodes(); j < ne; ++j) {
        out << node_ordinals[j] << ", ";
      }
      out << "\b\b)  " << std::endl;
    }
    out << std::setw(shiftwidth) << "num negative permutations: " << t.num_permutations() - t.num_positive_permutations() << std::endl;
    if (t.num_positive_permutations() < t.num_permutations()) {
      for (unsigned i=t.num_positive_permutations(), e=t.num_permutations(); i<e; ++i) {
        out << std::setw(shiftwidth) << i << ": (";
        t.permutation_node_ordinals(i,node_ordinals);
        for (unsigned j=0, ne = t.num_nodes(); j < ne; ++j) {
          out << node_ordinals[j] << ", ";
        }
        out << "\b\b)  " << std::endl;
      }
    }
  }

  out << std::endl;
  out << std::noboolalpha;

}

} //namespace stk


