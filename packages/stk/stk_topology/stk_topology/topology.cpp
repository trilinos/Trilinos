// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
  if (m_value < END_TOPOLOGY) {
    return char_name();
  }

  std::ostringstream oss;
  oss << char_name() << "_";
  if ( is_superelement() )
    oss << (static_cast<unsigned>(m_value) - topology::SUPERELEMENT_START);
  else if ( is_superface() )
    oss << (static_cast<unsigned>(m_value) - topology::SUPERFACE_START);
  else if ( is_superedge() )
    oss << (static_cast<unsigned>(m_value) - topology::SUPEREDGE_START);
  else
    oss << (static_cast<unsigned>(m_value));

  return oss.str();
}

STK_FUNCTION
const char * topology::char_name() const
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
  case QUAD_6:           return "QUADRILATERAL_6";
  case QUAD_8:           return "QUADRILATERAL_8";
  case QUAD_9:           return "QUADRILATERAL_9";
  case PARTICLE:         return "PARTICLE";
  case LINE_2_1D:        return "LINE_2_1D";
  case LINE_3_1D:        return "LINE_3_1D";
  case BEAM_2:           return "BEAM_2";
  case BEAM_3:           return "BEAM_3";
  case SHELL_LINE_2:     return "SHELL_LINE_2";
  case SHELL_LINE_3:     return "SHELL_LINE_3";
  case SPRING_2:         return "SPRING_2";
  case SPRING_3:         return "SPRING_3";
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
  case WEDGE_12:         return "WEDGE_12";
  case WEDGE_15:         return "WEDGE_15";
  case WEDGE_18:         return "WEDGE_18";
  case HEX_8:            return "HEXAHEDRON_8";
  case HEX_20:           return "HEXAHEDRON_20";
  case HEX_27:           return "HEXAHEDRON_27";
  default: break;
  }

  if ( is_superelement() ) {
    return "SUPERELEMENT_TOPOLOGY";
  }
  else if ( is_superface() ) {
    return "SUPERFACE_TOPOLOGY";
  }
  else if ( is_superedge() ) {
    return "SUPEREDGE_TOPOLOGY";
  }

  return "UNKNOWN_TOPOLOGY";
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

bool isTriangleElement (topology topo)
{
    return ((topo == topology::TRI_3_2D) || (topo == topology::TRI_4_2D) || (topo == topology::TRI_6_2D));
}

bool isQuadrilateralElement (topology topo)
{
    return ((topo == topology::QUAD_4_2D) || (topo == topology::QUAD_8_2D) || (topo == topology::QUAD_9_2D));
}

bool isTetrahedronElement (topology topo)
{
    return ((topo == topology::TET_4) || (topo == topology::TET_8) || (topo == topology::TET_10) || (topo == topology::TET_11));
}

bool isHexahedronElement (topology topo)
{
    return ((topo == topology::HEX_8) || (topo == topology::HEX_20) || (topo == topology::HEX_27));
}

bool is_quad_side(topology topo)
{
    return ((topo == topology::QUAD_4) || (topo == topology::QUAD_6) || (topo == topology::QUAD_8) || (topo == topology::QUAD_9));
}

bool is_tri_side(topology topo)
{
    return ((topo == topology::TRI_3) || (topo == topology::TRI_4) || (topo == topology::TRI_6));
}

void verbose_print_topology(std::ostream &out, topology t)
{
  unsigned shiftwidth = 34;
  // jvo: Is 27 the max number of nodes for any topology? What about superelements?
  unsigned node_ordinals[27] = {0};

  out << std::boolalpha;

  out << t << std::endl;
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

  unsigned numEdges = t.num_edges();
  out << std::setw(shiftwidth) << "num edges: " << numEdges << std::endl;
  if (numEdges > 0) {
    for (unsigned i=0; i<numEdges; ++i) {
      const unsigned num_edge_nodes = t.edge_topology(i).num_nodes();
      out << std::setw(shiftwidth) << " " << t.edge_topology(i);
      out << std::setw(shiftwidth) << " " << i << ": (";
      t.edge_node_ordinals(i,node_ordinals);
      for (unsigned j=0, ne = num_edge_nodes; j < ne; ++j) {
        out << node_ordinals[j] << ", ";
      }
      out << "\b\b)  " << std::endl;
    }
  }

  unsigned numFaces = t.num_faces();
  out << std::setw(shiftwidth) << "num faces: " << numFaces << std::endl;
  if (numFaces > 0) {
    for (unsigned i=0; i<numFaces; ++i) {
      out << std::setw(shiftwidth) << t.face_topology(i) << " " << i << ": (";
      t.face_node_ordinals(i,node_ordinals);
      for (unsigned j=0, ne = t.face_topology(i).num_nodes(); j < ne; ++j) {
        out << node_ordinals[j] << ", ";
      }
      out << "\b\b)  " << std::endl;
    }
  }

  // jvo: is positive permutation according to right-hand-rule?
  unsigned numPermutations = t.num_permutations();
  unsigned numPositivePermutations = t.num_positive_permutations();
  out << std::setw(shiftwidth) << "num permutations: " << numPermutations << std::endl;
  out << std::setw(shiftwidth) << "num positive permutations: " << numPositivePermutations << std::endl;
  if (numPermutations > 0) {
    for (unsigned i=0; i<numPositivePermutations; ++i) {
      out << std::setw(shiftwidth) << i << ": (";
      t.permutation_node_ordinals(i,node_ordinals);
      for (unsigned j=0, ne = t.num_nodes(); j < ne; ++j) {
        out << node_ordinals[j] << ", ";
      }
      out << "\b\b)  " << std::endl;
    }
    out << std::setw(shiftwidth) << "num negative permutations: " << numPermutations - numPositivePermutations << std::endl;
    if (numPositivePermutations < numPermutations) {
      for (unsigned i=numPositivePermutations; i<numPermutations; ++i) {
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

bool is_solid_element(stk::topology t)
{
    return t.rank()==stk::topology::ELEM_RANK && !t.is_shell() && t.dimension()==3;
}
} //namespace stk


