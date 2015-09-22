// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#ifndef hex_refine_info_hpp
#define hex_refine_info_hpp

namespace stk {
namespace performance_tests {

class HexRefineInfo
{
public:
  typedef unsigned HRF_EntityId;

  const unsigned m_level;
  const unsigned                m_nx;
  const unsigned                m_ny;
  const unsigned                m_nz;
  const unsigned m_one_based_offset;

  // nx, ny, nz are the *original* sizes of the mesh
  HexRefineInfo(unsigned level, unsigned nx, unsigned ny, unsigned nz, unsigned one_based_offset=1)
    : m_level(level), m_nx(nx), m_ny(ny), m_nz(nz), m_one_based_offset(one_based_offset)
  {

  }

#define POW3(x) ((x)*(x)*(x))

  HRF_EntityId elem_id_offset(unsigned level) const
  {
    if (level == 0)
      return 0;
    else
      return elem_id_offset(level-1) +  m_nx*m_ny*m_nz*POW3(1u << (level-1));
  }

  HRF_EntityId elem_id_offset() const
  {
    return elem_id_offset(m_level);
  }

  HRF_EntityId node_id_offset(unsigned level) const
  {
    if (level == 0)
      return 0;
    else
      return node_id_offset(level-1) +  (m_nx+1)*(m_ny+1)*(m_nz+1)*POW3(1u << (level-1));
  }

  HRF_EntityId node_id_offset() const
  {
    return node_id_offset(m_level);
  }

  unsigned lnx() const { return m_nx*(1u<<m_level); }
  unsigned lny() const { return m_ny*(1u<<m_level); }
  unsigned lnz() const { return m_nz*(1u<<m_level); }
  unsigned num_nodes() const { return (1+lnx())*(1+lny())*(1+lnz()); }
  unsigned num_elems() const { return lnx()*lny()*lnz(); }
  /**
   * Thinking in terms of a 3D grid of nodes, get the id of the node in
   * the (x, y, z) position.
   */
  HRF_EntityId node_id( unsigned x , unsigned y , unsigned z ) const  {
    return node_id_offset() + m_one_based_offset + x + ( lnx() + 1 ) * ( y + ( lny() + 1 ) * z );
  }

  /**
   * Thinking in terms of a 3D grid of elements, get the id of the
   * element in the (x, y, z) position.
   */
  HRF_EntityId elem_id( unsigned x , unsigned y , unsigned z ) const  {
    return elem_id_offset() + m_one_based_offset + x + lnx() * ( y + lny() * z );
  }

  /**
   * Thinking in terms of a 3D grid of nodes, compute the (x, y, z) position
   * of a node given it's id.
   */
  void node_x_y_z( HRF_EntityId entity_id, unsigned &x , unsigned &y , unsigned &z ) const
  {
    entity_id -= m_one_based_offset;
    entity_id -= node_id_offset();

    x = entity_id % (lnx()+1);
    entity_id /= (lnx()+1);

    y = entity_id % (lny()+1);
    entity_id /= (lny()+1);

    z = entity_id;
  }

  /**
   * Thinking in terms of a 3D grid of elements, compute the (x, y, z) position
   * of an element given it's id.
   */
  void elem_x_y_z( HRF_EntityId entity_id, unsigned &x , unsigned &y , unsigned &z ) const
  {
    entity_id -= m_one_based_offset;
    entity_id -= elem_id_offset();

    x = entity_id % lnx();
    entity_id /= lnx();

    y = entity_id % lny();
    entity_id /= lny();

    z = entity_id;
  }


};

}
}

#endif
