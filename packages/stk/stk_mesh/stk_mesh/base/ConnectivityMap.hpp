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

#ifndef STK_STK_CONNECTIVITY_MAP_HPP
#define STK_STK_CONNECTIVITY_MAP_HPP

#include <stk_mesh/base/Types.hpp>
#include <array>

namespace stk  {
namespace mesh {

// TODO - enforce invariants (such as only makes sense to have fixed connectivity to an element from other elements)

struct ConnectivityMap
{
  typedef std::array<std::array<ConnectivityType, 4>, 4> map_type;

  static ConnectivityMap const& classic_stk_mesh()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), dynamic(),  dynamic(), dynamic()}},
      /*edge*/  {{ fixed()  , invalid(),  dynamic(), dynamic()}},
      /*face*/  {{ fixed()  , dynamic(),  invalid(), dynamic()}},
      /*elem*/  {{ fixed()  , dynamic(),  dynamic(), invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& classic_stk_mesh_2d()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), dynamic(),  invalid(), dynamic()}},
      /*edge*/  {{ fixed()  , invalid(),  invalid(), dynamic()}},
      /*face*/  {{ invalid(), invalid(),  invalid(), invalid()}},
      /*elem*/  {{ fixed()  , dynamic(),  invalid(), invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& default_map()
  {
    return classic_stk_mesh();
  }

  static ConnectivityMap const& default_map_2d()
  {
    return classic_stk_mesh_2d();
  }

  static ConnectivityMap const& minimal_upward_connectivity_map()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), invalid(),  invalid(), dynamic()}},
      /*edge*/  {{ fixed()  , invalid(),  invalid(), invalid()}},
      /*face*/  {{ fixed()  , dynamic(),  invalid(), invalid()}},
      /*elem*/  {{ fixed()  , dynamic(),  dynamic(), invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& fixed_edges_map()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), dynamic(),  dynamic(), dynamic()}},
      /*edge*/  {{ fixed()  , invalid(),  dynamic(), dynamic()}},
      /*face*/  {{ fixed()  , fixed()  ,  invalid(), dynamic()}},
      /*elem*/  {{ fixed()  , fixed()  ,  dynamic(), invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& fixed_edges_map_2d()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), dynamic(),  invalid(), dynamic()}},
      /*edge*/  {{ fixed()  , invalid(),  invalid(), dynamic()}},
      /*face*/  {{ invalid(), invalid(),  invalid(), invalid()}},
      /*elem*/  {{ fixed()  , fixed()  ,  invalid(), invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& fixed_downward_map()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), dynamic(),  dynamic(), dynamic()}},
      /*edge*/  {{ fixed()  , invalid(),  dynamic(), dynamic()}},
      /*face*/  {{ fixed()  , fixed()  ,  invalid(), dynamic()}},
      /*elem*/  {{ fixed()  , fixed()  ,  fixed()  , invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& element_node()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), invalid(),  invalid(), dynamic()}},
      /*edge*/  {{ invalid(), invalid(),  invalid(), invalid()}},
      /*face*/  {{ invalid(), invalid(),  invalid(), invalid()}},
      /*elem*/  {{ fixed()  , invalid(),  invalid(), invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& none()
  {
    static const map_type map =
    {{
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& none_2d()
  {
    static const map_type map =
    {{
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& all_dynamic()
  {
    static const map_type map =
    {{
      {{dynamic(),dynamic(),dynamic(),dynamic()}},
      {{dynamic(),dynamic(),dynamic(),dynamic()}},
      {{dynamic(),dynamic(),dynamic(),dynamic()}},
      {{dynamic(),dynamic(),dynamic(),dynamic()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& all_dynamic_2d()
  {
    static const map_type map =
    {{
      {{dynamic(),dynamic(),invalid(),dynamic()}},
      {{dynamic(),dynamic(),invalid(),dynamic()}},
      {{invalid(),invalid(),invalid(),invalid()}},
      {{dynamic(),dynamic(),invalid(),dynamic()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityType invalid() { return INVALID_CONNECTIVITY_TYPE; }
  static ConnectivityType fixed()   { return FIXED_CONNECTIVITY;   }
  static ConnectivityType dynamic() { return DYNAMIC_CONNECTIVITY; }

  ConnectivityType & operator()(EntityRank from, EntityRank to)
  { return m_map[from][to]; }

  const ConnectivityType & operator()(EntityRank from, EntityRank to) const
  { return m_map[from][to]; }

  bool valid(EntityRank from, EntityRank to) const
  {
    // You can always have dynamic connectivity to "beyond-element" entities
    return (from > stk::topology::ELEMENT_RANK || to > stk::topology::ELEMENT_RANK) ? true : m_map[from][to] != INVALID_CONNECTIVITY_TYPE;
  }

  map_type m_map;
};

inline std::ostream & operator<<(std::ostream &out, const ConnectivityMap & map)
{
   out << "From\\To:  NODE     EDGE     FACE     ELEM\n";
   out << "NODE   :  "; for (int i=0 ; i<4 ; ++i) { out << map.m_map[0][i] << "  "; } out << "\n";
   out << "EDGE   :  "; for (int i=0 ; i<4 ; ++i) { out << map.m_map[1][i] << "  "; } out << "\n";
   out << "FACE   :  "; for (int i=0 ; i<4 ; ++i) { out << map.m_map[2][i] << "  "; } out << "\n";
   out << "ELEM   :  "; for (int i=0 ; i<4 ; ++i) { out << map.m_map[3][i] << "  "; } out << "\n";
   return out;
}

} //namespace mesh
} //namespace stk

#endif //STK_STK_CONNECTIVITY_MAP_HPP
