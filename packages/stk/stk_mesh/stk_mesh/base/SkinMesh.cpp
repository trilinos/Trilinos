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

#include <stk_mesh/base/SkinMesh.hpp>
#include <stddef.h>                     // for size_t
#include <algorithm>                    // for sort, set_intersection
#include <iostream>                     // for operator<<, cerr, ostream
#include <iterator>                     // for back_insert_iterator, etc
#include <set>                          // for _Rb_tree_const_iterator, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/Selector.hpp>   // for Selector, operator&
#include <utility>                      // for pair, make_pair
#include <vector>                       // for vector
#include <boost/tuple/tuple.hpp>
#include "boost/tuple/detail/tuple_basic.hpp"  // for tie, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, get_cell_topology
#include "stk_mesh/base/Types.hpp"      // for EntityVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_topology/topology.tcc"    // for topology::num_nodes, etc
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowErrorMsgIf
namespace stk { namespace mesh { class Part; } }


namespace stk { namespace mesh {

namespace {

void get_common_elements( BulkData const& mesh, EntityVector const& nodes, EntityVector & elements)
{
  elements.clear();

  if (nodes.size() >= 1u) {
    elements.insert(elements.end(), mesh.begin_elements(nodes[0]), mesh.end_elements(nodes[0]));
    std::sort(elements.begin(), elements.end());

    EntityVector tmp_out, tmp ;
    tmp_out.reserve(elements.size());
    tmp.reserve(elements.size());

    for (size_t i=1, ie=nodes.size(); i<ie; ++i) {
      tmp.insert(tmp.end(), mesh.begin_elements(nodes[i]), mesh.end_elements(nodes[i]));
      std::sort(tmp.begin(),tmp.end());

      std::set_intersection(  elements.begin(),elements.end()
                            , tmp.begin(), tmp.end()
                            , std::back_inserter(tmp_out)
                           );

      elements.swap(tmp_out);
      tmp_out.clear();
      tmp.clear();
    }
  }
}

typedef std::set< std::pair<Entity, unsigned> > Boundary;

size_t skin_mesh_find_elements_with_external_sides(BulkData & mesh,
                                                   const BucketVector& element_buckets,
                                                   Boundary& boundary)
{
  const EntityRank side_rank = mesh.mesh_meta_data().side_rank();
  size_t num_sides_to_create = 0;
  for (size_t i=0, ie=element_buckets.size(); i<ie; ++i) {
    const Bucket & b = *element_buckets[i];
    stk::topology element_topology = b.topology();
    const unsigned num_sides = element_topology.num_sides();

    // skip over elements which do not have sides
    if (num_sides == 0u) continue;

    // TODO: skip shells and particles for now
    if (element_topology.is_shell() || (element_topology == stk::topology::PARTICLE) ) {
      std::cerr << "Skinning shells and particlesis currently not supported!";
      //ThrowErrorMsgIf(element_topology.is_shell(), "Skinnig meshes with shells is currently unsupported");
      continue;
    }

    for (size_t j=0, je=b.size(); j<je; ++j) {
      Entity elem = b[j];

      Entity const * const elem_nodes = mesh.begin_nodes(elem);

      for (unsigned k=0; k<num_sides; ++k) {

        stk::topology side_topology = element_topology.side_topology(k);
        // get the side nodes
        EntityVector side_nodes(side_topology.num_nodes());
        element_topology.side_nodes(elem_nodes,k, side_nodes.begin());

        // find elements that also use the side nodes
        EntityVector common_elements;
        get_common_elements(mesh, side_nodes, common_elements);

        bool found_adjacent_element = false;
        for (size_t l=0, le=common_elements.size(); !found_adjacent_element && l<le; ++l) {
          Entity potential_element = common_elements[l];
          stk::topology potential_element_topology = mesh.bucket(potential_element).topology();

          if (potential_element == elem) continue;              // skip self adjacency
          if (potential_element_topology.is_shell()) continue;  // skip shells

          const unsigned potential_num_sides = potential_element_topology.num_sides();
          Entity const * const potential_elem_nodes = mesh.begin_nodes(potential_element);

          for (unsigned m=0; m < potential_num_sides; ++m) {
            stk::topology potential_side_topology = potential_element_topology.side_topology(m);

            // the side topologies are different -- not a match
            if (side_topology != potential_side_topology) continue;

            EntityVector potential_side_nodes(potential_side_topology.num_nodes());
            potential_element_topology.side_nodes(potential_elem_nodes, m, potential_side_nodes.begin());

            bool equivalent = false;
            unsigned permutation_id = 0;

            boost::tie(equivalent,permutation_id) = side_topology.equivalent(side_nodes, potential_side_nodes);

            // the sides are not a match
            if (equivalent == false) continue;
            // if the permutation_id is to a positive permutation
            // the sides are not opposing, i.e. the elements are superimposed on each other
            ThrowErrorMsgIf( (permutation_id < side_topology.num_positive_permutations())
                , "Error: Non shell elements cannot be superimposed");

            found_adjacent_element = true;
            break;
          }
        }
        if (!found_adjacent_element) {
          ConnectivityOrdinal const* existing_side_ordinals = mesh.begin_ordinals(elem,side_rank);
          unsigned num_exisitng_sides = mesh.num_connectivity(elem,side_rank);

          bool side_exists = false;
          for ( unsigned s=0; !side_exists && s<num_exisitng_sides; ++s) {
            if (existing_side_ordinals[s] == k)
              side_exists=true;
          }
          if (!side_exists)
            ++num_sides_to_create;

          boundary.insert(std::make_pair(elem,k));;
        }
      }
    }
  }
  return num_sides_to_create;
}

void skin_mesh_attach_new_sides_to_connected_entities(BulkData & mesh,
                                                      Boundary& boundary,
                                                      EntityVector const& sides,
                                                      PartVector const& skin_parts)
{
  const unsigned spatial_dimension = mesh.mesh_meta_data().spatial_dimension();
  const bool consider_negative_permutations = spatial_dimension < 3u;
  const EntityRank side_rank = mesh.mesh_meta_data().side_rank();

  size_t side_index = 0;
  for (Boundary::const_iterator itr = boundary.begin(), end_itr = boundary.end();
       itr != end_itr;
       ++itr
      )
  {
    Entity elem;
    unsigned side_ordinal=0;

    boost::tie(elem,side_ordinal) = *itr;

    stk::topology element_topology = mesh.bucket(elem).topology();
    stk::topology side_topology = element_topology.side_topology(side_ordinal);

    Entity side = Entity();

    ConnectivityOrdinal const* existing_side_ordinals = mesh.begin_ordinals(elem,side_rank);
    Entity const* existing_sides = mesh.begin(elem,side_rank);
    unsigned num_exisitng_sides = mesh.num_connectivity(elem,side_rank);

    bool side_exists = false;
    for ( unsigned s=0; !side_exists && s<num_exisitng_sides; ++s) {
      if (existing_side_ordinals[s] == side_ordinal) {
        side_exists=true;
        side = existing_sides[s];
      }
    }

    if (!side_exists) {
      side = sides[side_index++];

      Entity const * const elem_nodes = mesh.begin_nodes(elem);

      EntityVector side_nodes(side_topology.num_nodes());
      element_topology.side_nodes(elem_nodes, side_ordinal, side_nodes.begin());

      unsigned perm_id = side_topology.lexicographical_smallest_permutation(side_nodes, consider_negative_permutations);

      EntityVector ordered_side_nodes(side_nodes.size());
      side_topology.permutation_nodes(side_nodes, perm_id, ordered_side_nodes.begin());

      // attach nodes to side
      for (size_t i=0, ie=ordered_side_nodes.size(); i<ie; ++i) {
        mesh.declare_relation( side, ordered_side_nodes[i], static_cast<RelationIdentifier>(i));
      }

      // attach side to element
      mesh.declare_relation( elem, side, static_cast<RelationIdentifier>(side_ordinal), static_cast<Permutation>(perm_id));
    }

    PartVector add_parts(skin_parts);
    add_parts.push_back( & mesh.mesh_meta_data().get_cell_topology_root_part( get_cell_topology( side_topology )));

    mesh.change_entity_parts(side,add_parts);
  }
}

} //end un-named namespace

void skin_mesh( BulkData & mesh, PartVector const& skin_parts)
{
  skin_mesh(mesh, mesh.mesh_meta_data().universal_part(), skin_parts);
}

void skin_mesh( BulkData & mesh, Selector const& element_selector, PartVector const& skin_parts)
{
  ThrowErrorMsgIf( mesh.synchronized_state() == BulkData::MODIFIABLE, "mesh is not SYNCHRONIZED" );

  const Part & locally_owned = mesh.mesh_meta_data().locally_owned_part();
  const EntityRank side_rank = mesh.mesh_meta_data().side_rank();
  const BucketVector& element_buckets = mesh.get_buckets(stk::topology::ELEMENT_RANK, element_selector & locally_owned);

  Boundary boundary;

  // find elements with external sides
  size_t num_sides_to_create = skin_mesh_find_elements_with_external_sides(mesh,
                                                                           element_buckets,
                                                                           boundary);
  // create the sides
  std::vector<size_t> requests(mesh.mesh_meta_data().entity_rank_count(), 0);
  requests[side_rank] = num_sides_to_create;

  mesh.modification_begin();

  EntityVector sides;
  mesh.generate_new_entities(requests, sides);

  // attach new sides to connected entities
  skin_mesh_attach_new_sides_to_connected_entities(mesh, boundary, sides, skin_parts);

  mesh.modification_end();
}

}} // namespace stk::mesh
