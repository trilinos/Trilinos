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

#include <cassert>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stddef.h>                     // for size_t
#include <algorithm>                    // for sort, set_intersection
#include <iostream>                     // for operator<<, cerr, ostream
#include <iterator>                     // for back_insert_iterator, etc
#include <set>                          // for _Rb_tree_const_iterator, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/FindPermutation.hpp>
#include <stk_mesh/base/Selector.hpp>   // for Selector, operator&
#include <utility>                      // for pair, make_pair
#include <vector>                       // for vector
#include <tuple>
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, get_cell_topology
#include "stk_mesh/base/Types.hpp"      // for EntityVector, etc
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowErrorMsgIf
namespace stk { namespace mesh { class Part; } }


namespace stk { namespace mesh {

namespace {

typedef std::set< std::pair<Entity, unsigned> > Boundary;

bool isTopologySupportedForSkinning(stk::topology topology)
{
    if (topology.is_shell())
    {
        return false;
    }
    else if ( topology.dimension() < 2 )
    {
        return false;
    }
    else if ( topology.num_sides() == 0 )
    {
        return false;
    }
    return true;
}
size_t skin_mesh_find_elements_with_external_sides(BulkData & mesh,
                                                   const BucketVector& element_buckets,
                                                   Boundary& boundary,
                                                   const Selector * secondary_selector)
{
  const EntityRank side_rank = mesh.mesh_meta_data().side_rank();
  EntityVector side_nodes;
  EntityVector potential_side_nodes;
  EntityVector common_elements;
  size_t num_sides_to_create = 0;
  for (size_t i=0, ie=element_buckets.size(); i<ie; ++i) {
    const Bucket & b = *element_buckets[i];
    stk::topology element_topology = b.topology();
    if (!isTopologySupportedForSkinning(element_topology))
    {
        static std::map<stk::topology, bool> warningReported;
        if ( warningReported.find(element_topology) == warningReported.end() )
        {
            std::cerr << "Skinning " << element_topology << " is currently not supported." << std::endl;
            warningReported[element_topology] = true;
        }
        continue;
    }
    const unsigned num_sides = element_topology.num_sides();

    for (size_t j=0, je=b.size(); j<je; ++j) {
      Entity elem = b[j];

      Entity const * const elem_nodes = mesh.begin_nodes(elem);

      for (unsigned k=0; k<num_sides; ++k) {

        stk::topology side_topology = element_topology.side_topology(k);
        side_nodes.resize(side_topology.num_nodes());
        element_topology.side_nodes(elem_nodes,k, side_nodes.data());

        impl::find_entities_these_nodes_have_in_common(mesh, stk::topology::ELEM_RANK, side_nodes.size(), side_nodes.data(), common_elements);

        bool found_adjacent_element = false;
        for (size_t l=0, le=common_elements.size(); !found_adjacent_element && l<le; ++l) {
          Entity potential_element = common_elements[l];
          if (potential_element == elem) continue;              // skip self adjacency
          stk::topology potential_element_topology = mesh.bucket(potential_element).topology();

          if (potential_element_topology.is_shell()) continue;  // skip shells
          if (secondary_selector && !((*secondary_selector)(mesh.bucket(potential_element))))
            continue;  // skip elements not in the secondary selector

          const unsigned potential_num_sides = potential_element_topology.num_sides();
          Entity const * const potential_elem_nodes = mesh.begin_nodes(potential_element);

          for (unsigned m=0; m < potential_num_sides; ++m) {
            stk::topology potential_side_topology = potential_element_topology.side_topology(m);

            // the side topologies are different -- not a match
            if (side_topology != potential_side_topology) continue;

            potential_side_nodes.resize(potential_side_topology.num_nodes());
            potential_element_topology.side_nodes(potential_elem_nodes, m, potential_side_nodes.data());

            stk::EquivalentPermutation result = stk::mesh::side_equivalent(mesh, elem, k, potential_side_nodes.data());

            // the sides are not a match
            if (result.is_equivalent == false) continue;
            // if the permutation_id is to a positive permutation
            // the sides are not opposing, i.e. the elements are superimposed on each other

            if ( result.permutation_number < side_topology.num_positive_permutations() )
            {
                stk::mesh::PartVector user_parts_elem1;
                stk::mesh::PartVector user_parts_elem2;

                const stk::mesh::PartVector &parts1 = mesh.bucket(elem).supersets();
                const stk::mesh::PartVector &parts2 = mesh.bucket(potential_element).supersets();
                for (size_t mm=0;mm<parts1.size();mm++)
                {
                    if ( !stk::mesh::is_auto_declared_part(*parts1[mm]) )
                    {
                        user_parts_elem1.push_back(parts1[mm]);
                    }
                }
                for (size_t mm=0;mm<parts2.size();mm++)
                {
                    if ( !stk::mesh::is_auto_declared_part(*parts2[mm]) )
                    {
                        user_parts_elem2.push_back(parts2[mm]);
                    }
                } 
                std::ostringstream os;
                os << "*** Warning: element " << mesh.identifier(elem) << " ( ";
                for (size_t mm=0;mm<user_parts_elem1.size();mm++)
                {
                    os << user_parts_elem1[mm]->name() << " ";
                } 
                os << ") and element " << mesh.identifier(potential_element) << " ( ";
                for (size_t mm=0;mm<user_parts_elem2.size();mm++)
                {
                    os << user_parts_elem2[mm]->name() << " ";
                } 
                os << ") have overlapping volumes. ***" << std::endl;
                std::cerr << os.str();
                continue;
            }
            
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

          boundary.insert(std::make_pair(elem,k));
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

    std::tie(elem,side_ordinal) = *itr;

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
      element_topology.side_nodes(elem_nodes, side_ordinal, side_nodes.data());

      unsigned lexmin_perm_id = side_topology.lexicographical_smallest_permutation(side_nodes.data(), consider_negative_permutations);

      EntityVector ordered_side_nodes(side_nodes.size());
      side_topology.permutation_nodes(side_nodes.data(), lexmin_perm_id, ordered_side_nodes.data());

      // attach nodes to side
      for (size_t i=0, ie=ordered_side_nodes.size(); i<ie; ++i) {
        mesh.declare_relation( side, ordered_side_nodes[i], static_cast<RelationIdentifier>(i));
      }

      // attach side to element
      Permutation permut =
              stk::mesh::find_permutation(mesh, element_topology, elem_nodes,
                                    side_topology, ordered_side_nodes.data(), side_ordinal);
      STK_ThrowRequireMsg(permut != INVALID_PERMUTATION, ":  skin_mesh_attach_new_sides_to_connected_entities could not find valid permutation to connect face to element");

      mesh.declare_relation(elem, side, static_cast<RelationIdentifier>(side_ordinal), permut);
    // mesh.declare_relation( elem, side, static_cast<RelationIdentifier>(side_ordinal), static_cast<Permutation>(perm_id));

    }

    bool need_to_add_side_to_skin_part_reason1 = side_exists && !mesh.bucket(side).shared();
    bool need_to_add_side_to_skin_part_reason2 = !side_exists;

    if (need_to_add_side_to_skin_part_reason1 || need_to_add_side_to_skin_part_reason2)
    {
        bool only_owner_can_change_parts = mesh.bucket(side).owned();
        if (only_owner_can_change_parts)
        {
          PartVector add_parts(skin_parts);
          add_parts.push_back( & mesh.mesh_meta_data().get_topology_root_part(side_topology) );
          mesh.change_entity_parts(side,add_parts);
        }
    }
  }
}

} //end un-named namespace

void skin_mesh( BulkData & mesh, PartVector const& skin_parts)
{
  skin_mesh(mesh, mesh.mesh_meta_data().universal_part(), skin_parts);
}

void skin_mesh( BulkData & mesh, Selector const& element_selector, PartVector const& skin_parts,
                const Selector * secondary_selector)
{
  STK_ThrowErrorMsgIf( mesh.in_modifiable_state(), "mesh is not SYNCHRONIZED" );

  const Part & locally_owned = mesh.mesh_meta_data().locally_owned_part();
  const EntityRank side_rank = mesh.mesh_meta_data().side_rank();
  const BucketVector& element_buckets = mesh.get_buckets(stk::topology::ELEMENT_RANK, element_selector & locally_owned);

  Boundary boundary;
  // find elements with external sides
  size_t num_sides_to_create = skin_mesh_find_elements_with_external_sides(mesh,
                                                                           element_buckets,
                                                                           boundary, secondary_selector);
  // create the sides
  std::vector<size_t> requests(mesh.mesh_meta_data().entity_rank_count(), 0);
  requests[side_rank] = num_sides_to_create;

  mesh.modification_begin();

  EntityVector sides;
  mesh.generate_new_entities(requests, sides);

  // attach new sides to connected entities
  skin_mesh_attach_new_sides_to_connected_entities(mesh, boundary, sides, skin_parts);

  mesh.internal_modification_end_for_skin_mesh(mesh.mesh_meta_data().side_rank(), impl::MeshModification::MOD_END_SORT, element_selector, secondary_selector);
}

}} // namespace stk::mesh
