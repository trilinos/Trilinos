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

#include <stk_mesh/base/CreateFaces.hpp> // for create_faces

#include <stddef.h>                     // for size_t
#include <algorithm>                    // for swap, lower_bound, max, etc
#include <functional>                   // for equal_to
#include <iterator>                     // for back_insert_iterator, etc
#include <vector>                       // for vector, etc
#include <unordered_map>
#include <array>

#include <stk_mesh/base/BulkData.hpp>   // for BulkData, EntityLess, etc
#include <stk_mesh/base/FindPermutation.hpp>
#include <stk_mesh/base/Entity.hpp>     // for Entity, hash_value
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, get_cell_topology
#include <stk_mesh/base/Selector.hpp>   // for operator&, Selector, etc
#include <stk_mesh/base/Types.hpp>      // for EntityVector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Part.hpp"       // for Part
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/NamedPair.hpp"  // for EntityCommInfo::operator=, etc
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/baseImpl/ConnectEdgesImpl.hpp>



namespace stk {
namespace mesh {

namespace {

struct shared_face_type
{
  stk::topology::topology_t topology;
  EntityKeyVector           nodes;
  EntityKey                 local_key;
  EntityKey                 global_key;

  shared_face_type(stk::topology my_topology) :
    topology(my_topology.value()),
    nodes(my_topology.num_nodes())
  {
  }

  shared_face_type(const shared_face_type & a) :
    topology(a.topology),
    nodes(a.nodes),
    local_key(a.local_key),
    global_key(a.global_key)
  {}

  //shared_face_type & operator = (const shared_face_type & a)
  //{
  //  nodes = a.nodes;
  //  topology = a.topology;
  //  local_key = a.local_key;
  //  global_key = a.global_key;
  //  return *this;
  //}
};

typedef std::unordered_map<EntityVector, Entity, stk::mesh::impl::HashValueForEntityVector> face_map_type;
typedef std::vector< shared_face_type > shared_face_map_type;

struct create_face_impl
{
  typedef void result_type;

  create_face_impl(   size_t & count_faces,
        std::vector<stk::mesh::EntityId>               & available_ids
      , face_map_type        & face_map
      , Bucket               & bucket
      , FaceCreationBehavior faceCreationBehavior
  )
  : m_count_faces(count_faces)
  , m_available_ids(available_ids)
  , m_face_map(face_map)
  , m_bucket(bucket)
  , m_face_creation_behavior(faceCreationBehavior)
  {}

  template <typename Topology>
  void operator()(Topology)
  {
    typedef topology::topology_type< Topology::value> ElemTopology;

    ElemTopology elemTopology;

    BulkData & mesh = m_bucket.mesh();

    std::array<EntityId,Topology::num_nodes> elem_node_ids;

    for (size_t ielem=0, eelem=m_bucket.size(); ielem<eelem; ++ielem) {
      Entity const *elem_nodes = m_bucket.begin_nodes(ielem);
      STK_ThrowRequire(m_bucket.num_nodes(ielem) == Topology::num_nodes);
      for (unsigned n=0; n != Topology::num_nodes; ++n) {
        elem_node_ids[n] = mesh.identifier(elem_nodes[n]);
      }

      std::vector<bool> face_exists(Topology::num_faces,false);
      {
        const int num_existing_faces = m_bucket.num_faces(ielem);
        STK_ThrowRequire(num_existing_faces <= static_cast<int>(Topology::num_faces));

        Entity const *face_entity = m_bucket.begin_faces(ielem);
        ConnectivityOrdinal const *face_ords = m_bucket.begin_face_ordinals(ielem);
        for (int side_ordinal=0 ; side_ordinal < num_existing_faces ; ++side_ordinal) {
          face_exists[face_ords[side_ordinal]] = mesh.is_valid(face_entity[side_ordinal]);
        }
      }

      for (unsigned side_ordinal=0; side_ordinal != Topology::num_faces; ++side_ordinal) {

          if (!face_exists[side_ordinal]) {
              if (m_face_creation_behavior == FaceCreationBehavior::CREATE_FACES_FACE_CREATION_CLASSIC) {
                  topology faceTopology = elemTopology.face_topology(side_ordinal);
                  stk::mesh::Entity element = m_bucket[ielem];
                  EntityVector permuted_face_nodes(faceTopology.num_nodes());
                  stk::mesh::impl::find_side_nodes(mesh, element, side_ordinal, permuted_face_nodes);
                  Entity face;

                  typename face_map_type::iterator iface = m_face_map.find(permuted_face_nodes);
                  if (iface == m_face_map.end()) {
                      STK_ThrowRequireMsg(m_count_faces < m_available_ids.size(), "Error: face generation exhausted available identifier list. Report to sierra-help");
                      EntityId face_id = m_available_ids[m_count_faces];
                      m_count_faces++;

                      PartVector add_parts;
                      add_parts.push_back( & mesh.mesh_meta_data().get_topology_root_part(faceTopology) );

                      face = mesh.declare_solo_side(face_id, add_parts);
                      m_face_map[permuted_face_nodes] = face;

                      const int num_face_nodes = faceTopology.num_nodes();
                      for (int n=0; n<num_face_nodes; ++n) {
                          Entity & node = permuted_face_nodes[n];
                          mesh.declare_relation(face,node,n);
                      }

                      Permutation permut = stk::mesh::find_permutation(mesh, elemTopology, elem_nodes,
                                                                                       faceTopology, &permuted_face_nodes[0], side_ordinal);
                      mesh.declare_relation(m_bucket[ielem], face, side_ordinal, permut);
                  }
                  else {
                      face = iface->second;
                      Permutation permut = stk::mesh::find_permutation(mesh, elemTopology, elem_nodes,
                                                                 faceTopology, &permuted_face_nodes[0], side_ordinal);
                      STK_ThrowRequireMsg(permut != INVALID_PERMUTATION, "CreateFaces:  could not find valid permutation to connect face to element");
                      mesh.declare_relation(m_bucket[ielem], face, side_ordinal, permut);
                  }
              }
              else { //
                  topology faceTopology = elemTopology.face_topology(side_ordinal);
                  stk::mesh::Entity elem = m_bucket[ielem];
                  stk::mesh::Part & face_topology_part = mesh.mesh_meta_data().get_topology_root_part(faceTopology);
                  stk::mesh::Entity new_face = stk::mesh::impl::get_or_create_face_at_element_side(
                          mesh, elem, side_ordinal, m_available_ids[m_count_faces], stk::mesh::PartVector(1,&face_topology_part));
                  if (mesh.identifier(new_face) == m_available_ids[m_count_faces]) {
                      stk::mesh::impl::connect_face_to_other_elements(mesh, new_face, elem,side_ordinal);
                      m_count_faces++;
                  }
              }
          }
      }
    }
  }

  //members
  size_t                                          & m_count_faces;
  std::vector<stk::mesh::EntityId>                & m_available_ids;
  face_map_type         & m_face_map;
  Bucket                & m_bucket;
  FaceCreationBehavior   m_face_creation_behavior;
};

} //namespace

namespace experimental {
void create_faces( BulkData & mesh )
{
    stk::mesh::create_all_sides(mesh, mesh.mesh_meta_data().universal_part(), stk::mesh::PartVector(), false);
}

void create_faces( BulkData & mesh, const Selector & element_selector )
{
    stk::mesh::create_all_sides(mesh, element_selector, stk::mesh::PartVector(), false);
}

void create_faces( BulkData & mesh, const Selector & element_selector, Part *part_to_insert_new_faces)
{
    stk::mesh::PartVector parts = {part_to_insert_new_faces};
    stk::mesh::create_all_sides(mesh, element_selector, parts, false);
}

void create_faces( BulkData & mesh, bool connect_faces_to_edges)
{
    stk::mesh::create_all_sides(mesh, mesh.mesh_meta_data().universal_part(), stk::mesh::PartVector(), connect_faces_to_edges);
}

void create_faces( BulkData & mesh, const Selector & element_selector, bool connect_faces_to_edges)
{
    stk::mesh::create_all_sides(mesh, element_selector, stk::mesh::PartVector(), connect_faces_to_edges);
}
}


void internal_create_faces( BulkData & mesh, const Selector & element_selector, bool connect_faces_to_edges, FaceCreationBehavior faceCreationBehavior);

void create_faces( BulkData & mesh )
{
  create_faces(mesh, mesh.mesh_meta_data().universal_part());
}

void create_faces( BulkData & mesh, const Selector & element_selector)
{
    create_faces(mesh, element_selector, false);
}

void create_faces( BulkData & mesh, bool connect_faces_to_edges)
{
  create_faces(mesh, mesh.mesh_meta_data().universal_part(), connect_faces_to_edges);
}

void create_faces( BulkData & mesh, const Selector & element_selector, bool connect_faces_to_edges)
{
    internal_create_faces(mesh, element_selector, connect_faces_to_edges, FaceCreationBehavior::CREATE_FACES_FACE_CREATION_CLASSIC);
}


void internal_create_faces( BulkData & mesh, const Selector & element_selector, bool connect_faces_to_edges, FaceCreationBehavior faceCreationBehavior)
{
  std::vector<stk::mesh::EntityId> ids_requested;

  std::vector<size_t> localEntityCounts;
  stk::mesh::count_entities(element_selector, mesh, localEntityCounts);
  unsigned guessMultiplier = 6;
  unsigned numRequested = localEntityCounts[stk::topology::ELEMENT_RANK] * guessMultiplier;
  mesh.generate_new_ids(stk::topology::FACE_RANK, numRequested, ids_requested);
  size_t count_faces = 0;

  impl::edge_map_type edge_map;
  if (connect_faces_to_edges) {
      //populate the edge_map with existing edges

      BucketVector const & edge_buckets = mesh.buckets(stk::topology::EDGE_RANK);

      for (size_t i=0, ie=edge_buckets.size(); i<ie; ++i) {
          Bucket &b = *edge_buckets[i];

          const unsigned num_nodes = b.topology().num_nodes();
          EntityVector edge_nodes(num_nodes);

          for (size_t j=0, je=b.size(); j<je; ++j) {
              Entity edge = b[j];
              Entity const *nodes_rel = b.begin_nodes(j);

              for (unsigned n=0; n<num_nodes; ++n) {
                  edge_nodes[n] = nodes_rel[n];
              }

              edge_map[edge_nodes] = edge;
          }
      }
  }

  bool i_started = mesh.modification_begin();

  face_map_type        face_map;

  //populate the face_map with existing faces
  BucketVector const & face_buckets = mesh.buckets(stk::topology::FACE_RANK);

  for (size_t i=0, ie=face_buckets.size(); i<ie; ++i) {
    Bucket &b = *face_buckets[i];

    const unsigned num_nodes = b.topology().num_nodes();
    EntityVector face_nodes(num_nodes);

    for (size_t j=0, je=b.size(); j<je; ++j) {
      Entity face = b[j];
      Entity const *nodes_rel = b.begin_nodes(j);

      for (unsigned n=0; n<num_nodes; ++n) {
        face_nodes[n] = nodes_rel[n];
      }

      face_map[face_nodes] = face;
    }
  }

  // create faces and connect them to elements
  BucketVector const& element_buckets =
      mesh.get_buckets(stk::topology::ELEMENT_RANK, element_selector & mesh.mesh_meta_data().locally_owned_part());

  BucketVector shells_first_element_buckets;
  shells_first_element_buckets.reserve(element_buckets.size());
  //create the faces for the elements in each bucket
  for (size_t i=0, e=element_buckets.size(); i<e; ++i) {
      Bucket *b = element_buckets[i];
      if (b->topology().is_shell()) {
          shells_first_element_buckets.push_back(b);
      }
  }
  for (size_t i=0, e=element_buckets.size(); i<e; ++i) {
      Bucket *b = element_buckets[i];
      if ( ! b->topology().is_shell()) {
          shells_first_element_buckets.push_back(b);
      }
  }
  for (size_t i=0, e=shells_first_element_buckets.size(); i<e; ++i) {
     Bucket &this_bucket = *shells_first_element_buckets[i];
     create_face_impl functor( count_faces, ids_requested, face_map, this_bucket, faceCreationBehavior);
     stk::topology::apply_host_functor< create_face_impl > apply(functor);
     apply( this_bucket.topology() );
  }
  if (connect_faces_to_edges) {
      // connect pre-existing edges to new faces
      impl::connect_faces_to_edges(mesh, element_selector, edge_map);
  }

  if (i_started) {
    bool oldOption = mesh.use_entity_ids_for_resolving_sharing();
    mesh.set_use_entity_ids_for_resolving_sharing(false);
    std::vector<EntityRank> entity_rank_vector = {stk::topology::FACE_RANK};
    mesh.modification_end_for_entity_creation( entity_rank_vector );
    mesh.set_use_entity_ids_for_resolving_sharing(oldOption);
  }
}

}
}
