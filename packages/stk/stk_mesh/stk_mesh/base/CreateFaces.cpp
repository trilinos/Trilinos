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

#include <stk_mesh/base/CreateFaces.hpp> // for create_faces

#include <stddef.h>                     // for size_t
#include <algorithm>                    // for swap, lower_bound, max, etc
#include <functional>                   // for equal_to
#include <iterator>                     // for back_insert_iterator, etc
#include <vector>                       // for vector, etc

#include <boost/array.hpp>              // for array
#include "boost/unordered/unordered_map.hpp"
#include "boost/utility/enable_if.hpp"  // for enable_if_c

#include <stk_mesh/base/BulkData.hpp>   // for BulkData, EntityLess, etc
#include <stk_mesh/base/Entity.hpp>     // for Entity, hash_value
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, get_cell_topology
#include <stk_mesh/base/Selector.hpp>   // for operator&, Selector, etc
#include <stk_mesh/base/Types.hpp>      // for EntityVector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Part.hpp"       // for Part
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

#include "stk_topology/apply_functor.tcc"  // for topology::apply_functor
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_topology/topology.tcc"    // for topology::num_nodes
#include "stk_topology/topology_type.tcc"  // for topology::topology_type

#include <stk_util/parallel/ParallelComm.hpp>  // for CommBuffer, CommAll
#include "stk_util/util/NamedPair.hpp"  // for EntityCommInfo::operator=, etc
#include <stk_mesh/base/CreateEdges.hpp>

namespace stk {
namespace mesh {

namespace {

typedef std::vector<EntityKey> EntityKeyVector;
typedef std::vector<EntityId>  EntityIdVector;

struct shared_face_type
{
  stk::topology::topology_t topology;
  EntityKeyVector           nodes;
  EntityKey                 local_key;
  EntityKey                 global_key;

  shared_face_type(stk::topology my_topology) :
    topology(my_topology.value())
  {
    nodes.resize(my_topology.num_nodes());
  }

  shared_face_type(const shared_face_type & a) :
    topology(a.topology),
    nodes(a.nodes),
    local_key(a.local_key),
    global_key(a.global_key)
  {}

  shared_face_type & operator = (const shared_face_type & a)
  {
    nodes = a.nodes;
    topology = a.topology;
    local_key = a.local_key;
    global_key = a.global_key;

    return *this;

  }
};

typedef boost::unordered_map<EntityVector,Entity> face_map_type;
typedef std::vector< shared_face_type > shared_face_map_type;

struct create_face_impl
{
  typedef void result_type;

  create_face_impl(   size_t & count_faces,
        std::vector<stk::mesh::EntityId>               & available_ids
      , face_map_type        & face_map
      , shared_face_map_type & shared_face_map
      , Bucket               & bucket
  )
  : m_count_faces(count_faces)
  , m_available_ids(available_ids)
  , m_face_map(face_map)
  , m_shared_face_map(shared_face_map)
  , m_bucket(bucket)
  {}

  template <typename Topology>
  typename boost::enable_if_c< (Topology::num_faces > 0u), void>::type
  operator()(Topology t)
  {
    typedef topology::topology_type< Topology::value> ElemTopology;

    ElemTopology elemTopology;

    BulkData & mesh = m_bucket.mesh();

    boost::array<EntityId,Topology::num_nodes> elem_node_ids;

    for (size_t ielem=0, eelem=m_bucket.size(); ielem<eelem; ++ielem) {
      Entity const *elem_nodes = m_bucket.begin_nodes(ielem);
      ThrowRequire(m_bucket.num_nodes(ielem) == Topology::num_nodes);
      for (size_t n=0; n<Topology::num_nodes; ++n) {
        elem_node_ids[n] = mesh.identifier(elem_nodes[n]);
      }

      std::vector<bool> face_exists(Topology::num_faces,false);
      {
        const int num_existing_faces = m_bucket.num_faces(ielem);
        ThrowRequire(num_existing_faces <= static_cast<int>(Topology::num_faces));

        Entity const *face_entity = m_bucket.begin_faces(ielem);
        ConnectivityOrdinal const *face_ords = m_bucket.begin_face_ordinals(ielem);
        for (int side_ordinal=0 ; side_ordinal < num_existing_faces ; ++side_ordinal) {
          face_exists[face_ords[side_ordinal]] = mesh.is_valid(face_entity[side_ordinal]);
        }
      }

      for (unsigned side_ordinal=0; side_ordinal < Topology::num_faces; ++side_ordinal) {
        if (face_exists[side_ordinal]) {
          continue;
        }
        topology faceTopology = elemTopology.face_topology(side_ordinal);

        // Use node identifier instead of node local_offset for cross-processor consistency.
        EntityIdVector side_node_ids(faceTopology.num_nodes());
        Topology::face_nodes(elem_node_ids, side_ordinal, side_node_ids.begin());
        unsigned smallest_permutation;
        EntityVector permuted_face_nodes(faceTopology.num_nodes());
        //if this is a shell OR these nodes are connected to a shell
        EntityVector side_nodes(faceTopology.num_nodes());
        for (unsigned count=0 ; count<faceTopology.num_nodes() ; ++count) {
            side_nodes[count] = mesh.get_entity(stk::topology::NODE_RANK,side_node_ids[count]);
        }
        bool is_connected_to_shell = stk::mesh::impl::do_these_nodes_have_any_shell_elements_in_common(mesh,faceTopology.num_nodes(),&side_nodes[0]);

        if (elemTopology.is_shell || is_connected_to_shell) {

            EntityIdVector element_node_id_vector(faceTopology.num_nodes());
            EntityIdVector element_node_ordinal_vector(faceTopology.num_nodes());
            EntityVector element_node_vector(faceTopology.num_nodes());
            elemTopology.face_node_ordinals(side_ordinal, &element_node_ordinal_vector[0]);
            for (unsigned count = 0; count < faceTopology.num_nodes(); ++count) {
                element_node_vector[count] = m_bucket.begin_nodes(ielem)[element_node_ordinal_vector[count]];
                element_node_id_vector[count] = mesh.identifier(element_node_vector[count]);
            }
            smallest_permutation = faceTopology.lexicographical_smallest_permutation_preserve_polarity(side_node_ids, element_node_id_vector);
            faceTopology.permutation_nodes(&element_node_vector[0], smallest_permutation, permuted_face_nodes.begin());
        }
        else {
            smallest_permutation = faceTopology.lexicographical_smallest_permutation(side_node_ids);
            EntityVector face_nodes(faceTopology.num_nodes());
            Topology::face_nodes(elem_nodes, side_ordinal, face_nodes.begin());
            faceTopology.permutation_nodes(face_nodes, smallest_permutation, permuted_face_nodes.begin());
        }
        Entity face;

        typename face_map_type::iterator iface = m_face_map.find(permuted_face_nodes);
        if (iface == m_face_map.end()) {
          ThrowRequireMsg(m_count_faces < m_available_ids.size(), "Error: face generation exhausted available identifier list. Report to sierra-help");
          EntityId face_id = m_available_ids[m_count_faces];
          m_count_faces++;

          PartVector add_parts;
          add_parts.push_back( & mesh.mesh_meta_data().get_cell_topology_root_part( get_cell_topology( faceTopology)));

          face = mesh.declare_entity( stk::topology::FACE_RANK, face_id, add_parts);
          m_face_map[permuted_face_nodes] = face;
          bool shared_face = true;
          const int num_face_nodes = faceTopology.num_nodes();
          for (int n=0; n<num_face_nodes; ++n) {
            Entity & node = permuted_face_nodes[n];
            mesh.declare_relation(face,node,n);
            shared_face = shared_face && mesh.bucket(node).shared();
          }
          if (shared_face) {
            shared_face_type sface(faceTopology);
            for (int n=0; n<num_face_nodes; ++n) {
              sface.nodes[n] = mesh.entity_key(permuted_face_nodes[n]);
            }
            const EntityKey &face_key = mesh.entity_key(face);
            sface.local_key   = face_key;
            sface.global_key  = face_key;
            m_shared_face_map.push_back( sface );
          }
        }
        else {
          face = iface->second;
        }
        mesh.declare_relation(m_bucket[ielem], face, side_ordinal);
      }
    }
  }

  template <typename Topology>
  typename boost::enable_if_c< (Topology::num_faces == 0u), void>::type
  operator()(Topology t)
  {}


  //members
  size_t                                          & m_count_faces;
  std::vector<stk::mesh::EntityId>                & m_available_ids;
  face_map_type         & m_face_map;
  shared_face_map_type  & m_shared_face_map;
  Bucket                & m_bucket;
};

} //namespace

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
  std::vector<stk::mesh::EntityId> ids_requested;

  std::vector<unsigned> localEntityCounts;
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

  shared_face_map_type shared_face_map;
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
     Bucket &b = *shells_first_element_buckets[i];
     create_face_impl functor( count_faces, ids_requested, face_map, shared_face_map, b);
     stk::topology::apply_functor< create_face_impl > apply(functor);
     apply( b.topology() );
  }
  if (connect_faces_to_edges) {
      // connect pre-existing edges to new faces
      impl::connect_faces_to_edges(mesh, element_selector, edge_map);
  }

  if (i_started) {
    bool oldOption = mesh.use_entity_ids_for_resolving_sharing();
    mesh.set_use_entity_ids_for_resolving_sharing(false);
    mesh.modification_end_for_entity_creation( stk::topology::FACE_RANK, BulkData::MOD_END_COMPRESS_AND_SORT );
    mesh.set_use_entity_ids_for_resolving_sharing(oldOption);
  }
}

}
}
