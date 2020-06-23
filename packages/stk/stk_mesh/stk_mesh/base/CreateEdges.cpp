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

#include "stk_mesh/base/CreateEdges.hpp"

#include <stddef.h>                     // for size_t
#include <algorithm>                    // for swap, lower_bound, max, etc
#include <functional>                   // for equal_to
#include <iterator>                     // for back_insert_iterator, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, EntityLess, etc
#include <stk_mesh/base/Entity.hpp>     // for Entity, hash_value
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, get_cell_topology
#include <stk_mesh/base/Selector.hpp>   // for operator&, Selector, etc
#include <stk_mesh/base/Types.hpp>      // for EntityVector, etc
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <vector>                       // for vector, etc
#include <unordered_map>
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_topology/apply_functor.hpp"  // for topology::apply_host_functor
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_topology/topology_utils.hpp"    // for topology::num_nodes
#include "stk_topology/topology_type.hpp"  // for topology::topology_type
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert
#include "stk_util/util/NamedPair.hpp"  // for EntityCommInfo::operator=, etc

namespace stk {
namespace mesh {

namespace impl {
  typedef std::unordered_map<EntityVector,Entity,stk::mesh::impl::HashValueForEntityVector> edge_map_type;
} //impl namespace

namespace {

typedef std::vector<EntityKey> EntityKeyVector;

struct create_single_edge_impl
{
  typedef void result_type;

  create_single_edge_impl(   size_t & count_edges
                    , std::vector<stk::mesh::EntityId>               & available_ids
                    , impl::edge_map_type        & edge_map
                    , BulkData             & mesh
                    , Entity               & element
                    , unsigned               edge_ordinal
                    , unsigned               num_edge_nodes
                    , const Entity*          edge_nodes
                    , Part * part_to_insert_new_edges
                  )
    : m_count_edges(count_edges)
     ,m_available_ids(available_ids)
    , m_edge_map(edge_map)
    , m_mesh(mesh)
    , m_element(element)
    , m_edge_ordinal(edge_ordinal)
    , m_num_edge_nodes(num_edge_nodes)
    , m_edge_nodes(edge_nodes)
    , m_part_to_insert_new_edges(part_to_insert_new_edges)
  {}

  template <typename Topology>
//  enable if?  check intel
  void operator()(Topology t)
  {
    Topology     elem_topo;
    auto edge_topo = Topology::edge_topology(m_edge_ordinal);

    BulkData & mesh = m_mesh;
    PartVector add_parts;

    add_parts.push_back( & mesh.mesh_meta_data().get_topology_root_part( edge_topo ));
    if (m_part_to_insert_new_edges)
      add_parts.push_back(m_part_to_insert_new_edges);

    Entity elem_nodes[(Topology::num_nodes > 0) ? Topology::num_nodes : 1];
    EntityVector edge_nodes(m_edge_nodes, m_edge_nodes+m_num_edge_nodes);
    OrdinalVector scratch1, scratch2, scratch3;

    Entity ielem = m_element;

    Entity const *nodes = m_mesh.begin_nodes(ielem);
    const int num_nodes = Topology::num_nodes;
    for (int n=0; n<num_nodes; ++n) {
      elem_nodes[n] = nodes[n];
    }

    std::vector<bool> edge_exist(Topology::num_edges,false);
    const int num_edges = m_mesh.num_edges(ielem);
    Entity const *edge_entity = m_mesh.begin_edges(ielem);
    ConnectivityOrdinal const *edge_ords = m_mesh.begin_edge_ordinals(ielem);
    for (int i=0 ; i < num_edges ; ++i)
    {
      if (mesh.is_valid(edge_entity[i]))
      {
        edge_exist[edge_ords[i]] = true;
      }
    }

    if (edge_exist[m_edge_ordinal]) {
      return;
    }

    //sort edge nodes into lexicographical smallest permutation
    if (EntityLess(mesh)(edge_nodes[1], edge_nodes[0])) {
      std::swap(edge_nodes[0], edge_nodes[1]);
    }

    typename impl::edge_map_type::iterator iedge = m_edge_map.find(edge_nodes);

    Entity side;
    Permutation perm = stk::mesh::Permutation::INVALID_PERMUTATION;
    if (iedge == m_edge_map.end()) {
      ThrowRequireMsg(m_count_edges < m_available_ids.size(), "Error: edge generation exhausted available identifier list. Report to sierra-help");
      EntityId edge_id = m_available_ids[m_count_edges];
      m_count_edges++;

      if(mesh.mesh_meta_data().spatial_dimension() == 2)
          side = mesh.declare_solo_side(edge_id, add_parts);
      else
          side = mesh.declare_edge(edge_id, add_parts);

      m_edge_map[edge_nodes] = side;
      for (unsigned n=0u; n<m_num_edge_nodes; ++n)
      {
          Entity node = edge_nodes[n];
          mesh.declare_relation(side,node,n, perm, scratch1, scratch2, scratch3);
      }
    }
    else {
      side = iedge->second;
    }
    perm = mesh.find_permutation(elem_topo, elem_nodes, edge_topo, edge_nodes.data(), m_edge_ordinal);
    ThrowRequireMsg(perm != INVALID_PERMUTATION, "CreateEdges:  could not find valid permutation to connect face to element");
    mesh.declare_relation(ielem, side, m_edge_ordinal, perm, scratch1, scratch2, scratch3);
  }

  //members
  size_t                                          & m_count_edges;
  std::vector<stk::mesh::EntityId>                & m_available_ids;
  impl::edge_map_type         & m_edge_map;
  BulkData              & m_mesh;
  Entity                  m_element;
  unsigned                m_edge_ordinal;
  unsigned                m_num_edge_nodes;
  const Entity*           m_edge_nodes;
  Part                  * m_part_to_insert_new_edges;
};

struct create_edge_impl
{
  typedef void result_type;

  create_edge_impl(   size_t & count_edges
                    , std::vector<stk::mesh::EntityId>               & available_ids
                    , impl::edge_map_type        & edge_map
                    , Bucket               & bucket
                    , Part * part_to_insert_new_edges
                  )
    : m_count_edges(count_edges)
    , m_available_ids(available_ids)
    , m_edge_map(edge_map)
    , m_bucket(bucket)
    , m_part_to_insert_new_edges(part_to_insert_new_edges)
  {}

  template <typename Topology>
  void operator()(Topology t)
  {
    stk::topology elem_topo = m_bucket.topology();

    if (elem_topo.edge_topology(0) == stk::topology::INVALID_TOPOLOGY) {
        return;  // No edges defined for this topology
    }

    BulkData & mesh = m_bucket.mesh();
    PartVector add_parts;
    EntityVector edge_nodes;

    Entity elem_nodes[(Topology::num_nodes > 0) ? Topology::num_nodes : 1];
    OrdinalVector scratch1, scratch2, scratch3;

    for (size_t ielem=0, eelem=m_bucket.size(); ielem<eelem; ++ielem) {
      {
        Entity const *nodes = m_bucket.begin_nodes(ielem);
        const int num_nodes = Topology::num_nodes;
        for (int n=0; n<num_nodes; ++n) {
          elem_nodes[n] = nodes[n];
        }
      }

      std::vector<bool> edge_exist(Topology::num_edges,false);
      const int num_edges = m_bucket.num_edges(ielem);
      Entity const *edge_entity = m_bucket.begin_edges(ielem);
      ConnectivityOrdinal const *edge_ords = m_bucket.begin_edge_ordinals(ielem);
      for (int i=0 ; i < num_edges ; ++i)
      {
        if (mesh.is_valid(edge_entity[i]))
        {
          edge_exist[edge_ords[i]] = true;
        }
      }

      for (unsigned e=0; e != Topology::num_edges; ++e) {

        if (edge_exist[e]) continue;

        auto edge_topo = Topology::edge_topology(e);
        const int num_edge_nodes = edge_topo.num_nodes();
        edge_nodes.resize(num_edge_nodes);
        Topology::edge_nodes(elem_nodes, e, edge_nodes.data());

        add_parts.clear();
        add_parts.push_back( & mesh.mesh_meta_data().get_topology_root_part( edge_topo ));
        if (m_part_to_insert_new_edges)
          add_parts.push_back(m_part_to_insert_new_edges);

        //sort side nodes into lexicographical smallest permutation
        if (EntityLess(mesh)(edge_nodes[1], edge_nodes[0])) {
          std::swap(edge_nodes[0], edge_nodes[1]);
        }

        typename impl::edge_map_type::iterator iedge = m_edge_map.find(edge_nodes);

        Entity side;
        Permutation perm = stk::mesh::Permutation::INVALID_PERMUTATION;
        if (iedge == m_edge_map.end()) {
          ThrowRequireMsg(m_count_edges < m_available_ids.size(), "Error: edge generation exhausted available identifier list. Report to sierra-help");
          EntityId edge_id = m_available_ids[m_count_edges];
          m_count_edges++;

          if(mesh.mesh_meta_data().spatial_dimension() == 2)
              side = mesh.declare_solo_side(edge_id, add_parts);
          else
              side = mesh.declare_edge(edge_id, add_parts);

          m_edge_map[edge_nodes] = side;
          for (int n=0; n<num_edge_nodes; ++n) {
            Entity node = edge_nodes[n];
            mesh.declare_relation(side,node,n, perm, scratch1, scratch2, scratch3);
          }
        }
        else {
          side = iedge->second;
        }
        perm = mesh.find_permutation(elem_topo, elem_nodes, edge_topo, edge_nodes.data(), e);
        ThrowRequireMsg(perm != INVALID_PERMUTATION, "CreateEdges:  could not find valid permutation to connect face to element");
        mesh.declare_relation(m_bucket[ielem], side, e, perm, scratch1, scratch2, scratch3);
      }
    }
  }

  //members
  size_t                                          & m_count_edges;
  std::vector<stk::mesh::EntityId>                & m_available_ids;
  impl::edge_map_type         & m_edge_map;
  Bucket                & m_bucket;
  Part                  * m_part_to_insert_new_edges;
};

struct connect_face_impl
{
  typedef void result_type;

  connect_face_impl(  impl::edge_map_type & edge_map
                    , Bucket        & bucket
                  )
    : m_edge_map(edge_map)
    , m_bucket(bucket)
  {}

  template <typename Topology>
  void operator()(Topology t)
  {
    stk::topology face_topo = m_bucket.topology();

    BulkData & mesh = m_bucket.mesh();
    EntityVector edge_nodes;

    Entity face_nodes[(Topology::num_nodes > 0) ? Topology::num_nodes : 1];
    OrdinalVector scratch1, scratch2, scratch3;

    for (size_t iface=0, eface=m_bucket.size(); iface<eface; ++iface) {
      {
        Entity const *nodes = m_bucket.begin_nodes(iface);
        const int num_nodes = Topology::num_nodes;
        for (int n=0; n<num_nodes; ++n) {
          face_nodes[n] = nodes[n];
        }
      }

      std::vector<bool> edge_exist(Topology::num_edges,false);
      const int num_edges = m_bucket.num_edges(iface);
      Entity const * edge_entity = m_bucket.begin_edges(iface);
      ConnectivityOrdinal const *edge_ords = m_bucket.begin_edge_ordinals(iface);
      for (int i=0 ; i < num_edges ; ++i)
      {
        if (mesh.is_valid(edge_entity[i]))
        {
          edge_exist[edge_ords[i]] = true;
        }
      }

      for (unsigned e=0; e != Topology::num_edges; ++e) {

        if (edge_exist[e]) continue;

        auto edge_topo = Topology::edge_topology(e);
        edge_nodes.resize(edge_topo.num_nodes());
        Topology::edge_nodes(face_nodes, e, edge_nodes.data());

        //sort edge nodes into lexicographical smallest permutation
        if (EntityLess(mesh)(edge_nodes[1], edge_nodes[0])) {
          std::swap(edge_nodes[0], edge_nodes[1]);
        }

        //the edge should already exist
        typename impl::edge_map_type::iterator iedge = m_edge_map.find(edge_nodes);

        //if this fails, we don't have the correct edge made yet
        //which is fine
        if (iedge != m_edge_map.end()) {
          Entity edge = iedge->second;
          Permutation perm = mesh.find_permutation(face_topo, face_nodes, edge_topo, edge_nodes.data(), e);
          ThrowRequireMsg(perm != INVALID_PERMUTATION, "CreateEdges:  could not find valid permutation to connect face to element");
          mesh.declare_relation(m_bucket[iface], edge, e, perm, scratch1, scratch2, scratch3);
        }
      }
    }
  }

  //members
  impl::edge_map_type         & m_edge_map;
  Bucket                & m_bucket;
};

} //empty namespace

namespace impl {

  void connect_faces_to_edges(BulkData & mesh,
                              const Selector & element_selector,
                              impl::edge_map_type edge_map) {
      // connect existing faces to edges
      if (mesh.mesh_meta_data().spatial_dimension() == 3u) {

          BucketVector const& face_buckets = mesh.get_buckets(stk::topology::FACE_RANK, element_selector & (mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part()));

          for (size_t i=0, e=face_buckets.size(); i<e; ++i) {
              Bucket &b = *face_buckets[i];

              connect_face_impl functor(edge_map, b);
              stk::topology::apply_host_functor< connect_face_impl > apply(functor);
              apply( b.topology() );
          }
      }
  }
} //namespace impl

void create_edges( BulkData & mesh )
{
  create_edges(mesh, mesh.mesh_meta_data().universal_part(), nullptr );
}

void create_edges( BulkData & mesh, const Selector & element_selector, Part * part_to_insert_new_edges )
{
  std::vector<stk::mesh::EntityId> ids_requested;

  std::vector<size_t> localEntityCounts;
  stk::mesh::count_entities(element_selector, mesh, localEntityCounts);
  unsigned guessMultiplier = 12;
  unsigned numRequested = localEntityCounts[stk::topology::ELEMENT_RANK] * guessMultiplier;
  mesh.generate_new_ids(stk::topology::EDGE_RANK, numRequested, ids_requested);
  size_t count_edges = 0;

  bool i_started = mesh.modification_begin();

  {
    {
      impl::edge_map_type        edge_map;
      //populate the edge_map with existing edges
      {
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

      // create edges and connect them to elements
      {
        BucketVector const& element_buckets = mesh.get_buckets(stk::topology::ELEMENT_RANK, element_selector & mesh.mesh_meta_data().locally_owned_part());

        //create the edges for the elements in each bucket
        for (size_t i=0, e=element_buckets.size(); i<e; ++i) {
          Bucket &b = *element_buckets[i];

          create_edge_impl functor( count_edges, ids_requested, edge_map, b, part_to_insert_new_edges);
          stk::topology::apply_host_functor< create_edge_impl > apply(functor);
          apply( b.topology() );
        }
      }

      // create any edges that are shared between locally-owned elements and elements that are
      // in element_selector but are not locally-owned
      {
        BucketVector const& element_buckets = mesh.get_buckets(stk::topology::ELEMENT_RANK, element_selector & !mesh.mesh_meta_data().locally_owned_part());

        std::vector<Entity> elements;
        for(size_t i=0, e=element_buckets.size(); i<e; ++i) {
          Bucket& b = *element_buckets[i];
          stk::topology elemTopology = b.topology();
          ThrowRequireMsg(elemTopology != stk::topology::INVALID_TOPOLOGY, "create_edges ERROR, element bucket with invalid topology.");
          const unsigned numEdgesPerElem = elemTopology.num_edges();
          for(size_t j=0, jend=b.size(); j<jend; ++j) {
            const Entity* elemNodes = b.begin_nodes(j);
            Entity edgeNodes[3];
            Entity localElemEdgeNodes[3];
            for(unsigned edgeIndex=0; edgeIndex<numEdgesPerElem; ++edgeIndex) {
              const unsigned numNodesPerEdge = elemTopology.edge_topology(edgeIndex).num_nodes();
              elemTopology.edge_nodes(elemNodes, edgeIndex, edgeNodes);
              stk::mesh::impl::find_locally_owned_elements_these_nodes_have_in_common(mesh, numNodesPerEdge, edgeNodes, elements);
              if (!elements.empty()) {
                for(size_t el=0; el<elements.size(); ++el) {
                  Entity localElem = elements[el];
                  unsigned localElemEdgeOrdinal = 10000;
                  stk::mesh::impl::find_element_edge_ordinal_and_equivalent_nodes(mesh, localElem, numNodesPerEdge, edgeNodes, localElemEdgeOrdinal, localElemEdgeNodes);
                  create_single_edge_impl functor( count_edges, ids_requested, edge_map, mesh, localElem, localElemEdgeOrdinal, numNodesPerEdge, localElemEdgeNodes, part_to_insert_new_edges);
                  stk::topology::apply_host_functor< create_single_edge_impl > apply(functor);
                  apply( mesh.bucket(localElem).topology() );
                }
                elements.clear();
              }
            }
          }
        }
      }

      impl::connect_faces_to_edges(mesh, mesh.mesh_meta_data().universal_part(), edge_map);
    }
  }

  if(i_started)
  {
    bool oldOption = mesh.use_entity_ids_for_resolving_sharing();
    mesh.set_use_entity_ids_for_resolving_sharing(false);
    std::vector<EntityRank> entity_rank_vector = {stk::topology::EDGE_RANK};
    mesh.modification_end_for_entity_creation( entity_rank_vector );
    mesh.set_use_entity_ids_for_resolving_sharing(oldOption);
  }
}

}
}
