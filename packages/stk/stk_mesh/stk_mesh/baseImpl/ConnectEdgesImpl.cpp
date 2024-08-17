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

#include <stk_mesh/baseImpl/ConnectEdgesImpl.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FindPermutation.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include "stk_topology/topology.hpp"    // for topology, etc

namespace stk {
namespace mesh {

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
          Entity const* original_edge_nodes = mesh.begin_nodes(edge);
          Permutation perm = stk::mesh::find_permutation(mesh, face_topo, face_nodes, edge_topo, original_edge_nodes, e);
          STK_ThrowRequireMsg(perm != INVALID_PERMUTATION, "CreateEdges:  could not find valid permutation to connect face to edge");
          mesh.declare_relation(m_bucket[iface], edge, e, perm, scratch1, scratch2, scratch3);
        }
      }
    }
  }

  //members
  impl::edge_map_type         & m_edge_map;
  Bucket                & m_bucket;
};

struct connect_face_entity_impl
{
  typedef void result_type;

  connect_face_entity_impl(impl::edge_map_type & edge_map,
                           BulkData& bulk,
                           stk::mesh::Entity face)
    : m_edge_map(edge_map)
    , m_bulk(bulk)
    , m_face(face)
  {}

  template<typename Topology>
  void operator()(Topology t)
  {
    stk::topology face_topo = m_bulk.bucket(m_face).topology();

    EntityVector edge_nodes;

    Entity face_nodes[(Topology::num_nodes > 0) ? Topology::num_nodes : 1];
    OrdinalVector scratch1, scratch2, scratch3;

    Entity const *nodes = m_bulk.begin_nodes(m_face);
    const int num_nodes = Topology::num_nodes; 
    for (int n=0; n<num_nodes; ++n) {
      face_nodes[n] = nodes[n];
    }

    std::vector<bool> edge_exist(Topology::num_edges,false);
    const int num_edges = m_bulk.num_edges(m_face);
    Entity const * edge_entity = m_bulk.begin_edges(m_face);
    ConnectivityOrdinal const *edge_ords = m_bulk.begin_edge_ordinals(m_face);
    for (int i=0 ; i < num_edges ; ++i)
    {
      if (m_bulk.is_valid(edge_entity[i]))
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
      if (EntityLess(m_bulk)(edge_nodes[1], edge_nodes[0])) {
        std::swap(edge_nodes[0], edge_nodes[1]);
      }

      //the edge should already exist
      typename impl::edge_map_type::iterator iedge = m_edge_map.find(edge_nodes);

      //if this fails, we don't have the correct edge made yet
      //which is fine
      if (iedge != m_edge_map.end()) {
        Entity edge = iedge->second;
        Entity const* original_edge_nodes = m_bulk.begin_nodes(edge);

        Permutation perm = stk::mesh::find_permutation(m_bulk, face_topo, face_nodes, edge_topo, original_edge_nodes, e);
        STK_ThrowRequireMsg(perm != INVALID_PERMUTATION, "Connect face to edge:  could not find valid permutation to connect face to edge");
        m_bulk.declare_relation(m_face, edge, e, perm, scratch1, scratch2, scratch3);
      }
    }
  }

  //members
  impl::edge_map_type         & m_edge_map;
  BulkData& m_bulk;
  Entity m_face;
};

namespace impl {

  void fill_edge_map_for_face(BulkData& bulk, Entity face, impl::edge_map_type& edge_map) {
    if(bulk.mesh_meta_data().spatial_dimension() != 3) { return; }
    if(bulk.entity_rank(face) != stk::topology::FACE_RANK) { return; }

    edge_map.clear();

    stk::topology sideTopo = bulk.bucket(face).topology();
    unsigned numEdges = sideTopo.num_edges();
    const Entity* sideNodes = bulk.begin_nodes(face);
    EntityVector edgeNodes;
    EntityVector commonEdges;

    for(unsigned i = 0; i < numEdges; i++) {
      stk::topology edgeTopo = sideTopo.edge_topology(i);
      edgeNodes.resize(edgeTopo.num_nodes());
      sideTopo.sub_topology_nodes(sideNodes, stk::topology::EDGE_RANK, i, edgeNodes.data());
      impl::find_entities_these_nodes_have_in_common(bulk, stk::topology::EDGE_RANK, 2, edgeNodes.data(), commonEdges);

      //sort edge nodes into lexicographical smallest permutation
      if (EntityLess(bulk)(edgeNodes[1], edgeNodes[0])) {
        std::swap(edgeNodes[0], edgeNodes[1]);
      }

      for(auto edge : commonEdges) {
        edge_map[edgeNodes] = edge;
      }
    }
  }

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

  void connect_face_to_edges(BulkData & mesh,
                             Entity face,
                             impl::edge_map_type edge_map) {

      if (mesh.mesh_meta_data().spatial_dimension() == 3u) {

        Bucket& b = mesh.bucket(face);

        connect_face_entity_impl functor(edge_map, mesh, face);
        stk::topology::apply_host_functor< connect_face_entity_impl > apply(functor);
        apply( b.topology() );
      }
  }

  void connect_face_to_edges(BulkData& bulk, Entity face)
  {
    impl::edge_map_type edge_map;
    fill_edge_map_for_face(bulk, face, edge_map);
    connect_face_to_edges(bulk, face, edge_map);
  }
} //namespace impl
}
}
