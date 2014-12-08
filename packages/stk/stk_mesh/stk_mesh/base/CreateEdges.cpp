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

#include "stk_mesh/base/CreateEdges.hpp"

#include <stddef.h>                     // for size_t
#include <algorithm>                    // for swap, lower_bound, max, etc
#include <boost/array.hpp>              // for array
#include <functional>                   // for equal_to
#include <iterator>                     // for back_insert_iterator, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, EntityLess, etc
#include <stk_mesh/base/Entity.hpp>     // for Entity, hash_value
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, get_cell_topology
#include <stk_mesh/base/Selector.hpp>   // for operator&, Selector, etc
#include <stk_mesh/base/Types.hpp>      // for EntityVector, etc
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_util/parallel/ParallelComm.hpp>  // for CommBuffer, CommAll
#include <vector>                       // for vector, etc
#include "boost/functional/hash/extensions.hpp"  // for hash
#include "boost/tuple/detail/tuple_basic.hpp"  // for get
#include "boost/unordered/detail/buckets.hpp"  // for iterator, etc
#include "boost/unordered/unordered_map.hpp"
#include "boost/utility/enable_if.hpp"  // for enable_if_c
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_topology/apply_functor.tcc"  // for topology::apply_functor
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_topology/topology.tcc"    // for topology::num_nodes
#include "stk_topology/topology_type.tcc"  // for topology::topology_type
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowAssert
#include "stk_util/util/NamedPair.hpp"  // for EntityCommInfo::operator=, etc

namespace stk {
namespace mesh {

namespace {

typedef boost::unordered_map<EntityVector,Entity> edge_map_type;
typedef std::vector<EntityKey> EntityKeyVector;

struct create_single_edge_impl
{
  typedef void result_type;

  create_single_edge_impl(   size_t & count_edges
                    , std::vector<stk::mesh::EntityId>               & available_ids
                    , edge_map_type        & edge_map
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
  typename boost::enable_if_c< (Topology::num_edges > 0u), void>::type
  operator()(Topology t)
  {
    typedef topology::topology_type< Topology::edge_topology> EdgeTopology;

    BulkData & mesh = m_mesh;
    PartVector add_parts;

    add_parts.push_back( & mesh.mesh_meta_data().get_cell_topology_root_part( get_cell_topology( EdgeTopology::value )));
    if (m_part_to_insert_new_edges)
      add_parts.push_back(m_part_to_insert_new_edges);

    boost::array<Entity,Topology::num_nodes> elem_nodes;
    EntityVector edge_nodes(m_edge_nodes, m_edge_nodes+m_num_edge_nodes);
    OrdinalVector ordinal_scratch;
    ordinal_scratch.reserve(64);
    PartVector part_scratch;
    part_scratch.reserve(64);

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

    typename edge_map_type::iterator iedge = m_edge_map.find(edge_nodes);

    Entity edge;
    Permutation perm = static_cast<Permutation>(0);
    if (iedge == m_edge_map.end()) {
      ThrowRequireMsg(m_count_edges < m_available_ids.size(), "Error: edge generation exhausted available identifier list. Report to sierra-help");
      EntityId edge_id = m_available_ids[m_count_edges];
      m_count_edges++;

      edge = mesh.declare_entity( stk::topology::EDGE_RANK, edge_id, add_parts);
      m_edge_map[edge_nodes] = edge;
      const int num_edge_nodes = EdgeTopology::num_nodes;
      for (int n=0; n<num_edge_nodes; ++n) {
        Entity node = edge_nodes[n];
        mesh.declare_relation(edge,node,n, perm, ordinal_scratch, part_scratch);
      }
    }
    else {
      edge = iedge->second;
    }
    mesh.declare_relation(ielem, edge, m_edge_ordinal, perm, ordinal_scratch, part_scratch);
  }

  template <typename Topology>
  typename boost::enable_if_c< (Topology::num_edges == 0u), void>::type
  operator()(Topology t)
  {}


  //members
  size_t                                          & m_count_edges;
  std::vector<stk::mesh::EntityId>                & m_available_ids;
  edge_map_type         & m_edge_map;
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
                    , edge_map_type        & edge_map
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
  typename boost::enable_if_c< (Topology::num_edges > 0u), void>::type
  operator()(Topology t)
  {
    typedef topology::topology_type< Topology::edge_topology> EdgeTopology;

    BulkData & mesh = m_bucket.mesh();
    PartVector add_parts;

    add_parts.push_back( & mesh.mesh_meta_data().get_cell_topology_root_part( get_cell_topology( EdgeTopology::value )));
    if (m_part_to_insert_new_edges)
      add_parts.push_back(m_part_to_insert_new_edges);

    boost::array<Entity,Topology::num_nodes> elem_nodes;
    EntityVector edge_nodes(EdgeTopology::num_nodes);
    EntityKeyVector edge_node_keys(EdgeTopology::num_nodes);
    OrdinalVector ordinal_scratch;
    ordinal_scratch.reserve(64);
    PartVector part_scratch;
    part_scratch.reserve(64);

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

      for (unsigned e=0; e < Topology::num_edges; ++e) {

        if (edge_exist[e]) continue;

        Topology::edge_nodes(elem_nodes, e, edge_nodes.begin());

        //sort edge nodes into lexicographical smallest permutation
        if (EntityLess(mesh)(edge_nodes[1], edge_nodes[0])) {
          std::swap(edge_nodes[0], edge_nodes[1]);
        }

        typename edge_map_type::iterator iedge = m_edge_map.find(edge_nodes);

        Entity edge;
        Permutation perm = static_cast<Permutation>(0);
        if (iedge == m_edge_map.end()) {
          ThrowRequireMsg(m_count_edges < m_available_ids.size(), "Error: edge generation exhausted available identifier list. Report to sierra-help");
          EntityId edge_id = m_available_ids[m_count_edges];
          m_count_edges++;

          edge = mesh.declare_entity( stk::topology::EDGE_RANK, edge_id, add_parts);
          m_edge_map[edge_nodes] = edge;
          const int num_edge_nodes = EdgeTopology::num_nodes;
          for (int n=0; n<num_edge_nodes; ++n) {
            Entity node = edge_nodes[n];
            mesh.declare_relation(edge,node,n, perm, ordinal_scratch, part_scratch);
          }
        }
        else {
          edge = iedge->second;
        }
        mesh.declare_relation(m_bucket[ielem], edge, e, perm, ordinal_scratch, part_scratch);
      }
    }
  }

  template <typename Topology>
  typename boost::enable_if_c< (Topology::num_edges == 0u), void>::type
  operator()(Topology t)
  {}


  //members
  size_t                                          & m_count_edges;
  std::vector<stk::mesh::EntityId>                & m_available_ids;
  edge_map_type         & m_edge_map;
  Bucket                & m_bucket;
  Part                  * m_part_to_insert_new_edges;
};

struct connect_face_impl
{
  typedef void result_type;

  connect_face_impl(  edge_map_type & edge_map
                    , Bucket        & bucket
                  )
    : m_edge_map(edge_map)
    , m_bucket(bucket)
  {}

  template <typename Topology>
  typename boost::enable_if_c< (Topology::num_edges > 0u), void>::type
  operator()(Topology t)
  {
    typedef topology::topology_type< Topology::edge_topology> EdgeTopology;

    BulkData & mesh = m_bucket.mesh();

    boost::array<Entity,Topology::num_nodes> face_nodes;
    EntityVector edge_nodes(EdgeTopology::num_nodes);
    OrdinalVector ordinal_scratch;
    ordinal_scratch.reserve(64);
    PartVector part_scratch;
    part_scratch.reserve(64);

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

      for (unsigned e=0; e < Topology::num_edges; ++e) {

        if (edge_exist[e]) continue;

        Topology::edge_nodes(face_nodes, e, edge_nodes.begin());

        //sort edge nodes into lexicographical smallest permutation
        if (EntityLess(mesh)(edge_nodes[1], edge_nodes[0])) {
          std::swap(edge_nodes[0], edge_nodes[1]);
        }

        //the edge should already exist
        typename edge_map_type::iterator iedge = m_edge_map.find(edge_nodes);

        ThrowAssert(iedge != m_edge_map.end());

        Entity edge = iedge->second;
        Permutation perm = static_cast<Permutation>(0);
        mesh.declare_relation(m_bucket[iface], edge, e, perm, ordinal_scratch, part_scratch);
      }
    }
  }

  template <typename Topology>
  typename boost::enable_if_c< (Topology::num_edges == 0u), void>::type
  operator()(Topology t)
  {}

  //members
  edge_map_type         & m_edge_map;
  Bucket                & m_bucket;
};

} //namespace

void create_edges( BulkData & mesh )
{
  create_edges(mesh, mesh.mesh_meta_data().universal_part(), 0 );
}

void create_edges( BulkData & mesh, const Selector & element_selector, Part * part_to_insert_new_edges )
{

  //  static size_t next_edge = static_cast<size_t>(mesh.parallel_rank()+1) << 32;
  // NOTE: This is a workaround to eliminate some bad behavior with the equation above when
  //       the #proc is a power of two.  The 256 below is the bin size of the Distributed Index.

  std::vector<stk::mesh::EntityId> ids_requested;

  std::vector<unsigned> localEntityCounts;
  stk::mesh::count_entities(element_selector, mesh, localEntityCounts);
  unsigned guessMultiplier = 12;
  unsigned numRequested = localEntityCounts[stk::topology::ELEMENT_RANK] * guessMultiplier;
  mesh.generate_new_ids(stk::topology::EDGE_RANK, numRequested, ids_requested);
  size_t count_edges = 0;

  bool i_started = mesh.modification_begin();

  {
    {
      edge_map_type        edge_map;
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
          stk::topology::apply_functor< create_edge_impl > apply(functor);
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
          const unsigned numNodesPerEdge = elemTopology.sub_topology(stk::topology::EDGE_RANK).num_nodes();
          for(size_t j=0, jend=b.size(); j<jend; ++j) {
            const Entity* elemNodes = b.begin_nodes(j);
            Entity edgeNodes[3];
            Entity localElemEdgeNodes[3];
            for(unsigned edgeIndex=0; edgeIndex<numEdgesPerElem; ++edgeIndex) {
              elemTopology.edge_nodes(elemNodes, edgeIndex, edgeNodes);
              stk::mesh::impl::find_locally_owned_elements_these_nodes_have_in_common(mesh, numNodesPerEdge, edgeNodes, elements);
              if (!elements.empty()) {
                for(size_t el=0; el<elements.size(); ++el) {
                  Entity localElem = elements[el];
                  unsigned localElemEdgeOrdinal = 10000;
                  stk::mesh::impl::find_element_edge_ordinal_and_equivalent_nodes(mesh, localElem, numNodesPerEdge, edgeNodes, localElemEdgeOrdinal, localElemEdgeNodes);
                  create_single_edge_impl functor( count_edges, ids_requested, edge_map, mesh, localElem, localElemEdgeOrdinal, numNodesPerEdge, localElemEdgeNodes, part_to_insert_new_edges);
                  stk::topology::apply_functor< create_single_edge_impl > apply(functor);
                  apply( b.topology() );
                }
                elements.clear();
              }
            }
          }
        }
      }

      // connect existing faces to edges
      if (mesh.mesh_meta_data().spatial_dimension() == 3u) {

        BucketVector const& face_buckets = mesh.get_buckets(stk::topology::FACE_RANK, element_selector & (mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part()));

        //create the edges for the faces in each bucket
        for (size_t i=0, e=face_buckets.size(); i<e; ++i) {
          Bucket &b = *face_buckets[i];

          connect_face_impl functor(edge_map, b);
          stk::topology::apply_functor< connect_face_impl > apply(functor);
          apply( b.topology() );
        }
      }
    }
  }

  if(i_started)
  {
    bool oldOption = mesh.use_entity_ids_for_resolving_sharing();
    mesh.set_use_entity_ids_for_resolving_sharing(false);
    mesh.modification_end_for_entity_creation( stk::topology::EDGE_RANK, BulkData::MOD_END_COMPRESS_AND_SORT );
    mesh.set_use_entity_ids_for_resolving_sharing(oldOption);
  }
}

}
}
