/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

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

struct shared_edge_type
{
  stk::topology::topology_t topology;
  EntityKey                 nodes[2];
  EntityKey                 local_key;
  EntityKey                 global_key;

  friend inline bool operator < (shared_edge_type const& l, shared_edge_type const& r)
  {
    if (l.topology < r.topology)   return true;
    if (l.topology > r.topology)   return false;
    if ( l.nodes[0] < r.nodes[0] ) return true;
    if ( l.nodes[0] > r.nodes[0] ) return false;

    return l.nodes[1] < r.nodes[1];
  }

  friend inline bool operator == (shared_edge_type const& l, shared_edge_type const& r)
  {
    return    (l.topology == r.topology)
           && (l.nodes[0] == r.nodes[0])
           && (l.nodes[1] == r.nodes[1]);
  }
};

typedef boost::unordered_map<EntityVector,Entity> edge_map_type;
typedef std::vector<EntityKey> EntityKeyVector;
typedef std::vector< shared_edge_type > shared_edge_map_type;

struct create_edge_impl
{
  typedef void result_type;

  create_edge_impl(   size_t               & next_edge
                    , edge_map_type        & edge_map
                    , shared_edge_map_type & shared_edge_map
                    , Bucket               & bucket
                  )
    : m_next_edge(next_edge)
    , m_edge_map(edge_map)
    , m_shared_edge_map(shared_edge_map)
    , m_bucket(bucket)
  {}

  template <typename Topology>
  typename boost::enable_if_c< (Topology::num_edges > 0u), void>::type
  operator()(Topology t)
  {
    typedef topology::topology_type< Topology::edge_topology> EdgeTopology;

    BulkData & mesh = m_bucket.mesh();
    PartVector add_parts;

    add_parts.push_back( & mesh.mesh_meta_data().get_cell_topology_root_part( get_cell_topology( EdgeTopology::value )));

    boost::array<Entity,Topology::num_nodes> elem_nodes;
    EntityVector edge_nodes(EdgeTopology::num_nodes);
    EntityKeyVector edge_node_keys(EdgeTopology::num_nodes);

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
        if (iedge == m_edge_map.end()) {
          EntityId edge_id = m_next_edge++;
          edge = mesh.declare_entity( stk::topology::EDGE_RANK, edge_id, add_parts);
          m_edge_map[edge_nodes] = edge;
          bool shared_edge = true;
          const int num_edge_nodes = EdgeTopology::num_nodes;
          for (int n=0; n<num_edge_nodes; ++n) {
            Entity node = edge_nodes[n];
            mesh.declare_relation(edge,node,n);
            shared_edge = shared_edge && mesh.bucket(node).shared();
          }
          if (shared_edge) {
            shared_edge_type sedge;
            sedge.topology = EdgeTopology::value;
            for (int n=0; n<num_edge_nodes; ++n) {
              sedge.nodes[n] = mesh.entity_key(edge_nodes[n]);
            }
            const EntityKey &edge_key = mesh.entity_key(edge);
            sedge.local_key   = edge_key;
            sedge.global_key  = edge_key;
            m_shared_edge_map.push_back( sedge );
          }
        }
        else {
          edge = iedge->second;
        }
        mesh.declare_relation(m_bucket[ielem], edge, e);
      }
    }
  }

  template <typename Topology>
  typename boost::enable_if_c< (Topology::num_edges == 0u), void>::type
  operator()(Topology t)
  {}


  //members
  size_t                & m_next_edge;
  edge_map_type         & m_edge_map;
  shared_edge_map_type  & m_shared_edge_map;
  Bucket                & m_bucket;
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
        mesh.declare_relation(m_bucket[iface], edge, e);
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

void update_shared_edges_global_ids( BulkData & mesh, shared_edge_map_type & shared_edge_map)
{
  //sort the edges for ease of lookup
  std::sort(shared_edge_map.begin(), shared_edge_map.end());

  std::vector< shared_edge_map_type >  shared_edges( mesh.parallel_size());

  for (shared_edge_map_type::const_iterator itr = shared_edge_map.begin(),
                                            end = shared_edge_map.end();
                                            itr != end;
                                            ++itr )
  {
    //find process that share this edge;

    PairIterEntityComm left_shared = mesh.entity_comm_sharing( itr->nodes[0] );
    PairIterEntityComm right_shared = mesh.entity_comm_sharing( itr->nodes[1] );

    EntityCommInfoVector shared_processes;

    std::set_intersection( left_shared.first, left_shared.second,
                           right_shared.first, right_shared.second,
                           std::back_inserter(shared_processes)
                         );

    for (EntityCommInfoVector::const_iterator comm_itr = shared_processes.begin(),
                                        comm_end = shared_processes.end();
                                        comm_itr != comm_end;
                                        ++comm_itr )
    {
      if (comm_itr->proc > mesh.parallel_rank())
        shared_edges[comm_itr->proc].push_back(*itr);
    }
  }

  CommAll comm( mesh.parallel() );

  //pack send buffers
  for (int allocation_pass=0; allocation_pass<2; ++allocation_pass) {
    if (allocation_pass==1) {
     comm.allocate_buffers( mesh.parallel_size() /4, 0);
    }

    for (int proc=mesh.parallel_rank()+1, parallel_size = mesh.parallel_size(); proc<parallel_size; ++proc) {
      for (size_t e=0, num_shared = shared_edges[proc].size(); e < num_shared; ++e) {
        shared_edge_type const & sedge = shared_edges[proc][e];
        comm.send_buffer(proc).pack<stk::topology::topology_t>(sedge.topology)
        .pack<EntityKey>(sedge.nodes[0])
        .pack<EntityKey>(sedge.nodes[1])
        .pack<EntityKey>(sedge.local_key);
      }
    }
  }

  comm.communicate();

  for ( unsigned ip = mesh.parallel_rank() ; ip > 0 ; ) {
    --ip;
    CommBuffer & buf = comm.recv_buffer( ip );
    while ( buf.remaining() ) {
      shared_edge_type sedge;

      buf.unpack<stk::topology::topology_t>(sedge.topology)
         .unpack<EntityKey>(sedge.nodes[0])
         .unpack<EntityKey>(sedge.nodes[1])
         .unpack<EntityKey>(sedge.global_key);

      shared_edge_map_type::iterator shared_itr = std::lower_bound(shared_edge_map.begin(), shared_edge_map.end(), sedge);

      //update the global global_key
      if (shared_itr != shared_edge_map.end() && *shared_itr == sedge) {
        shared_itr->global_key = sedge.global_key;
      }
    }
  }

  //update the entity keys
  for (size_t i=0, e=shared_edge_map.size(); i<e; ++i) {
    if (shared_edge_map[i].global_key != shared_edge_map[i].local_key) {
      mesh.internal_change_entity_key(shared_edge_map[i].local_key, shared_edge_map[i].global_key, mesh.get_entity(shared_edge_map[i].local_key));
    }
  }
}


void create_edges( BulkData & mesh )
{
  create_edges(mesh, mesh.mesh_meta_data().universal_part());
}

void create_edges( BulkData & mesh, const Selector & element_selector )
{

  //  static size_t next_edge = static_cast<size_t>(mesh.parallel_rank()+1) << 32;
  // NOTE: This is a workaround to eliminate some bad behavior with the equation above when
  //       the #proc is a power of two.  The 256 below is the bin size of the Distributed Index.
  static size_t next_edge = (static_cast<size_t>(mesh.parallel_rank()+1) << 32) + 256 * mesh.parallel_rank();

  mesh.modification_begin();

  {
    shared_edge_map_type shared_edge_map;

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

          create_edge_impl functor( next_edge, edge_map, shared_edge_map, b);
          stk::topology::apply_functor< create_edge_impl > apply(functor);
          apply( b.topology() );
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

    //update global ids of shared edges
    update_shared_edges_global_ids( mesh, shared_edge_map );
  }

  mesh.modification_end( BulkData::MOD_END_COMPRESS_AND_SORT );
}

}
}
