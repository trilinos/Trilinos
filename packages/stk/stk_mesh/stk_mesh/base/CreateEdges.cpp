#include <map>
#include <set>
#include <algorithm>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/CellTopology.hpp>
#include <stk_mesh/base/CreateAdjacentEntities.hpp>
#include <stk_mesh/base/BoundaryAnalysis.hpp>

#include <stk_util/parallel/ParallelComm.hpp>


#include <boost/array.hpp>
#include <boost/unordered_map.hpp>

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

  friend inline bool operator > (shared_edge_type const& l, shared_edge_type const& r)
  {
    return r < l;
  }

  friend inline bool operator == (shared_edge_type const& l, shared_edge_type const& r)
  {
    return    (l.topology == r.topology)
           && (l.nodes[0] == r.nodes[0])
           && (l.nodes[1] == r.nodes[1]);
  }

  friend inline bool operator != (shared_edge_type const& l, shared_edge_type const& r)
  {
    return !(l==r);
  }

  friend inline bool operator <= (shared_edge_type const& l, shared_edge_type const& r)
  {
    return (l<r) || (l==r);
  }

  friend inline bool operator >= (shared_edge_type const& l, shared_edge_type const& r)
  {
    return (l>r) || (l==r);
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
  void operator()(Topology t)
  {
    typedef topology::topology_type< Topology::edge_topology> EdgeTopology;

    BulkData & mesh = m_bucket.mesh();
    PartVector add_parts;

    add_parts.push_back( & mesh.mesh_meta_data().get_cell_topology_root_part( get_cell_topology( EdgeTopology::value )));


    boost::array<Entity,Topology::num_nodes> elem_nodes;
    EntityVector edge_nodes(EdgeTopology::num_nodes);
    EntityKeyVector edge_node_keys(EdgeTopology::num_nodes);

    for (size_t ielem=0, eelem=m_bucket.size(); ielem<eelem; ++ielem) {
      Entity elem = m_bucket[ielem];
      {
        PairIterRelation nodes = elem.node_relations();
        for (int n=0; n<Topology::num_nodes; ++n) {
          elem_nodes[n] = nodes[n].entity();
        }
      }

      PairIterRelation elem_edges = elem.relations(stk::topology::EDGE_RANK);

      for (int e=0; e < Topology::num_edges; ++e) {

        while(!elem_edges.empty() && elem_edges->relation_ordinal() < static_cast<RelationIdentifier>(e) ) {
          ++elem_edges;
        }

        //element already has this edge defined
        if (!elem_edges.empty() && elem_edges->relation_ordinal() == static_cast<RelationIdentifier>(e) ) {
          continue;
        }

        Topology::edge_nodes(elem_nodes, e, edge_nodes.begin());

        //sort edge nodes into lexicographical smallest permutation
        if (edge_nodes[1].key() < edge_nodes[0].key()) {
          std::swap(edge_nodes[0], edge_nodes[1]);
        }

        typename edge_map_type::iterator iedge = m_edge_map.find(edge_nodes);

        Entity edge;
        if (iedge == m_edge_map.end()) {
          EntityId edge_id = m_next_edge++;
          edge = mesh.declare_entity( stk::topology::EDGE_RANK, edge_id, add_parts);
          m_edge_map[edge_nodes] = edge;
          bool shared_edge = true;
          for (int n=0; n<EdgeTopology::num_nodes; ++n) {
            Entity node = edge_nodes[n];
            mesh.declare_relation(edge,node,n);
            shared_edge = shared_edge && node.bucket().shared();
          }
          if (shared_edge) {
            shared_edge_type sedge;
            sedge.topology = EdgeTopology::value;
            for (int n=0; n<EdgeTopology::num_nodes; ++n) {
              sedge.nodes[n] = edge_nodes[n].key();
            }
            sedge.local_key   = edge.key();
            sedge.global_key  = edge.key();
            m_shared_edge_map.push_back( sedge );
          }
        }
        else {
          edge = iedge->second;
        }
        mesh.declare_relation(elem,edge,e);
      }
    }


  }

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
  void operator()(Topology t)
  {
    typedef topology::topology_type< Topology::edge_topology> EdgeTopology;

    BulkData & mesh = m_bucket.mesh();

    boost::array<Entity,Topology::num_nodes> face_nodes;
    EntityVector edge_nodes(EdgeTopology::num_nodes);

    for (size_t iface=0, eface=m_bucket.size(); iface<eface; ++iface) {
      Entity face = m_bucket[iface];
      {
        PairIterRelation nodes = face.node_relations();
        for (int n=0; n<Topology::num_nodes; ++n) {
          face_nodes[n] = nodes[n].entity();
        }
      }

      PairIterRelation face_edges = face.relations(stk::topology::EDGE_RANK);

      for (int e=0; e < Topology::num_edges; ++e) {

        while(!face_edges.empty() && face_edges->relation_ordinal() < static_cast<RelationIdentifier>(e) ) {
          ++face_edges;
        }

        //faceent already has this edge defined
        if (!face_edges.empty() && face_edges->relation_ordinal() == static_cast<RelationIdentifier>(e) ) {
          continue;
        }

        Topology::edge_nodes(face_nodes, e, edge_nodes.begin());

        //sort edge nodes into lexicographical smallest permutation
        if (edge_nodes[1] < edge_nodes[0]) {
          std::swap(edge_nodes[0], edge_nodes[1]);
        }

        //the edge should already exist
        typename edge_map_type::iterator iedge = m_edge_map.find(edge_nodes);

        ThrowAssert(iedge != m_edge_map.end());

        Entity edge = iedge->second;
        mesh.declare_relation(face,edge,e);
      }
    }


  }

  //members
  edge_map_type         & m_edge_map;
  Bucket                & m_bucket;
};

} //namespace

void update_shared_edges_global_ids( BulkData & mesh, shared_edge_map_type & shared_edge_map, impl::EntityRepository & entity_repo)
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
      entity_repo.update_entity_key(shared_edge_map[i].global_key, shared_edge_map[i].local_key);
    }
  }

}



void create_edges( BulkData & mesh, const Selector & element_selector )
{

  static size_t next_edge = static_cast<size_t>(mesh.parallel_rank()+1) << 32;

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

          const int num_nodes = b.topology().num_nodes();
          EntityVector edge_nodes(num_nodes);

          for (size_t j=0, je=b.size(); j<je; ++j) {
            Entity edge = b[j];
            PairIterRelation nodes_rel = edge.node_relations();

            for (int n=0; n<num_nodes; ++n) {
              edge_nodes[n] = nodes_rel[n].entity();
            }

            edge_map[edge_nodes] = edge;
          }

        }
      }

      // create edges and connect them to elements
      {
        BucketVector element_buckets;

        get_buckets( element_selector & mesh.mesh_meta_data().locally_owned_part(), mesh.buckets(stk::topology::ELEMENT_RANK), element_buckets);

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

        BucketVector face_buckets;

        get_buckets( element_selector & (mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part()), mesh.buckets(stk::topology::FACE_RANK), face_buckets);

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
    update_shared_edges_global_ids( mesh, shared_edge_map, mesh.m_entity_repo );
  }

  mesh.modification_end();
}

}
}
