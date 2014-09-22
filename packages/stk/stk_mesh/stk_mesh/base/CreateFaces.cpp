/*------------------------------------------------------------------------*/
/*  Copyright (c) 2013, Sandia Corporation.
/*  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*  the U.S. Governement retains certain rights in this software.
/*  
/*  Redistribution and use in source and binary forms, with or without
/*  modification, are permitted provided that the following conditions are
/*  met:
/*  
/*      * Redistributions of source code must retain the above copyright
/*        notice, this list of conditions and the following disclaimer.
/*  
/*      * Redistributions in binary form must reproduce the above
/*        copyright notice, this list of conditions and the following
/*        disclaimer in the documentation and/or other materials provided
/*        with the distribution.
/*  
/*      * Neither the name of Sandia Corporation nor the names of its
/*        contributors may be used to endorse or promote products derived
/*        from this software without specific prior written permission.
/*  
/*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*  
/*------------------------------------------------------------------------*/
/*
 * CreateFaces.cpp created on Feb 20, 2014 by tcfishe
 */

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

#include "stk_topology/apply_functor.tcc"  // for topology::apply_functor
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_topology/topology.tcc"    // for topology::num_nodes
#include "stk_topology/topology_type.tcc"  // for topology::topology_type

#include <stk_util/parallel/ParallelComm.hpp>  // for CommBuffer, CommAll
#include "stk_util/util/NamedPair.hpp"  // for EntityCommInfo::operator=, etc

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

  friend inline bool operator < (shared_face_type const& l, shared_face_type const& r)
  {
    if (l.topology < r.topology)   return true;
    if (l.topology > r.topology)   return false;

    int num_nodes = static_cast<int>(l.nodes.size());
    for (int k = 0; k < num_nodes-1; ++k) {
      if (l.nodes[k] < r.nodes[k]) return true;
      if (l.nodes[k] > r.nodes[k]) return false;
    }
    const int j = num_nodes-1;
    return l.nodes[j] < r.nodes[j];
  }

  friend inline bool operator == (shared_face_type const& l, shared_face_type const& r)
	    {
    bool equal = l.topology == r.topology;
    int num_nodes = static_cast<int>(l.nodes.size());
    for (int k = 0; k < num_nodes; ++k) {
      bool node_found = false;
      for (int j = 0; j < num_nodes; ++j) {
        node_found = node_found || (l.nodes[k] == r.nodes[j]);
      }
      equal = equal && node_found;
    }
    return equal;
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

  create_face_impl(   size_t               & next_face
      , face_map_type        & face_map
      , shared_face_map_type & shared_face_map
      , Bucket               & bucket
  )
  : m_next_face(next_face)
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
        for (int f=0 ; f < num_existing_faces ; ++f) {
          face_exists[face_ords[f]] = mesh.is_valid(face_entity[f]);
        }
      }

      for (unsigned f=0; f < Topology::num_faces; ++f) {
        if (face_exists[f]) {
          continue;
        }

        topology faceTopology = elemTopology.face_topology(f);

        // Use node identifier instead of node local_offset for cross-processor consistency.
        EntityIdVector face_node_ids(faceTopology.num_nodes());
        Topology::face_nodes(elem_node_ids, f, face_node_ids.begin());
        const unsigned smallest_permutation = faceTopology.lexicographical_smallest_permutation(face_node_ids);

        EntityVector face_nodes(faceTopology.num_nodes());
        Topology::face_nodes(elem_nodes, f, face_nodes.begin());

        EntityVector permuted_face_nodes(faceTopology.num_nodes());
        faceTopology.permutation_nodes(face_nodes, smallest_permutation, permuted_face_nodes.begin());

        Entity face;

        typename face_map_type::iterator iface = m_face_map.find(permuted_face_nodes);
        if (iface == m_face_map.end()) {
          EntityId face_id = m_next_face++;
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
        mesh.declare_relation(m_bucket[ielem], face, f);
      }
    }
  }

  template <typename Topology>
  typename boost::enable_if_c< (Topology::num_faces == 0u), void>::type
  operator()(Topology t)
  {}


  //members
  size_t                & m_next_face;
  face_map_type         & m_face_map;
  shared_face_map_type  & m_shared_face_map;
  Bucket                & m_bucket;
};

void update_shared_faces_global_ids( BulkData & mesh, shared_face_map_type & shared_face_map)
{
  //sort the faces for ease of lookup
  std::sort(shared_face_map.begin(), shared_face_map.end());

  std::vector< shared_face_map_type >  shared_faces( mesh.parallel_size());

  for (shared_face_map_type::const_iterator itr = shared_face_map.begin(),
      end = shared_face_map.end(); itr != end; ++itr ) {
    //find process that share this face;

    PairIterEntityComm left_shared = mesh.entity_comm_map_shared(itr->nodes[0]);

    EntityCommInfoVector shared_processes(left_shared.first, left_shared.second);

    for (int i = 1; i < static_cast<int>(itr->nodes.size()); ++i) {
      PairIterEntityComm right_shared = mesh.entity_comm_map_shared( itr->nodes[i] );

      std::set_intersection( shared_processes.begin(), shared_processes.end(),
          right_shared.first, right_shared.second,
          shared_processes.begin()
      );

    }

    for (EntityCommInfoVector::const_iterator comm_itr = shared_processes.begin(),
        comm_end = shared_processes.end(); comm_itr != comm_end; ++comm_itr ) {
      if (comm_itr->proc > mesh.parallel_rank()) {
        shared_faces[comm_itr->proc].push_back(*itr);
      }
    }
  }

  CommAll comm( mesh.parallel() );

  //pack send buffers
  for (int allocation_pass=0; allocation_pass<2; ++allocation_pass) {
    if (allocation_pass >= 1) {
      comm.allocate_buffers( mesh.parallel_size() /4, 0);
    }

    for (int proc=mesh.parallel_rank()+1, parallel_size = mesh.parallel_size(); proc<parallel_size; ++proc) {
      for (size_t e=0, num_shared = shared_faces[proc].size(); e < num_shared; ++e) {
        shared_face_type const & sface = shared_faces[proc][e];
        CommBuffer & commbuff = comm.send_buffer(proc).pack<stk::topology::topology_t>(sface.topology);
        for (int i = 0; i < static_cast<int>(sface.nodes.size()); ++i) {
          commbuff.pack<EntityKey>(sface.nodes[i]);
        }
        commbuff.pack<EntityKey>(sface.local_key);
      }
    }
  }

  comm.communicate();

  for ( unsigned ip = mesh.parallel_rank() ; ip > 0 ; ) {
    --ip;
    CommBuffer & buf = comm.recv_buffer( ip );
    while ( buf.remaining() ) {
      stk::topology::topology_t topology;
      buf.unpack<stk::topology::topology_t>(topology);
      shared_face_type sface(topology);
      for (int i = 0; i < static_cast<int>(sface.nodes.size()); ++i) {
        buf.unpack<EntityKey>(sface.nodes[i]);
      }

      buf.unpack<EntityKey>(sface.global_key);

      shared_face_map_type::iterator shared_itr = std::lower_bound(shared_face_map.begin(), shared_face_map.end(), sface);

      //update the global global_key
      if (shared_itr != shared_face_map.end() && *shared_itr == sface) {
        shared_itr->global_key = sface.global_key;
      }
    }
  }

  //update the entity keys
  for (size_t i=0, e=shared_face_map.size(); i<e; ++i) {
    if (shared_face_map[i].global_key != shared_face_map[i].local_key) {
      mesh.internal_change_entity_key(shared_face_map[i].local_key, shared_face_map[i].global_key, mesh.get_entity(shared_face_map[i].local_key));
    }
  }
}

} //namespace

void create_faces( BulkData & mesh )
{
  create_faces(mesh, mesh.mesh_meta_data().universal_part());
}

void create_faces( BulkData & mesh, const Selector & element_selector )
{
  // static size_t next_face = static_cast<size_t>(mesh.parallel_rank()+1) << 32;
  // NOTE: This is a workaround to eliminate some bad behavior with the equation above when
  //       the #proc is a power of two.  The 256 below is the bin size of the Distributed Index.
  static size_t next_face = (static_cast<size_t>(mesh.parallel_rank()+1) << 32) + 256 * mesh.parallel_rank();

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

  //create the faces for the elements in each bucket
  for (size_t i=0, e=element_buckets.size(); i<e; ++i) {
    Bucket &b = *element_buckets[i];

    create_face_impl functor( next_face, face_map, shared_face_map, b);
    stk::topology::apply_functor< create_face_impl > apply(functor);
    apply( b.topology() );
  }

  //update global ids of shared faces
  update_shared_faces_global_ids( mesh, shared_face_map );

  if (i_started) {
    mesh.modification_end_for_entity_creation( stk::topology::FACE_RANK, BulkData::MOD_END_COMPRESS_AND_SORT );
  }
}

}
}
