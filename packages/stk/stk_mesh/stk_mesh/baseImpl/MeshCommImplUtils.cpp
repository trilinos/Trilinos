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

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/SameType.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/util/StaticAssert.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/MeshCommImplUtils.hpp>

#include <vector>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

bool shared_with_proc(const EntityCommListInfo& info, int proc) {
    if(info.entity_comm != nullptr)
    {    
        const EntityCommInfoVector& comm_vec = info.entity_comm->comm_map;
        for(size_t i=0; i<comm_vec.size(); ++i) {
            if (comm_vec[i].ghost_id!=BulkData::SHARED) {
                return false;
            }
            if (comm_vec[i].proc == proc) {
                return true;
            }
        }
    }    
    return false;
}

void pack_induced_memberships_for_entities_less_than_element_rank(
         const BulkData& bulk_data,
         stk::CommSparse& comm,
         const EntityCommListInfoVector & entityCommListInfoVec )
{
    OrdinalVector induced;
    for(size_t i = 0; i < entityCommListInfoVec.size(); ++i) 
    {
        stk::mesh::Entity entity = entityCommListInfoVec[i].entity;
        const int owner = bulk_data.parallel_owner_rank(entity);

        if (bulk_data.entity_rank(entity) < stk::topology::ELEM_RANK && shared_with_proc(entityCommListInfoVec[i], owner) )
        {
            const EntityState state = bulk_data.state(entity);
            if(state == stk::mesh::Modified || state == stk::mesh::Created)
            {
                induced.clear();

                induced_part_membership(bulk_data, entityCommListInfoVec[i].entity, induced);

                CommBuffer & buf = comm.send_buffer(owner);

                unsigned tmp = induced.size();

                buf.pack<stk::mesh::EntityKey>(entityCommListInfoVec[i].key);
                buf.pack<unsigned>(tmp);

                for(size_t j = 0; j < induced.size(); ++j) 
                {
                    buf.pack<unsigned>(induced[j]);
                }    
            }    
        }    
    }    
}

void append_parts_from_sharer_to_owner(const BulkData& bulk,
                                       Entity entity,
                                       OrdinalVector &partOrdinals)
{
    if(bulk.state(entity)==stk::mesh::Created)
    {    
        const stk::mesh::PartVector& all_parts = bulk.bucket(entity).supersets();
        partOrdinals.reserve(all_parts.size());
        for(stk::mesh::Part* part : all_parts)
        {
            bool isPartSameRankAsEntity = part->primary_entity_rank()==bulk.entity_rank(entity);
            bool isPartParallelConsistent = part->entity_membership_is_parallel_consistent();
            bool isRootTopologyPart = stk::mesh::is_topology_root_part(*part);
            bool isAutoDeclaredPart = stk::mesh::is_auto_declared_part(*part);
            if(isPartSameRankAsEntity && isPartParallelConsistent && !isRootTopologyPart && !isAutoDeclaredPart)
            {
                partOrdinals.push_back(part->mesh_meta_data_ordinal());
            }
        }
    }
}

void pack_induced_memberships( const BulkData& bulk_data,
                               stk::CommSparse & comm ,
                               const EntityCommListInfoVector & entityCommList )
{
  OrdinalVector induced ;
  const int myProc = bulk_data.parallel_rank();
  for ( size_t i=0; i<entityCommList.size(); ++i) {

    const int owner = bulk_data.parallel_owner_rank(entityCommList[i].entity);
    if ( owner != myProc && shared_with_proc( entityCommList[i] , owner ) ) {
      // Is shared with owner, send to owner.

      induced.clear();

      induced_part_membership(bulk_data, entityCommList[i].entity , induced );
      append_parts_from_sharer_to_owner(bulk_data, entityCommList[i].entity, induced);

      CommBuffer & buf = comm.send_buffer( bulk_data.parallel_owner_rank(entityCommList[i].entity) );

      const unsigned tmp = induced.size();

      buf.pack<stk::mesh::EntityKey>(entityCommList[i].key);
      buf.pack<unsigned>( tmp );

      for ( size_t j=0; j<induced.size(); ++j) {
        buf.pack<unsigned>( induced[j] );
      }
    }
  }
}

void pack_part_memberships( const BulkData& meshbulk, stk::CommSparse & comm ,
                            const std::vector<EntityProc> & send_list )
{
  const unsigned part_ord_shared = meshbulk.mesh_meta_data().globally_shared_part().mesh_meta_data_ordinal();
  for ( std::vector<EntityProc>::const_iterator
        i = send_list.begin() ; i != send_list.end() ; ++i ) {

    Entity entity = i->first;

    std::pair<const unsigned *, const unsigned *>
      part_ord = meshbulk.bucket(entity).superset_part_ordinals();

    // I am the owner.  The first two memberships are
    // universal_part and locally_owned_part.  The third
    // membership may be globally_shared_part ;
    //

    const unsigned count_all  = part_ord.second - part_ord.first ;
    const unsigned count_skip =
      ( 2 < count_all && part_ord.first[2] == part_ord_shared ) ? 3 : 2 ;

    const unsigned count_send = count_all - count_skip ;

    const unsigned * const start_send = part_ord.first + count_skip ;

    comm.send_buffer( i->second ).pack<EntityKey>( meshbulk.entity_key(entity) )
                                 .pack<unsigned>( count_send )
                                 .pack<unsigned>( start_send , count_send );
  }
}

void unpack_induced_parts_from_sharers(OrdinalVector& induced_parts,
                                   const EntityCommInfoVector& entity_comm_info,
                                   stk::CommSparse& comm,
                                   EntityKey expected_key)
{
//  Mesh entity membership changes must be synchronized among
//  processes that share mesh entities and propagated to
//  processes that ghost copies of the mesh entities.
//
//  Precondition: correct shared and ghosting lists.
//
//  Part memberships may have been added or removed
//  either explicitly or indirectly via entity relationships
//  being added or removed.q

    for(PairIterEntityComm ec = shared_comm_info_range(entity_comm_info); !ec.empty(); ++ec)
    {
        CommBuffer & buf = comm.recv_buffer(ec->proc);

        unsigned count = 0;
        stk::mesh::EntityKey key;
        buf.unpack<stk::mesh::EntityKey>(key);
        ThrowAssertMsg(key == expected_key, "Program error. Contact sierra-help@sandia.gov for support. Key mismatch!" << key << " not same as " << expected_key);

        buf.unpack<unsigned>(count);
        for(unsigned j = 0; j < count; ++j)
        {
            unsigned part_ord = 0;
            buf.unpack<unsigned>(part_ord);
            stk::util::insert_keep_sorted_and_unique(part_ord, induced_parts);
        }
    }
}

void pack_and_send_induced_parts_from_sharers_to_owners(const BulkData& bulkData, stk::CommSparse& comm, EntityCommListInfoVector& entity_comm_list)
{
    pack_and_communicate(comm,[&bulkData, &comm, &entity_comm_list]()
          { pack_induced_memberships(bulkData, comm, entity_comm_list); });
}

void filter_out_unneeded_induced_parts(const BulkData& bulkData, stk::mesh::Entity entity,
        const OrdinalVector& induced_parts, stk::mesh::OrdinalVector& remove_parts)
{
    remove_parts.clear();
    stk::mesh::EntityRank rank = bulkData.entity_rank(entity);
    const stk::mesh::PartVector& current_parts = bulkData.bucket(entity).supersets();

    OrdinalVector::const_iterator induced_parts_begin = induced_parts.begin(),
            induced_parts_end = induced_parts.end();

    for(stk::mesh::Part* part : current_parts)
    {
        if(part->was_induced(rank) &&
           !contains_ordinal(induced_parts_begin, induced_parts_end, part->mesh_meta_data_ordinal()))
        {
            remove_parts.push_back(part->mesh_meta_data_ordinal());
        }
    }
}

void communicate_shared_entity_info(const BulkData &mesh,
                  stk::CommSparse &comm,
                  std::vector<std::vector<shared_entity_type> > &shared_entities)
{
  stk::pack_and_communicate(comm, [&comm, &mesh, &shared_entities]() {
    const int parallelSize = mesh.parallel_size();
    for(int proc = 0; proc < parallelSize; ++proc)
    {
      if(proc != mesh.parallel_rank())
      {
        const size_t numShared = shared_entities[proc].size();
        for (size_t e = 0; e < numShared; ++e) {
          shared_entity_type const & sentity = shared_entities[proc][e];
          const size_t num_nodes_on_entity = sentity.nodes.size();
          comm.send_buffer(proc).pack<stk::topology::topology_t>(sentity.topology);
          for (size_t i = 0; i < num_nodes_on_entity; ++i ) {
            comm.send_buffer(proc).pack(sentity.nodes[i]);
          }
          comm.send_buffer(proc).pack(sentity.local_key);
        }
      }
    }
  });
}

void communicateSharingInfoToProcsThatShareEntity(const int numProcs, const int myProcId, stk::CommSparse& commStage2, stk::mesh::EntityToDependentProcessorsMap &entityKeySharing)
{
    for(int phase = 0; phase < 2; ++phase)
    {    
        stk::mesh::EntityToDependentProcessorsMap::iterator iter = entityKeySharing.begin();
        for(; iter != entityKeySharing.end(); iter++)
        {
            std::vector<int> sharingProcs(iter->second.begin(), iter->second.end());
            for(size_t j = 0; j < sharingProcs.size(); j++)
            {
                if(sharingProcs[j] == myProcId) { continue; }
                commStage2.send_buffer(sharingProcs[j]).pack<stk::mesh::EntityKey>(iter->first);
                commStage2.send_buffer(sharingProcs[j]).pack<size_t>(sharingProcs.size());
                for(size_t k = 0; k < sharingProcs.size(); k++)
                {
                    commStage2.send_buffer(sharingProcs[j]).pack<int>(sharingProcs[k]);
                }
            }
        }

        if(phase == 0)
        {
            commStage2.allocate_buffers();
        }
        else
        {
            commStage2.communicate();
        }
    }
}

void unpackCommunicationsAndStoreSharedEntityToProcPair(const int numProcs, const int myProcId, stk::CommSparse& commStage2, std::vector<std::pair<stk::mesh::EntityKey, int> >& sharedEntities)
{
    for(int procIndex = 0; procIndex < numProcs; procIndex++)
    {
        if(myProcId == procIndex) { continue; }
        stk::CommBuffer & dataFromAnotherProc = commStage2.recv_buffer(procIndex);
        while(dataFromAnotherProc.remaining())
        {
            EntityKey key;
            size_t numSharingProcs = 0;
            dataFromAnotherProc.unpack<stk::mesh::EntityKey>(key);
            dataFromAnotherProc.unpack<size_t>(numSharingProcs);
            for(size_t j = 0; j < numSharingProcs; j++)
            {
                int sharingProc = -1;
                dataFromAnotherProc.unpack<int>(sharingProc);
                if(sharingProc != myProcId)
                {
                    sharedEntities.emplace_back(key, sharingProc);
                }
            }
        }
    }
}

bool is_received_entity_in_local_shared_entity_list(
          bool use_entity_ids_for_resolving_sharing,
          const std::vector<shared_entity_type>::iterator &shared_itr,
          const std::vector<shared_entity_type>& shared_entities_this_proc,
          const shared_entity_type &shared_entity_from_other_proc)
{
  bool entitiesHaveSameNodes = shared_itr != shared_entities_this_proc.end() && *shared_itr == shared_entity_from_other_proc;
  bool entitiesAreTheSame = false;

  if ( use_entity_ids_for_resolving_sharing ) {
    entitiesAreTheSame = entitiesHaveSameNodes && shared_itr->local_key == shared_entity_from_other_proc.local_key;
  }
  else {
    entitiesAreTheSame = entitiesHaveSameNodes;
  }

  return entitiesAreTheSame;
}

bool ghost_id_is_found_in_comm_data(const PairIterEntityComm& comm_data,
                                    int entity_owner,
                                    int ghost_id)
{
  bool found_ghost_id = false;
  for (size_t i = 0; i < comm_data.size(); ++i) {
    if ((comm_data[i].ghost_id == static_cast<unsigned>(ghost_id)) &&
        (comm_data[i].proc == entity_owner)) {
          found_ghost_id = true;
          break;
    }
  }   
  return found_ghost_id;
}   

bool all_ghost_ids_are_found_in_comm_data(const PairIterEntityComm& comm_data,
                                          int entity_owner,
                                          const std::vector<int>& recvd_ghost_ids)
{
  bool found_all_ghost_ids = true;
  for (int ghost_id : recvd_ghost_ids) {
    if (!ghost_id_is_found_in_comm_data(comm_data, entity_owner, ghost_id)) {
      found_all_ghost_ids = false;
      break;
    }   
  }   
  return found_all_ghost_ids;
}

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

