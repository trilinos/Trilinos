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
#ifndef STK_TOOLS_CustomAura_hpp_
#define STK_TOOLS_CustomAura_hpp_

namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { class Ghosting; } }

#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_util/parallel/CommSparse.hpp"

namespace stk {
namespace tools {

using SymmCommData = std::pair<stk::mesh::EntityKey,int>;
using SymmCommDataVector = std::vector<SymmCommData>;
using ProcSymmCommDataVector = std::vector<SymmCommDataVector>;
using SymmCommMap = std::map<stk::mesh::EntityKey, std::set<int>>;

inline void populate_symm_comm_map(stk::mesh::BulkData& bulk,
                                   SymmCommMap& symmCommMap,
                                   stk::topology::rank_t rank) {
  symmCommMap.clear();
  int numProcs = bulk.parallel_size();
  int myProc = bulk.parallel_rank();

  // Populate Local Comm Map
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(rank, bulk.mesh_meta_data().locally_owned_part());
  for(const stk::mesh::Bucket *bucket : buckets) {
    for(stk::mesh::Entity node : *bucket) {
      std::vector<int> commProcs;
      bulk.comm_procs(node,commProcs);
      if ( commProcs.size() > 0 ) {
        commProcs.push_back(myProc);
        std::set<int> cProcs(commProcs.begin(),commProcs.end());
        symmCommMap[bulk.entity_key(node)] = cProcs;
      }
    }
  }

  // Extract
  ProcSymmCommDataVector procCommData(numProcs);
  for( auto& entry : symmCommMap ) {
    std::set<int>& procs = entry.second;
    for( int ip : procs ) {
      for ( int ij : procs ) {
        procCommData[ip].push_back(std::make_pair(entry.first,ij));
      }
    }
  }

  // Pack
  stk::CommSparse commSparse(bulk.parallel());
  for(int phase = 0; phase < 2; ++phase) {
    for (int proc=0;proc<numProcs;proc++) {
      if ( proc != myProc ) {
        SymmCommDataVector& dataVec = procCommData[proc];
        if (!dataVec.empty()) {
          stk::pack_vector_to_proc(commSparse, dataVec, proc);
        }
      }
    }
    if(phase == 0) {
      commSparse.allocate_buffers();
    } else {
      commSparse.communicate();
    }
  }

  // Unpack
  for (int proc=0;proc<numProcs;proc++) {
    if ( proc != myProc ) {
      while(commSparse.recv_buffer(proc).remaining()) {
        SymmCommDataVector commDataVector;
        stk::unpack_vector_from_proc(commSparse, commDataVector, proc);
        for( SymmCommData& a : commDataVector ) {
          symmCommMap[a.first].insert(a.second);
        }
      }
    }
  }

}

inline void nodes_in_rings(const stk::mesh::BulkData& bulk,
                           const stk::mesh::EntityVector& nodesToSearchAround,
                           stk::mesh::EntityVector& nodesInRings,
                           const int num_rings,
                           const stk::mesh::Selector& insideSelection,
                           int currentRing = 0){
  stk::mesh::EntityVector newlyAddedNodes;
  for (unsigned int iSource=0; iSource<nodesToSearchAround.size(); ++iSource){
    auto& sourceNode = nodesToSearchAround[iSource];
    unsigned int numConnectedElements = bulk.num_elements(sourceNode);
    const stk::mesh::Entity* connectedElems = bulk.begin_elements(sourceNode);
    for (unsigned int iElem=0; iElem<numConnectedElements; ++iElem){
      const stk::mesh::Entity& elem = connectedElems[iElem];
      unsigned int numConnNodes = bulk.num_nodes(elem);
      const stk::mesh::Entity* connectedNodes = bulk.begin_nodes(elem);
      for (unsigned int iNode=0; iNode<numConnNodes; ++iNode){
        stk::mesh::Entity connectedNode = connectedNodes[iNode];
        auto foundIter = std::find(nodesInRings.begin(), nodesInRings.end(), connectedNode);
        if (foundIter != nodesInRings.end()){ continue;}  // Already found
        if ( !insideSelection(bulk.bucket(connectedNode)) ) { continue; } // Not insideSelection
        newlyAddedNodes.emplace_back(connectedNode);
        nodesInRings.emplace_back(connectedNode);
      }
    }
  }

  ++currentRing;

  if ( currentRing >= num_rings || newlyAddedNodes.size() == 0 ) { return; }

  nodes_in_rings(bulk, newlyAddedNodes, nodesInRings, num_rings, insideSelection, currentRing );

}

inline std::set<int> procs_in_rings(stk::mesh::BulkData& bulkData,
                                    stk::mesh::Entity& node,
                                    const SymmCommMap& symmNodeMap,
                                    int numRings,
                                    const stk::mesh::Selector& insideSelection) {
  stk::mesh::EntityVector nodesInRings;
  nodes_in_rings(bulkData,{node},nodesInRings,numRings,insideSelection);
  std::set<int> procs;
  for(auto& ringNode: nodesInRings) {
    stk::mesh::EntityKey ringKey = bulkData.entity_key(ringNode);
    if (symmNodeMap.count(ringKey) > 0 ) {
      const std::set<int>& proc_list = symmNodeMap.at(ringKey);
      for(int p : proc_list) {
        procs.emplace(p);
      }
    } else {
      std::vector<int> proc_list;
      bulkData.comm_procs(ringNode, proc_list);
      for(int p : proc_list) {
        procs.emplace(p);
      }
    }
  }
  return procs;
}

inline void fill_list_of_entities_to_send_for_ring_aura_ghosting(stk::mesh::BulkData& bulkData,
                                                                 stk::mesh::Selector const& selector,
                                                                 stk::mesh::EntityProcVec &entitiesToGhost,
                                                                 const SymmCommMap& symmNodeMap,
                                                                 int numRings)
{
  const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::NODE_RANK, selector);
  for(const stk::mesh::Bucket *bucket : buckets)
  {
    for(stk::mesh::Entity node : *bucket)
    {
      std::set<int> proc_list = procs_in_rings(bulkData,node,symmNodeMap,numRings,selector);
      unsigned num_elements = bulkData.num_elements(node);
      const stk::mesh::Entity* elements = bulkData.begin_elements(node);
      for(unsigned i=0;i<num_elements;++i)
      {
        if(bulkData.bucket(elements[i]).owned())
        {
          for(int proc : proc_list )
          {
            entitiesToGhost.push_back(stk::mesh::EntityProc(elements[i], proc));
          }
        }
      }
    }
  }

}

inline stk::mesh::Ghosting * create_ringed_aura(stk::mesh::BulkData & stkMeshBulkData,
                                                stk::mesh::Selector const& selector,
                                                const std::string & name,
                                                const int numRings)
{
  stk::mesh::Ghosting * customAura = nullptr;

  if (stkMeshBulkData.parallel_size() > 1) {
    if ( numRings == 1 && stkMeshBulkData.is_automatic_aura_on() ) {
      // Do Nothing
    } else {
      stkMeshBulkData.modification_begin();
      customAura = &stkMeshBulkData.create_ghosting(name);
      stkMeshBulkData.modification_end();
      for (int ring(0); ring < numRings; ++ring) {
        stkMeshBulkData.modification_begin();
        stk::mesh::EntityProcVec entitiesToGhost;
        SymmCommMap symmNodeMap;
        stk::tools::populate_symm_comm_map(stkMeshBulkData,symmNodeMap,stk::topology::NODE_RANK);
        fill_list_of_entities_to_send_for_ring_aura_ghosting(stkMeshBulkData, selector, entitiesToGhost,symmNodeMap,numRings);
        stkMeshBulkData.change_ghosting(*customAura, entitiesToGhost);
        stkMeshBulkData.modification_end();
      }
    }
  }

  return customAura;
}

inline void fill_list_of_entities_to_send_for_aura_like_ghosting(stk::mesh::BulkData& bulkData,
                                                                 stk::mesh::Selector const& selector,
                                                                 stk::mesh::EntityProcVec &entitiesToGhost)
{
  const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::NODE_RANK, selector);
  for(const stk::mesh::Bucket *bucket : buckets)
  {
    for(stk::mesh::Entity shared_node : *bucket)
    {
      unsigned num_elements = bulkData.num_elements(shared_node);
      const stk::mesh::Entity* elements = bulkData.begin_elements(shared_node);
      for(unsigned i=0;i<num_elements;++i)
      {
        if(bulkData.bucket(elements[i]).owned())
        {
          std::vector<int> comm_shared_procs;
          bulkData.comm_shared_procs(shared_node, comm_shared_procs);
          for(int proc : comm_shared_procs )
          {
            entitiesToGhost.push_back(stk::mesh::EntityProc(elements[i], proc));
          }
        }
      }
    }
  }
}

inline stk::mesh::Ghosting * create_custom_aura(stk::mesh::BulkData & stkMeshBulkData,
                                                stk::mesh::Selector const& selector,
                                                const std::string & name)
{
  stk::mesh::Ghosting * customAura = nullptr;
  if (!stkMeshBulkData.is_automatic_aura_on())
  {
    stkMeshBulkData.modification_begin();

    customAura = &stkMeshBulkData.create_ghosting(name);
    stk::mesh::EntityProcVec entitiesToGhost;
    fill_list_of_entities_to_send_for_aura_like_ghosting(stkMeshBulkData, selector, entitiesToGhost);
    stkMeshBulkData.change_ghosting(*customAura, entitiesToGhost);

    stkMeshBulkData.modification_end();
  }

  return customAura;
}

inline void destroy_custom_aura(stk::mesh::BulkData & stkMeshBulkData, stk::mesh::Ghosting * customAura)
{
  if (nullptr != customAura)
  {
    stkMeshBulkData.modification_begin();
    stkMeshBulkData.destroy_ghosting(*customAura);
    stkMeshBulkData.modification_end();
  }
}

}
}

#endif
