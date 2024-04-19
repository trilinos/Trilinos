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

#ifndef STK_IO_ProcessSetsOrBlocks_h
#define STK_IO_ProcessSetsOrBlocks_h

#include "Ioss_Region.h"                // for Region, NodeSetContainer, etc
#include "StkMeshIoBroker.hpp"
#include "IossBridge.hpp"
#include <stk_mesh/base/Types.hpp>
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/CreateEdges.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include "Ioss_EdgeBlock.h"
#include "Ioss_EdgeSet.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementSet.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_EntityBlock.h"
#include "Ioss_EntityType.h"
#include "Ioss_FaceBlock.h"
#include "Ioss_FaceSet.h"
#include "Ioss_Field.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_Property.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_CommSet.h"

#include <map>

namespace stk { namespace io {

template <typename INT>
void process_node_sharing(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  // Register node sharing information for all nodes on processor
  // boundaries.
  //

  stk::ParallelMachine myComm = bulk.parallel();
  int numProc = bulk.parallel_size();
  if (bulk.parallel_size() > 1)
  {
    Ioss::CommSet* io_cs = region.get_commset("commset_node");
    size_t num_sharings = io_cs->get_field("entity_processor").raw_count();

    // Check for corrupt incomplete nemesis information.  Some old
    // files are being used which do not have the correct nemesis
    // sharing data. They can be identified by an incorrect global
    // node count (typically equal to 1) in addition to an empty node sharing list.
    // Assume that if the node sharing list is non-empty, then no matter  what the
    // global node count is, the data is most likely ok.
    size_t global_node_count = region.get_property("global_node_count").get_int();
    size_t local_node_count = stk::mesh::count_entities(bulk, stk::topology::NODE_RANK, bulk.mesh_meta_data().locally_owned_part());

    STK_ThrowErrorMsgIf(num_sharings == 0 && global_node_count < local_node_count,
                    "ERROR: Invalid communication/node sharing information found in file '"
                     << region.get_database()->get_filename() << "'\n"
                     << "       There is no node sharing information and the "
                     << "global node count is  " << global_node_count
                     << " which is less than the node count on processor "
                     << stk::parallel_machine_rank(bulk.parallel())
                     << " which is " << local_node_count << ".\n"
                     << "       A possible work-around is to join (epu) and re-spread (decomp) the mesh files.");

    std::vector<INT> entity_proc;
    io_cs->get_field_data("entity_processor", entity_proc);

    std::set< std::pair<INT, int> > disconnectedNodes;

    for (size_t i = 0; i < num_sharings; ++i) {
        INT nodeId = entity_proc[i*2];
        stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, nodeId);
        if(!bulk.is_valid(node)) {
            int otherProc = entity_proc[i*2+1];
            disconnectedNodes.insert(std::make_pair(nodeId, otherProc));
        }
    }

    stk::CommSparse commSparse(myComm);
    stk::pack_and_communicate(commSparse, [&]()
                              {
                                 for (auto &nodeAndProc : disconnectedNodes) {
                                     INT nodeId = nodeAndProc.first;
                                     int otherProc = nodeAndProc.second;
                                     commSparse.send_buffer(otherProc).pack(nodeId);
                                 }
                              });

    for(int i=0; i<numProc; ++i) {
        stk::CommBuffer &buf = commSparse.recv_buffer(i);
        while(buf.remaining()) {
            INT nodeId;
            buf.unpack(nodeId);
            disconnectedNodes.insert(std::make_pair(nodeId, i));
        }
    }

    for (size_t i = 0; i < num_sharings; ++i) {
        INT nodeId = entity_proc[i*2];
        int otherProc = entity_proc[i*2+1];
        std::pair<INT, int> nodeAndProc = std::make_pair(nodeId, otherProc);
        if(disconnectedNodes.find(nodeAndProc) == disconnectedNodes.end()) {
            stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, nodeId);
            bulk.add_node_sharing(node, otherProc);
        }
    }
  }
}


void process_nodeblocks(Ioss::Region &region, stk::mesh::MetaData &meta);

template <typename INT>
void process_nodeblocks(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  // This must be called after the "process_element_blocks" call
  // since there may be nodes that exist in the database that are
  // not part of the analysis mesh due to subsetting of the element
  // blocks.

  // Currently, all nodes found in the finite element mesh are defined
  // as nodes in the stk_mesh database. If some of the element blocks
  // are omitted, then there will be disconnected nodes defined.
  // However, if we only define nodes that are connected to elements,
  // then we risk missing "free" nodes that the user may want to have
  // existing in the model.
  const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
  assert(node_blocks.size() == 1);

  Ioss::NodeBlock *nb = node_blocks[0];

  std::vector<INT> ids;
  nb->get_field_data("ids", ids);

  stk::mesh::Part& nodePart = bulk.mesh_meta_data().get_topology_root_part(stk::topology::NODE);
  stk::mesh::PartVector nodeParts = {&nodePart};

  std::vector<stk::mesh::Entity> nodes;
  nodes.reserve(ids.size());

  bulk.declare_entities(stk::topology::NODE_RANK, ids, nodeParts, nodes);

  for (size_t i=0; i < ids.size(); i++) {
    bulk.set_local_id(nodes[i], i);
  }

  process_node_sharing<INT>(region, bulk);
}

stk::mesh::Part* get_part_from_alias(const Ioss::Region &region, const stk::mesh::MetaData &meta, const std::string &name);
stk::mesh::Part* get_part_for_grouping_entity(const Ioss::Region &region, const stk::mesh::MetaData &meta, const Ioss::GroupingEntity *entity);

void process_elementblocks(Ioss::Region &region, stk::mesh::MetaData &meta, TopologyErrorHandler handler);
template <typename INT>
void process_elementblocks(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::mesh::EntityVector elems;
  stk::mesh::Permutation perm = stk::mesh::Permutation::INVALID_PERMUTATION;
  stk::mesh::OrdinalVector scratch1, scratch2, scratch3;

  stk::mesh::Part& nodePart = bulk.mesh_meta_data().get_topology_root_part(stk::topology::NODE);
  stk::mesh::PartVector nodeParts = {&nodePart};

  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
  for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
      it != elem_blocks.end(); ++it) {
    Ioss::ElementBlock *entity = *it;

    if (stk::io::include_entity(entity)) {
      stk::mesh::Part* part = get_part_for_grouping_entity(region, meta, entity);
      stk::mesh::PartVector elemParts = {part};
      if (part != nullptr) {
          std::vector<INT> elem_ids ;
          std::vector<INT> connectivity ;

          entity->get_field_data("ids", elem_ids);
          entity->get_field_data("connectivity", connectivity);

          if (elem_ids.empty()) {
            continue;
          }

          stk::topology topo = part->topology();
          STK_ThrowRequireMsg( topo != stk::topology::INVALID_TOPOLOGY, " INTERNAL_ERROR: Part " << part->name() << " has invalid topology");

          bulk.declare_entities(stk::topology::ELEM_RANK, elem_ids, elemParts, elems);

          size_t element_count = elem_ids.size();
          int nodes_per_elem = topo.num_nodes();

          size_t offset = entity->get_offset();
          for(size_t i=0; i<element_count; ++i) {
              INT *conn = &connectivity[i*nodes_per_elem];
              stk::mesh::Entity element = elems[i];

              bulk.set_local_id(element, offset + i);

              for(int j = 0; j < nodes_per_elem; ++j)
              {
                  stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, conn[j]);
                  if (!bulk.is_valid(node)) {
                      node = bulk.declare_node(conn[j], nodeParts);
                  }
                  bulk.declare_relation(element, node, j, perm, scratch1, scratch2, scratch3);
              }
          }
      }
    }
  }
}

void process_nodesets_without_distribution_factors(Ioss::Region &region, stk::mesh::MetaData &meta);
void process_nodesets(Ioss::Region &region, stk::mesh::MetaData &meta);
template <typename INT>
void process_nodesets(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
    // Should only process nodes that have already been defined via the element
    // blocks connectivity lists.
    const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
    const stk::mesh::MetaData &meta = bulk.mesh_meta_data();

    for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
            it != node_sets.end(); ++it) {
        Ioss::NodeSet *entity = *it;

        if (stk::io::include_entity(entity)) {
            stk::mesh::Part* part = get_part_for_grouping_entity(region, meta, entity);
            if(nullptr != part)
            {
                stk::mesh::PartVector add_parts( 1 , part );

                std::vector<INT> node_ids ;
                size_t node_count = entity->get_field_data("ids", node_ids);

                for(size_t i=0; i<node_count; ++i) {
                    stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, node_ids[i] );
                    if (bulk.is_valid(node)) {
                        bulk.change_entity_parts(node, add_parts);
                    }
                }
            }
        }
    }
}

struct NodesetCompare
{
    bool operator()(const Ioss::NodeSet* nodeset1, const Ioss::NodeSet* nodeset2) const
    {
        return nodeset1->hash() < nodeset2->hash();
    }
};
typedef std::map<Ioss::NodeSet*, stk::mesh::Part*, NodesetCompare> NodesetMap;
void populate_hidden_nodesets(Ioss::Region &io, const stk::mesh::MetaData & meta, NodesetMap &nodesetMap);

template <typename INT>
void process_hidden_nodesets(Ioss::Region &io, stk::mesh::BulkData & bulk)
{
    stk::mesh::MetaData &meta = bulk.mesh_meta_data();
    NodesetMap nodesetMap;

    populate_hidden_nodesets(io, meta, nodesetMap);

    for(auto pos = nodesetMap.begin(); pos != nodesetMap.end(); pos++)
    {
        Ioss::NodeSet* io_node_set = pos->first;
        stk::mesh::Part *mesh_part = pos->second;

        stk::mesh::PartVector add_parts( 1 , mesh_part );

        std::vector<INT> node_ids ;
        size_t node_count = io_node_set->get_field_data("ids", node_ids);

        for(size_t i=0; i<node_count; ++i) {
            stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, node_ids[i] );
            if (bulk.is_valid(node)) {
                bulk.change_entity_parts(node, add_parts);
            }
        }
    }
}

void process_sidesets(Ioss::Region &region, stk::mesh::BulkData &bulk, const stk::mesh::EntityIdProcMap &elemIdMovedToProc, stk::io::StkMeshIoBroker::SideSetFaceCreationBehavior behavior);
void process_sidesets(Ioss::Region &region, stk::mesh::MetaData &meta);
void process_face_blocks(Ioss::Region &region, stk::mesh::BulkData &bulk);
void process_face_blocks(Ioss::Region &region, stk::mesh::MetaData &meta, TopologyErrorHandler handler);
void process_edge_blocks(Ioss::Region &region, stk::mesh::BulkData &bulk);
void process_edge_blocks(Ioss::Region &region, stk::mesh::MetaData &meta, TopologyErrorHandler handler);
void process_assemblies(Ioss::Region &region, stk::mesh::MetaData &meta);
void build_assembly_hierarchies(Ioss::Region &region, stk::mesh::MetaData &meta);

}} // namespace stk io

#endif
