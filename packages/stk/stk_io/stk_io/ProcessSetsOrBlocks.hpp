#ifndef STK_IO_ProcessSetsOrBlocks_h
#define STK_IO_ProcessSetsOrBlocks_h

#include "Ioss_Region.h"                // for Region, NodeSetContainer, etc
#include "StkMeshIoBroker.hpp"
#include "Ioss_GroupingEntity.h"
#include "IossBridge.hpp"
#include <stk_mesh/base/Types.hpp>
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/BulkData.hpp"

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
#include "stk_mesh/base/FEMHelpers.hpp"

namespace stk { namespace io {

stk::mesh::EntityId get_side_entity_id(int64_t elem_id, int side_ordinal);

void process_nodeblocks(Ioss::Region &region, stk::mesh::MetaData &meta);

template <typename INT>
#ifdef STK_BUILT_IN_SIERRA
void process_nodeblocks(Ioss::Region &region, stk::mesh::BulkData &bulk)
#else
void process_nodeblocks(Ioss::Region &region, stk::mesh::BulkData &bulk, stk::ParallelMachine comm)
#endif
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

  stk::mesh::Part& nodePart = bulk.mesh_meta_data().get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::NODE));

  for (size_t i=0; i < ids.size(); i++) {
    stk::mesh::Entity node = bulk.declare_entity(stk::topology::NODE_RANK, ids[i], nodePart);
    bulk.set_local_id(node, i);
  }

  // Register node sharing information for all nodes on processor
  // boundaries.
  //
#ifdef STK_BUILT_IN_SIERRA
  if (bulk.parallel_size() > 1)
#else
  if (stk::parallel_machine_size(comm) > 1)
#endif
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
    ThrowErrorMsgIf (num_sharings == 0 && global_node_count < ids.size(),
                    "ERROR: Invalid communication/node sharing information found in file '"
                     << region.get_database()->get_filename() << "'\n"
                     << "       There is no node sharing information and the "
                     << "global node count is  " << global_node_count
                     << " which is less than the node count on processor "
                     << stk::parallel_machine_rank(bulk.parallel())
                     << " which is " << ids.size() << ".\n"
                     << "       A possible work-around is to join (epu) and re-spread (decomp) the mesh files.");

    std::vector<INT> entity_proc;
    io_cs->get_field_data("entity_processor", entity_proc);

    for (size_t i = 0; i < num_sharings; ++i) {
      stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, entity_proc[i*2]);
      bulk.add_node_sharing(node, entity_proc[i*2+1]);
    }
  }
}

void process_elementblocks(Ioss::Region &region, stk::mesh::MetaData &meta);
template <typename INT>
void process_elementblocks(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);

  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
  for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
      it != elem_blocks.end(); ++it) {
    Ioss::ElementBlock *entity = *it;

    if (stk::io::include_entity(entity)) {
      const std::string &name = entity->name();
      stk::mesh::Part* const part = meta.get_part(name);
      STKIORequire(part != NULL);

      stk::topology topo = part->topology();
      if (topo == stk::topology::INVALID_TOPOLOGY) {
        std::ostringstream msg ;
        msg << " INTERNAL_ERROR: Part " << part->name() << " has invalid topology";
        throw std::runtime_error( msg.str() );
      }

      std::vector<INT> elem_ids ;
      std::vector<INT> connectivity ;

      entity->get_field_data("ids", elem_ids);
      entity->get_field_data("connectivity", connectivity);

      size_t element_count = elem_ids.size();
      int nodes_per_elem = topo.num_nodes();

      stk::mesh::EntityIdVector id_vec(nodes_per_elem);

      size_t offset = entity->get_offset();
      for(size_t i=0; i<element_count; ++i) {
        INT *conn = &connectivity[i*nodes_per_elem];
        std::copy(&conn[0], &conn[0+nodes_per_elem], id_vec.begin());
        stk::mesh::Entity element = stk::mesh::declare_element(bulk, *part, elem_ids[i], id_vec);

        bulk.set_local_id(element, offset + i);
      }
    }
  }
}

void process_nodesets(Ioss::Region &region, stk::mesh::MetaData &meta);
template <typename INT>
void process_nodesets(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  // Should only process nodes that have already been defined via the element
  // blocks connectivity lists.
  const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
  const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(bulk);

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
      it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      const std::string & name = entity->name();
      stk::mesh::Part* const part = meta.get_part(name);
      STKIORequire(part != NULL);
      stk::mesh::PartVector add_parts( 1 , part );

      std::vector<INT> node_ids ;
      size_t node_count = entity->get_field_data("ids", node_ids);

      stk::mesh::EntityRank n_rank = stk::topology::NODE_RANK;
      for(size_t i=0; i<node_count; ++i) {
        stk::mesh::Entity node = bulk.get_entity(n_rank, node_ids[i] );
        if (!bulk.is_valid(node)) {
          node = bulk.declare_entity(n_rank, node_ids[i], add_parts );
        }
        else {
          bulk.change_entity_parts(node, add_parts);
        }
      }
    }
  }
}

void process_sidesets(Ioss::Region &region, stk::mesh::BulkData &bulk, const stk::mesh::EntityIdProcMap &elemIdMovedToProc, stk::io::StkMeshIoBroker::SideSetFaceCreationBehavior behavior);
void process_sidesets(Ioss::Region &region, stk::mesh::MetaData &meta);

}} // namespace stk io

#endif
