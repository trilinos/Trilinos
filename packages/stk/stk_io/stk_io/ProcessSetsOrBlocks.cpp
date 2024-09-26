
// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "ProcessSetsOrBlocks.hpp"
#include <cstdint>                                 // for int64_t
#include <map>                                     // for allocator, map<>::...
#include <stk_mesh/base/BulkData.hpp>              // for BulkData
#include <stk_mesh/base/MetaData.hpp>              // for MetaData, put_fiel...
#include <utility>                                 // for pair
#include "IossBridge.hpp"                          // for default_part_proce...
#include "Ioss_Assembly.h"                         // for Assembly, EntityCo...
#include "Ioss_EdgeBlock.h"                        // for EdgeBlock
#include "Ioss_ElementTopology.h"                  // for ElementTopology
#include "Ioss_EntityType.h"                       // for SIDESET, EDGEBLOCK
#include "Ioss_FaceBlock.h"                        // for FaceBlock
#include "Ioss_Field.h"                            // for Field, Field::MESH
#include "Ioss_GroupingEntity.h"                   // for GroupingEntity
#include "Ioss_NodeBlock.h"                        // for NodeBlock
#include "Ioss_SideBlock.h"                        // for SideBlock
#include "Ioss_SideSet.h"                          // for SideSet, SideBlock...
#include "StkIoUtils.hpp"                          // for part_primary_entit...
#include "StkMeshIoBroker.hpp"                     // for StkMeshIoBroker
#include "stk_mesh/base/Bucket.hpp"                // for Bucket
#include "stk_mesh/base/EntityKey.hpp"             // for EntityKey
#include "stk_mesh/base/FEMHelpers.hpp"            // for declare_element_edge
#include "stk_mesh/base/Field.hpp"                 // for Field
#include "stk_mesh/base/SideSetEntry.hpp"          // for SideSet
#include "stk_mesh/base/Types.hpp"                 // for EntityId, PartVector
#include "stk_mesh/baseImpl/ConnectEdgesImpl.hpp"  // for connect_face_to_edges
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"     // for connect_edge_to_el...
#include "stk_topology/topology.hpp"               // for topology, topology...
#include "stk_util/diag/StringUtil.hpp"            // for case_strcmp
#include "stk_util/environment/RuntimeWarning.hpp"
#include "stk_util/parallel/CommSparse.hpp"        // for CommSparse, pack_a...
#include "stk_util/parallel/ParallelComm.hpp"      // for CommBuffer
#include "stk_util/util/ReportHandler.hpp"         // for ThrowRequireMsg
#include "stk_util/util/SortAndUnique.hpp"         // for sort_and_unique

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk { namespace io {

struct ElemSidePartOrds
{
    stk::mesh::EntityId elem;
    int sideOrdinal;
    stk::mesh::OrdinalVector partOrdinals;
};

void process_nodeblocks(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const  Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
  assert(node_blocks.size() == 1);

  auto & coord_field = meta.declare_field<double>(stk::topology::NODE_RANK, meta.coordinate_field_name());
  stk::mesh::put_field_on_mesh(coord_field, meta.universal_part(), meta.spatial_dimension(), nullptr);
  stk::io::set_field_output_type(coord_field, stk::io::FieldOutputType::VECTOR_3D);
  stk::io::set_field_role(coord_field, Ioss::Field::MESH);
  meta.set_coordinate_field(&coord_field);

  Ioss::NodeBlock *nb = node_blocks[0];
  stk::io::define_io_fields(nb, Ioss::Field::ATTRIBUTE, meta.universal_part(), stk::topology::NODE_RANK);
}

void process_elementblocks(Ioss::Region &region, stk::mesh::MetaData &meta, TopologyErrorHandler handler)
{
  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
  stk::io::default_part_processing(elem_blocks, meta, handler);
}

void process_nodesets_without_distribution_factors(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
  stk::io::default_part_processing(node_sets, meta);
}

void process_nodesets(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
  stk::io::default_part_processing(node_sets, meta);

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
      it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      stk::mesh::Part* const part = meta.get_part(entity->name());

      STKIORequire(part !=  nullptr);
      STKIORequire(entity->field_exists("distribution_factors"));

      std::string nodesetName = part->name();
      std::string nodesetDistFieldName = "distribution_factors_" + nodesetName;

      stk::mesh::Field<double> & distribution_factors_field_per_nodeset =
        meta.declare_field<double>(stk::topology::NODE_RANK, nodesetDistFieldName);

      stk::io::set_field_role(distribution_factors_field_per_nodeset, Ioss::Field::MESH);
      stk::mesh::put_field_on_mesh(distribution_factors_field_per_nodeset, *part, nullptr);
    }
  }
}

// ========================================================================
void process_surface_entity(Ioss::SideSet *sset, stk::mesh::MetaData &meta)
{
  assert(sset->type() == Ioss::SIDESET);
  const Ioss::SideBlockContainer& blocks = sset->get_side_blocks();
  stk::io::default_part_processing(blocks, meta);
  stk::mesh::Part* const ss_part = meta.get_part(sset->name());
  STKIORequire(ss_part !=  nullptr);

  stk::mesh::FieldBase *distribution_factors_field = nullptr;
  bool surface_df_defined = false; // Has the surface df field been defined yet?

  size_t block_count = sset->block_count();
  for (size_t i=0; i < block_count; i++) {
    Ioss::SideBlock *sb = sset->get_block(i);
    if (stk::io::include_entity(sb)) {
      stk::mesh::Part * const sb_part = meta.get_part(sb->name());
      STKIORequire(sb_part != nullptr);
      if (ss_part->mesh_meta_data_ordinal() != sb_part->mesh_meta_data_ordinal()) {
        meta.declare_part_subset(*ss_part, *sb_part);
      }

      if (sb->field_exists("distribution_factors")) {
        if (!surface_df_defined) {
          stk::topology::rank_t side_rank = static_cast<stk::topology::rank_t>(stk::io::part_primary_entity_rank(*sb_part));
          std::string field_name = sset->name() + "_df";
          distribution_factors_field = &meta.declare_field<double>(side_rank, field_name);
          stk::io::set_field_role(*distribution_factors_field, Ioss::Field::MESH);
          stk::io::set_distribution_factor_field(*ss_part, *distribution_factors_field);
          surface_df_defined = true;
        }
        stk::io::set_distribution_factor_field(*sb_part, *distribution_factors_field);
        int side_node_count = sb->topology()->number_nodes();
        stk::mesh::put_field_on_mesh(*distribution_factors_field, *sb_part, side_node_count, nullptr);
      }
    }
  }
}

void process_sidesets(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::SideSetContainer& side_sets = region.get_sidesets();
  stk::io::default_part_processing(side_sets, meta);

  for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
      it != side_sets.end(); ++it) {
    Ioss::SideSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      process_surface_entity(entity, meta);
    }
  }
}

void create_processed_edge(stk::mesh::BulkData& bulk,
                           const stk::mesh::Entity elem,
                           int side_ordinal,
                           const stk::mesh::PartVector& add_parts,
                           stk::io::StkMeshIoBroker::SideSetFaceCreationBehavior behavior,
                           stk::mesh::EntityId side_id_for_classic_behavior)
{
    if(bulk.mesh_meta_data().spatial_dimension() == 2 &&
       behavior == stk::io::StkMeshIoBroker::STK_IO_SIDE_CREATION_USING_GRAPH_TEST)
    {
        bulk.declare_element_side(elem, side_ordinal, add_parts);
    }
    else
    {
        stk::mesh::Entity side = stk::mesh::declare_element_edge(bulk, side_id_for_classic_behavior, elem, side_ordinal);
        bulk.change_entity_parts(side, add_parts);
    }
}

void create_processed_face(stk::mesh::BulkData& bulk,
                          const stk::mesh::Entity elem,
                          int side_ordinal,
                          const stk::mesh::PartVector& add_parts,
                          stk::io::StkMeshIoBroker::SideSetFaceCreationBehavior behavior,
                          stk::mesh::EntityId side_id_for_classic_behavior,
                          std::vector<ElemSidePartOrds>& sidesToMove)
{
    if(behavior == stk::io::StkMeshIoBroker::STK_IO_SIDESET_FACE_CREATION_CLASSIC) {
        bulk.declare_element_side(elem, side_ordinal, add_parts);
    }
    else if(behavior == stk::io::StkMeshIoBroker::STK_IO_SIDESET_FACE_CREATION_CURRENT) {
        stk::mesh::Entity new_face =
                stk::mesh::impl::get_or_create_face_at_element_side(bulk, elem, side_ordinal,
                                                                    side_id_for_classic_behavior,
                                                                    add_parts);
        stk::mesh::impl::connect_face_to_other_elements(bulk, new_face, elem, side_ordinal);
    }
    else if(behavior == stk::io::StkMeshIoBroker::STK_IO_SIDE_CREATION_USING_GRAPH_TEST) {
        if(bulk.bucket(elem).owned()) {
            bulk.declare_element_side(elem, side_ordinal, add_parts);
        }
        else {
            stk::mesh::OrdinalVector ords(add_parts.size());
            for(size_t n = 0; n < add_parts.size(); ++n)
            {
                ords[n] = add_parts[n]->mesh_meta_data_ordinal();
            }
            sidesToMove.push_back( {bulk.identifier(elem), side_ordinal, ords});
        }
    }
}

Ioss::ElementTopology* get_side_block_topology_from_entries(Ioss::DatabaseIO* db, Ioss::SideBlock* sb)
{
  Ioss::Region* region = db->get_region();
  if(nullptr == region) return nullptr;

  Ioss::Int64Vector elementSide;
  if (db->int_byte_size_api() == 4) {
    Ioss::IntVector es32;
    sb->get_field_data("element_side", es32);
    elementSide.resize(es32.size());
    std::copy(es32.begin(), es32.end(), elementSide.begin());
  }
  else {
    sb->get_field_data("element_side", elementSide);
  }

  int heterogenousTopo = 0;
  Ioss::ElementTopology* blockSideTopo = nullptr;
  size_t              number_sides = elementSide.size() / 2;
  Ioss::ElementBlock *block        = nullptr;
  std::int64_t sideSetOffset = Ioss::Utils::get_side_offset(sb);
  for (size_t iel = 0; iel < number_sides; iel++) {
    int64_t elemId   = elementSide[2 * iel]; 
    int64_t elemSide = elementSide[2 * iel + 1] + sideSetOffset - 1;
    elemId           = db->element_global_to_local(elemId);
    if (block == nullptr || !block->contains(elemId)) {
      block = region->get_element_block(elemId);
      assert(block != nullptr);
    }

    int64_t oneBasedElemSide = elemSide + 1;
    Ioss::ElementTopology* sideTopo = block->topology()->boundary_type(oneBasedElemSide);

    if(nullptr != blockSideTopo && sideTopo != blockSideTopo) {
      blockSideTopo = nullptr;
      heterogenousTopo = 1;
      break;      
    }

    blockSideTopo = sideTopo;	
  }

  int topoId = (blockSideTopo != nullptr) ? Ioss::ElementTopology::get_unique_id(blockSideTopo->name()) : 0;

  if (db->is_parallel()) {
    std::vector<int> topoVec{topoId, heterogenousTopo};
    db->util().global_array_minmax(topoVec, Ioss::ParallelUtils::DO_MAX);
    topoId = topoVec[0];
    heterogenousTopo = topoVec[1];
  }

  blockSideTopo = heterogenousTopo ? nullptr : Ioss::ElementTopology::factory(topoId);

  return blockSideTopo;
}

bool set_sideset_topology(const Ioss::SideSet *ss, stk::mesh::Part* part,
          	          const Ioss::ElementTopology* sbTopo, stk::topology stkTopology,
			  bool printWarning = false)
{
  if (stkTopology == stk::topology::INVALID_TOPOLOGY) {
    if (printWarning) {
      std::ostringstream os;
      os<<"stk_io WARNING: failed to obtain sensible topology for sideset: " << ss->name()<<", iossTopology: "<<sbTopo->name()<<", stk-part: "<<part->name()<<", rank: "<<part->primary_entity_rank()<<", stk-topology: "<<stkTopology<<". Probably because this GroupingEntity is empty on this MPI rank. Unable to set correct stk topology hierarchy. Proceeding, but beware of unexpected behavior."<<std::endl;
      std::cerr<<os.str();
    }
  }
  else {
    for(auto subsetPart : part->subsets()) {
      if(subsetPart->topology() != stkTopology) {
	return false;
      }
    }
      
    stk::mesh::set_topology(*part, stkTopology);
    return true;
  }

  return false;
}

void set_sideset_topology(const Ioss::SideSet *ss, stk::mesh::Part *part, const stk::mesh::MetaData &meta)
{
  if(nullptr == ss) return;
  if(nullptr == part) return;
  if(ss->side_block_count() != 1) return;

  Ioss::DatabaseIO* db = ss->get_database();
  Ioss::Region* region = db->get_region();
  if(nullptr == region) return;

  Ioss::SideBlock* sb = ss->get_block(0);

  if(sb->name() != ss->name()) {
    stk::topology stkTopology = map_ioss_topology_to_stk(sb->topology(), meta.spatial_dimension());
    if(set_sideset_topology(ss, part, sb->topology(), stkTopology, true)) return;
  }

  Ioss::ElementTopology* blockSideTopo = get_side_block_topology_from_entries(db, sb);
  
  if(nullptr != blockSideTopo) {
    stk::topology stkTopology = map_ioss_topology_to_stk(blockSideTopo, meta.spatial_dimension());
    set_sideset_topology(ss, part, blockSideTopo, stkTopology);
  }
}

template <typename INT>
void process_surface_entity(const Ioss::SideSet* sset, stk::mesh::BulkData & bulk, std::vector<ElemSidePartOrds> &sidesToMove, stk::io::StkMeshIoBroker::SideSetFaceCreationBehavior behavior)
{
    assert(sset->type() == Ioss::SIDESET);

    const stk::mesh::MetaData &meta = bulk.mesh_meta_data();

    Ioss::Region *region = sset->get_database()->get_region();

    stk::mesh::SideSet *stkSideSet = nullptr;
    stk::mesh::Part *stkSideSetPart = get_part_for_grouping_entity(*region, meta, sset);
    bool internalSidesetFlag = true;

    if(nullptr != stkSideSetPart) {
        if(bulk.does_sideset_exist(*stkSideSetPart)) {
            stkSideSet = & bulk.get_sideset(*stkSideSetPart);
            stkSideSet->set_from_input(true);
        } else {
            stkSideSet = & bulk.create_sideset(*stkSideSetPart, true);
        }

        internalSidesetFlag = stkSideSet->get_accept_all_internal_non_coincident_entries();
        stkSideSet->set_accept_all_internal_non_coincident_entries(false);
    }

    size_t block_count = sset->block_count();
    for (size_t i=0; i < block_count; i++) {
        Ioss::SideBlock *block = sset->get_block(i);
        if (stk::io::include_entity(block)) {
            std::vector<INT> elem_side ;

            stk::mesh::Part *sbPart = get_part_for_grouping_entity(*region, meta, block);
            if (sbPart == nullptr) {
               sbPart = get_part_for_grouping_entity(*region, meta, sset);
            }

            stk::mesh::SideSet *sbSideSet = nullptr;
            if(nullptr != sbPart) {
                if(sbPart->id() != stkSideSetPart->id())
                  stk::RuntimeWarning() << "process_surface_entity: sideblock " << sbPart->name() << " with id " << sbPart->id()
                                        << " does not have the same id as parent sideset "
                                        << stkSideSetPart->name() << " with id " << stkSideSetPart->id();

                const stk::mesh::Part& sbParentPart = stk::mesh::get_sideset_parent(*sbPart);

                if(sbParentPart.mesh_meta_data_ordinal() != stkSideSetPart->mesh_meta_data_ordinal()) {
                  sbSideSet = & bulk.create_sideset(*sbPart);
                }
            }

            stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;

            // NOTE: Using the exodus sideset definition which is the
            // pair "element_id local_side_ordinal" will not correctly
            // identify embedded faces in a mesh.  For example, if
            // side 1 of element 2 is the same as side 3 of element 4,
            // this will result in two different faces being created
            // below instead of a single face since both faces share
            // the same nodal connectivity and should be the same
            // face.

            block->get_field_data("element_side", elem_side);
            stk::mesh::PartVector add_parts;

            if(nullptr != sbPart) {
                add_parts.push_back(sbPart);
            }

            // Get topology of the sides being defined to see if they
            // are 'faces' or 'edges'.  This is needed since for shell-type
            // elements, (and actually all elements) a sideset can specify either a face or an edge...
            // For a quad shell, sides 1,2 are faces and 3,4,5,6 are edges.

            // NOTE: This assumes that the sides within a side block are homogenous.  If the side_set
            //       is not split into homogenous side_blocks, then the topology will not necessarily
            //       be the same and this could fail (a sideset of mixed edges and faces)
            int par_dimen = block->topology()->parametric_dimension();
            // NOTE: Needed for shell side topology offset translation to IOSS
            std::int64_t sideSetOffset = Ioss::Utils::get_side_offset(block);

            size_t side_count = elem_side.size() / 2;
            for(size_t is=0; is<side_count; ++is) {
                stk::mesh::Entity const elem = bulk.get_entity(elem_rank, elem_side[is*2]);

                // If NULL, then the element was probably assigned to an
                // element block that appears in the database, but was
                // subsetted out of the analysis mesh. Only process if
                // non-null.
                if (bulk.is_valid(elem)) {
                    // Ioss uses 1-based side ordinal, stk::mesh uses 0-based.
                    int side_ordinal = elem_side[is*2+1] + sideSetOffset - 1;
                    stk::mesh::EntityId side_id_for_classic_behavior = stk::mesh::impl::side_id_formula(elem_side[is*2], side_ordinal);

                    if(par_dimen == 0)
                    {
                        stk::topology elemTopo = bulk.bucket(elem).topology();
                        stk::topology faceTopo = elemTopo.sub_topology(elemTopo.side_rank(), side_ordinal);

                        Ioss::ElementTopology *ioss_topo = Ioss::ElementTopology::factory(faceTopo.name(), false);
                        par_dimen = ioss_topo->parametric_dimension();
                    }

                    STK_ThrowRequireMsg((par_dimen == 1) || (par_dimen == 2), "Invalid value for par_dimen:" << par_dimen);

                    if(nullptr != sbSideSet) {
                        sbSideSet->add({elem, side_ordinal});
                    } else {
                         if(nullptr != stkSideSet) {
                           stkSideSet->add({elem, side_ordinal});
                         }
                    }

                    if (par_dimen == 1) {
                      stk::topology elemTopo = bulk.bucket(elem).topology();
                      // conversion from face ordinal to edge ordinal for shells
                      if (elemTopo.is_shell()) {
                        side_ordinal -= elemTopo.num_faces();
                      }
                      create_processed_edge(bulk, elem, side_ordinal, add_parts, behavior, side_id_for_classic_behavior);
                    }
                    else if (par_dimen == 2) {
                        create_processed_face(bulk, elem, side_ordinal, add_parts, behavior,
                                              side_id_for_classic_behavior, sidesToMove);
                    }
                }
            }
        }
    }
    if (stkSideSet != nullptr)
    {
        stkSideSet->set_accept_all_internal_non_coincident_entries(internalSidesetFlag);
        stk::util::sort_and_unique(*stkSideSet);
    }
}

void process_surface_entity(const Ioss::SideSet* sset, stk::mesh::BulkData & bulk, std::vector<ElemSidePartOrds> &sidesToMove, stk::io::StkMeshIoBroker::SideSetFaceCreationBehavior behavior)
{
  if (stk::io::db_api_int_size(sset) == 4) {
    process_surface_entity<int>(sset, bulk, sidesToMove, behavior);
  }
  else {
    process_surface_entity<int64_t>(sset, bulk, sidesToMove, behavior);
  }
}

template<typename INT>
void process_faces(stk::mesh::BulkData& bulk, std::vector<INT>& face_ids, std::vector<INT>& connectivity,
                  int nodes_per_face, size_t offset, stk::mesh::PartVector& faceParts)
{
    stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    stk::mesh::EntityVector faces;
    bulk.declare_entities(stk::topology::FACE_RANK, face_ids, faceParts, faces);
    stk::mesh::Part& nodePart = meta.get_topology_root_part(stk::topology::NODE);
    stk::mesh::PartVector nodeParts = {&nodePart};
    stk::mesh::Permutation perm = stk::mesh::Permutation::INVALID_PERMUTATION;
    stk::mesh::OrdinalVector scratch1, scratch2, scratch3;

    size_t face_count = face_ids.size();

    for(size_t i=0; i<face_count; ++i) {
        INT *conn = &connectivity[i*nodes_per_face];
        stk::mesh::Entity face = faces[i];

        bulk.set_local_id(face, offset + i);

        for(int j = 0; j < nodes_per_face; ++j)
        {
            stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, conn[j]);
            if (!bulk.is_valid(node)) {
                node = bulk.declare_node(conn[j], nodeParts);
            }
            bulk.declare_relation(face, node, j, perm, scratch1, scratch2, scratch3);
        }
        stk::mesh::impl::connect_face_to_elements(bulk, face);
        stk::mesh::impl::connect_face_to_edges(bulk, face);
    }
}

template<typename INT>
void process_edges(stk::mesh::BulkData& bulk, std::vector<INT>& edge_ids, std::vector<INT>& connectivity,
                  int nodes_per_edge, size_t offset, stk::mesh::PartVector& edgeParts)
{
    stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    stk::mesh::EntityVector edges;
    bulk.declare_entities(stk::topology::EDGE_RANK, edge_ids, edgeParts, edges);
    stk::mesh::Part& nodePart = meta.get_topology_root_part(stk::topology::NODE);
    stk::mesh::PartVector nodeParts = {&nodePart};
    stk::mesh::Permutation perm = stk::mesh::Permutation::INVALID_PERMUTATION;
    stk::mesh::OrdinalVector scratch1, scratch2, scratch3;

    size_t edge_count = edge_ids.size();

    for(size_t i=0; i<edge_count; ++i) {
        INT *conn = &connectivity[i*nodes_per_edge];
        stk::mesh::Entity edge = edges[i];

        bulk.set_local_id(edge, offset + i);

        for(int j = 0; j < nodes_per_edge; ++j)
        {
            stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, conn[j]);
            if (!bulk.is_valid(node)) {
                node = bulk.declare_node(conn[j], nodeParts);
            }
            bulk.declare_relation(edge, node, j, perm, scratch1, scratch2, scratch3);
        }
        bool connectedElems = stk::mesh::impl::connect_edge_to_elements(bulk, edge);
        if (!connectedElems) {
          bulk.destroy_entity(edge);
        }
    }
}

void create_shared_edges(stk::mesh::BulkData& bulk, stk::mesh::Part* part)
{
    stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    if(meta.side_rank() == stk::topology::EDGE_RANK) { return; }

    stk::mesh::EntityVector edges;
    stk::mesh::get_entities(bulk, stk::topology::EDGE_RANK, *part, edges);
    stk::mesh::EntityVector edgeNodes;
    std::vector<int> sharedProcs;
    stk::CommSparse commSparse(bulk.parallel());

    pack_and_communicate(commSparse, [&commSparse, &bulk, &sharedProcs, &edgeNodes, &edges]()
                        {
                            for(auto edge : edges) {
                                bool isShared = true;
                                const stk::mesh::Entity* nodes = bulk.begin_nodes(edge);
                                unsigned numNodes = bulk.num_nodes(edge);
                                edgeNodes.resize(numNodes);

                                for(unsigned i = 0; i < numNodes; i++) {
                                    isShared |= bulk.bucket(nodes[i]).shared();
                                    edgeNodes[i] = nodes[i];
                                }

                                if(isShared) {
                                    bulk.shared_procs_intersection(edgeNodes, sharedProcs);

                                    for(int proc : sharedProcs) {
                                        stk::CommBuffer& buf = commSparse.send_buffer(proc);
                                        buf.pack<stk::mesh::EntityId>(bulk.identifier(edge));
                                        buf.pack<unsigned>(numNodes);

                                        for(auto edgeNode : edgeNodes) {
                                            buf.pack<stk::mesh::EntityKey>(bulk.entity_key(edgeNode));
                                        }
                                    }
                                }
                            }  
                        });
                    
    std::vector<stk::mesh::EntityId> edgeIds;
    std::vector<stk::mesh::EntityId> connectivity;
    unsigned nodesPerEdge = 0;

    unpack_communications(commSparse, [&commSparse, &nodesPerEdge, &edgeIds, &connectivity](int procId)
                         {
                             stk::mesh::EntityKey nodeKey;
                             stk::mesh::EntityId edgeId;
                             unsigned numNodes;
    
                             commSparse.recv_buffer(procId).unpack<stk::mesh::EntityId>(edgeId);
                             edgeIds.push_back(edgeId);
                             commSparse.recv_buffer(procId).unpack<unsigned>(numNodes);
    
                             if(nodesPerEdge == 0) {
                                 nodesPerEdge = numNodes;
                             }
                             STK_ThrowAssert(numNodes == nodesPerEdge);
    
                             for(unsigned i = 0; i < numNodes; i++) {
                                 commSparse.recv_buffer(procId).unpack<stk::mesh::EntityKey>(nodeKey);
                                 connectivity.push_back(nodeKey.id());
                             }
                         });

    stk::mesh::PartVector partVector{part};
    process_edges<stk::mesh::EntityId>(bulk, edgeIds, connectivity, nodesPerEdge, 0, partVector);
}

template <typename INT>
void process_face_entity(const Ioss::FaceBlock* fb, stk::mesh::BulkData & bulk)
{
    assert(fb->type() == Ioss::FACEBLOCK);

    const stk::mesh::MetaData &meta = bulk.mesh_meta_data();

    Ioss::Region *region = fb->get_database()->get_region();

    stk::mesh::Part *part = get_part_for_grouping_entity(*region, meta, fb);
    STK_ThrowRequireMsg(part != nullptr, " INTERNAL_ERROR: Part for edge block: " << fb->name() << " does not exist");
    stk::mesh::PartVector faceParts = {part};

    stk::topology topo = part->topology();
    STK_ThrowRequireMsg( topo != stk::topology::INVALID_TOPOLOGY, " INTERNAL_ERROR: Part " << part->name() << " has invalid topology");

    std::vector<INT> face_ids;
    std::vector<INT> connectivity ;

    fb->get_field_data("ids", face_ids);
    fb->get_field_data("connectivity", connectivity);

    int nodes_per_face = topo.num_nodes();
    size_t offset = fb->get_offset();
    process_faces<INT>(bulk, face_ids, connectivity, nodes_per_face, offset, faceParts);
}

template <typename INT>
void process_edge_entity(const Ioss::EdgeBlock* eb, stk::mesh::BulkData & bulk)
{
    assert(eb->type() == Ioss::EDGEBLOCK);

    const stk::mesh::MetaData &meta = bulk.mesh_meta_data();

    Ioss::Region *region = eb->get_database()->get_region();

    stk::mesh::Part *part = get_part_for_grouping_entity(*region, meta, eb);
    STK_ThrowRequireMsg(part != nullptr, " INTERNAL_ERROR: Part for edge block: " << eb->name() << " does not exist");
    stk::mesh::PartVector edgeParts = {part};

    stk::topology topo = part->topology();
    STK_ThrowRequireMsg( topo != stk::topology::INVALID_TOPOLOGY, " INTERNAL_ERROR: Part " << part->name() << " has invalid topology");

    std::vector<INT> edge_ids;
    std::vector<INT> connectivity ;

    eb->get_field_data("ids", edge_ids);
    eb->get_field_data("connectivity", connectivity);

    int nodes_per_edge = topo.num_nodes();
    size_t offset = eb->get_offset();
    process_edges<INT>(bulk, edge_ids, connectivity, nodes_per_edge, offset, edgeParts);

    create_shared_edges(bulk, part);
}


void process_face_entity(const Ioss::FaceBlock* fb, stk::mesh::BulkData & bulk)
{
  if (stk::io::db_api_int_size(fb) == 4) {
    process_face_entity<int>(fb, bulk);
  }
  else {
    process_face_entity<int64_t>(fb, bulk);
  }
}

void process_edge_entity(const Ioss::EdgeBlock* eb, stk::mesh::BulkData & bulk)
{
  if (stk::io::db_api_int_size(eb) == 4) {
    process_edge_entity<int>(eb, bulk);
  }
  else {
    process_edge_entity<int64_t>(eb, bulk);
  }
}

void send_element_side_to_element_owner(stk::CommSparse &comm,
                                        const stk::mesh::EntityIdProcMap &elemIdMovedToProc,
                                        const std::vector<ElemSidePartOrds> &sidesToMove)
{
    for(const ElemSidePartOrds& sideToMove : sidesToMove)
    {
        const stk::mesh::EntityIdProcMap::const_iterator iter = elemIdMovedToProc.find(sideToMove.elem);
        STK_ThrowRequireWithSierraHelpMsg(iter!=elemIdMovedToProc.end());
        int destProc = iter->second;
        comm.send_buffer(destProc).pack(sideToMove.elem);
        comm.send_buffer(destProc).pack(sideToMove.sideOrdinal);
        pack_vector_to_proc(comm, sideToMove.partOrdinals, destProc);
    }
}

void unpack_and_declare_element_side(stk::CommSparse & comm, stk::mesh::BulkData & bulk, int procId)
{
    stk::mesh::EntityId elemId;
    unsigned sideOrdinal;
    comm.recv_buffer(procId).unpack(elemId);
    comm.recv_buffer(procId).unpack(sideOrdinal);
    stk::mesh::OrdinalVector partOrdinals;
    unpack_vector_from_proc(comm, partOrdinals, procId);

    stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, elemId);
    STK_ThrowRequireWithSierraHelpMsg(bulk.is_valid(elem));
    stk::mesh::PartVector add_parts;
    stk::mesh::impl::convert_part_ordinals_to_parts(bulk.mesh_meta_data(), partOrdinals, add_parts);
    bulk.declare_element_side(elem, sideOrdinal, add_parts);
}

void move_sideset_to_follow_element(stk::mesh::BulkData &bulk, const stk::mesh::EntityIdProcMap &elemIdMovedToProc, const std::vector<ElemSidePartOrds> &sidesToMove)
{
    stk::CommSparse comm(bulk.parallel());
    stk::pack_and_communicate(comm, [&comm, &elemIdMovedToProc, &sidesToMove]() { send_element_side_to_element_owner(comm, elemIdMovedToProc, sidesToMove); });
    stk::unpack_communications(comm, [&comm, &bulk](int procId) { unpack_and_declare_element_side(comm, bulk, procId); });
}

void process_sidesets(Ioss::Region &region, stk::mesh::BulkData &bulk, const stk::mesh::EntityIdProcMap &elemIdMovedToProc, stk::io::StkMeshIoBroker::SideSetFaceCreationBehavior behavior)
{
    const Ioss::SideSetContainer& side_sets = region.get_sidesets();

    std::vector<ElemSidePartOrds> sidesToMove;
    for(Ioss::SideSetContainer::const_iterator it = side_sets.begin(); it != side_sets.end(); ++it)
        if(stk::io::include_entity(*it))
            process_surface_entity(*it, bulk, sidesToMove, behavior);

    move_sideset_to_follow_element(bulk, elemIdMovedToProc, sidesToMove);
}

void process_face_blocks(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
    const Ioss::FaceBlockContainer& face_blocks = region.get_face_blocks();

    for(const Ioss::FaceBlock* faceBlock : face_blocks) {
        if(stk::io::include_entity(faceBlock)) {
            process_face_entity(faceBlock, bulk);
        }
    }
}

void process_edge_blocks(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
    const Ioss::EdgeBlockContainer& edge_blocks = region.get_edge_blocks();

    for(const Ioss::EdgeBlock* edgeBlock : edge_blocks) {
        if(stk::io::include_entity(edgeBlock)) {
            process_edge_entity(edgeBlock, bulk);
        }
    }
}

void process_face_blocks(Ioss::Region &region, stk::mesh::MetaData &meta, TopologyErrorHandler handler)
{
  const Ioss::FaceBlockContainer& face_blocks = region.get_face_blocks();
  stk::io::default_part_processing(face_blocks, meta, handler);
}

void process_edge_blocks(Ioss::Region &region, stk::mesh::MetaData &meta, TopologyErrorHandler handler)
{
  const Ioss::EdgeBlockContainer& edge_blocks = region.get_edge_blocks();
  stk::io::default_part_processing(edge_blocks, meta, handler);
}

void process_assemblies(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::AssemblyContainer& assemblies = region.get_assemblies();
  stk::io::default_part_processing(assemblies, meta);
}

void build_assembly_hierarchies(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::AssemblyContainer& assemblies = region.get_assemblies();
  for(const Ioss::Assembly* assembly : assemblies) {
    if(!include_entity(assembly)) {continue;}

    const std::string& assemblyName = assembly->name();
    stk::mesh::Part* assemblyPart = meta.get_part(assemblyName);
    const Ioss::EntityContainer& members = assembly->get_members();
    for(const Ioss::GroupingEntity* member : members) {
      if(!include_entity(member)) {continue;}

      stk::mesh::Part* memberPart = meta.get_part(member->name());
      meta.declare_part_subset(*assemblyPart, *memberPart);
    }
  }
}

void populate_hidden_nodesets(Ioss::Region &io, const stk::mesh::MetaData & meta, NodesetMap &nodesetMap)
{
    static const std::string s_nodeset_suffix("_n");

    for(auto mesh_part : meta.get_mesh_parts())
    {
        if(stk::io::is_part_io_part(*mesh_part))
        {
            const std::string& part_name(mesh_part->name());

            // Determine the type of io entity and load into the pairing vector
            auto entity = io.get_entity(part_name);
            if(entity == nullptr) {
                continue;
            }

            if(entity->type() == Ioss::SIDEBLOCK || entity->type() == Ioss::SIDESET)
            {
                std::string ns_name = part_name;
                ns_name += s_nodeset_suffix;
                auto io_node_set = io.get_nodeset(ns_name);
                if(io_node_set != nullptr) {
                    nodesetMap[io_node_set] = mesh_part;
                }
            }
        }
    }
}

stk::mesh::Part* get_part_from_alias(const Ioss::Region &region, const stk::mesh::MetaData &meta, const std::string &name, Ioss::EntityType type)
{
    stk::mesh::Part* part = nullptr;

    std::vector<std::string> aliases;
    region.get_aliases(name, type, aliases);

    for(std::string &alias : aliases)
    {
        if(sierra::case_strcmp(alias, name) != 0)
        {
            part = meta.get_part(alias);

            if(nullptr != part) {
                break;
            }
        }
    }

    return part;
}

stk::mesh::Part* get_part_for_grouping_entity(const Ioss::Region &region, const stk::mesh::MetaData &meta, const Ioss::GroupingEntity *entity)
{
    const std::string &name = entity->name();
    auto type = entity->type();
    stk::mesh::Part* part = meta.get_part(name);

    if(nullptr == part) {
        part = get_part_from_alias(region, meta, name, type);
    }
    return part;
}


}} // namespace stk io
