
// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "ProcessSetsOrBlocks.hpp"
#include <cstdint>                                 // for int64_t
#include <map>                                     // for allocator, etc
#include <stk_mesh/base/BulkData.hpp>              // for BulkData
#include <stk_mesh/base/MetaData.hpp>              // for MetaData, etc
#include <utility>                                 // for pair

#include "stk_util/util/SortAndUnique.hpp"
#include "IossBridge.hpp"                          // for include_entity, etc
#include "Ioss_ElementTopology.h"                  // for ElementTopology
#include "Ioss_Field.h"                            // for Field, etc
#include "Ioss_NodeBlock.h"                        // for NodeBlock
#include "Ioss_SideBlock.h"                        // for SideBlock
#include "Ioss_SideSet.h"                          // for SideSet, etc
#include "StkIoUtils.hpp"
#include "StkMeshIoBroker.hpp"                     // for StkMeshIoBroker, etc
#include "stk_mesh/base/Bucket.hpp"                // for Bucket
#include "stk_mesh/base/CoordinateSystems.hpp"     // for Cartesian
#include "stk_mesh/base/Field.hpp"                 // for Field
#include "stk_mesh/base/TopologyDimensions.hpp"    // for ElementNode
#include "stk_mesh/base/Types.hpp"                 // for OrdinalVector, etc
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_topology/topology.hpp"               // for topology, etc
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/diag/StringUtil.hpp"
#include "stk_util/parallel/CommSparse.hpp"        // for CommSparse, etc
#include "stk_util/parallel/ParallelComm.hpp"      // for CommBuffer
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

  stk::mesh::Field<double, stk::mesh::Cartesian>& coord_field =
    meta.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian> >(stk::topology::NODE_RANK, stk::io::CoordinateFieldName);
  stk::io::set_field_role(coord_field, Ioss::Field::MESH);

  meta.set_coordinate_field(&coord_field);

  Ioss::NodeBlock *nb = node_blocks[0];
  stk::mesh::put_field_on_mesh(coord_field, meta.universal_part(),
                               meta.spatial_dimension(),
                               (stk::mesh::FieldTraits<stk::mesh::Field<double, stk::mesh::Cartesian>>::data_type*) nullptr);
  stk::io::define_io_fields(nb, Ioss::Field::ATTRIBUTE, meta.universal_part(), stk::topology::NODE_RANK);
}



void process_elementblocks(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
  stk::io::default_part_processing(elem_blocks, meta);
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
        meta.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, nodesetDistFieldName);

      stk::io::set_field_role(distribution_factors_field_per_nodeset, Ioss::Field::MESH);
      stk::mesh::put_field_on_mesh(distribution_factors_field_per_nodeset, *part,
                                   (stk::mesh::FieldTraits<stk::mesh::Field<double>>::data_type*) nullptr);
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

  stk::mesh::Field<double, stk::mesh::ElementNode> *distribution_factors_field = nullptr;
  bool surface_df_defined = false; // Has the surface df field been defined yet?

  size_t block_count = sset->block_count();
  for (size_t i=0; i < block_count; i++) {
    Ioss::SideBlock *sb = sset->get_block(i);
    if (stk::io::include_entity(sb)) {
      stk::mesh::Part * const sb_part = meta.get_part(sb->name());
      STKIORequire(sb_part != nullptr);
      meta.declare_part_subset(*ss_part, *sb_part);

      if (sb->field_exists("distribution_factors")) {
        if (!surface_df_defined) {
          stk::topology::rank_t side_rank = static_cast<stk::topology::rank_t>(stk::io::part_primary_entity_rank(*sb_part));
          std::string field_name = sset->name() + "_df";
          distribution_factors_field =
            &meta.declare_field<stk::mesh::Field<double, stk::mesh::ElementNode> >(side_rank, field_name);
          stk::io::set_field_role(*distribution_factors_field, Ioss::Field::MESH);
          stk::io::set_distribution_factor_field(*ss_part, *distribution_factors_field);
          surface_df_defined = true;
        }
        stk::io::set_distribution_factor_field(*sb_part, *distribution_factors_field);
        int side_node_count = sb->topology()->number_nodes();
        stk::mesh::put_field_on_mesh(*distribution_factors_field,
                                     *sb_part, side_node_count,
                                     (stk::mesh::FieldTraits<stk::mesh::Field<double, stk::mesh::ElementNode>>::data_type*) nullptr);
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

template <typename INT>
void process_surface_entity(const Ioss::SideSet* sset, stk::mesh::BulkData & bulk, std::vector<ElemSidePartOrds> &sidesToMove, stk::io::StkMeshIoBroker::SideSetFaceCreationBehavior behavior)
{
    assert(sset->type() == Ioss::SIDESET);

    const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(bulk);

    Ioss::Region *region = sset->get_database()->get_region();
//    const std::string universalAlias = region->get_alias("universal_sideset");
//    if (sset->name() == universalAlias)
//        return;

    stk::mesh::SideSet *stkSideSet = nullptr;
    stk::mesh::Part *stkSideSetPart = get_part_for_grouping_entity(*region, meta, sset);

    if(nullptr != stkSideSetPart)
    {
        if(bulk.does_sideset_exist(*stkSideSetPart))
        {
            stkSideSet = & bulk.get_sideset(*stkSideSetPart);
        }
        else
        {
            stkSideSet = & bulk.create_sideset(*stkSideSetPart, true);
        }
    }

    size_t block_count = sset->block_count();
    for (size_t i=0; i < block_count; i++) {
        Ioss::SideBlock *block = sset->get_block(i);
        if (stk::io::include_entity(block)) {
            std::vector<INT> elem_side ;

            stk::mesh::Part *sb_part = get_part_for_grouping_entity(*region, meta, block);
            if (sb_part == nullptr)
            {
               sb_part = get_part_for_grouping_entity(*region, meta, sset);
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

            if(nullptr != sb_part) {
                add_parts.push_back(sb_part);
            }

            // Get topology of the sides being defined to see if they
            // are 'faces' or 'edges'.  This is needed since for shell-type
            // elements, (and actually all elements) a sideset can specify either a face or an edge...
            // For a quad shell, sides 1,2 are faces and 3,4,5,6 are edges.

            // NOTE: This assumes that the sides within a side block are homogenous.  If the side_set
            //       is not split into homogenous side_blocks, then the topology will not necessarily
            //       be the same and this could fail (a sideset of mixed edges and faces)
            int par_dimen = block->topology()->parametric_dimension();

            size_t side_count = elem_side.size() / 2;
            for(size_t is=0; is<side_count; ++is) {
                stk::mesh::Entity const elem = bulk.get_entity(elem_rank, elem_side[is*2]);

                // If NULL, then the element was probably assigned to an
                // element block that appears in the database, but was
                // subsetted out of the analysis mesh. Only process if
                // non-null.
                if (bulk.is_valid(elem)) {
                    // Ioss uses 1-based side ordinal, stk::mesh uses 0-based.
                    int side_ordinal = elem_side[is*2+1] - 1;
                    stk::mesh::EntityId side_id_for_classic_behavior = stk::mesh::impl::side_id_formula(elem_side[is*2], side_ordinal);

                    if(par_dimen == 0)
                    {
                        stk::topology elemTopo = bulk.bucket(elem).topology();
                        stk::topology faceTopo = elemTopo.sub_topology(elemTopo.side_rank(), side_ordinal);

                        Ioss::ElementTopology *ioss_topo = Ioss::ElementTopology::factory(faceTopo.name(), false);
                        par_dimen = ioss_topo->parametric_dimension();
                    }

                    ThrowRequireMsg((par_dimen == 1) || (par_dimen == 2), "Invalid value for par_dimen:" << par_dimen);

                    if(nullptr != stkSideSet)
                    {
                        stkSideSet->add({elem, side_ordinal});
                    }

                    if (par_dimen == 1) {
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


void send_element_side_to_element_owner(stk::CommSparse &comm,
                                        const stk::mesh::EntityIdProcMap &elemIdMovedToProc,
                                        const std::vector<ElemSidePartOrds> &sidesToMove)
{
    for(const ElemSidePartOrds& sideToMove : sidesToMove)
    {
        const stk::mesh::EntityIdProcMap::const_iterator iter = elemIdMovedToProc.find(sideToMove.elem);
        ThrowRequireWithSierraHelpMsg(iter!=elemIdMovedToProc.end());
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
    ThrowRequireWithSierraHelpMsg(bulk.is_valid(elem));
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

stk::mesh::Part* get_part_from_alias(const Ioss::Region &region, const stk::mesh::MetaData &meta, const std::string &name)
{
    stk::mesh::Part* part = nullptr;

    std::vector<std::string> aliases;
    region.get_aliases(name, aliases);

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
    stk::mesh::Part* part = meta.get_part(name);

    if(nullptr == part) {
        part = get_part_from_alias(region, meta, name);
    }
    return part;
}


}} // namespace stk io
