#include "ProcessSetsOrBlocks.hpp"
#include <iosfwd>
#include <map>
#include "Ioss_Field.h"
#include "Ioss_SideSet.h"
#include "IossBridge.hpp"

#include "StkMeshIoBroker.hpp"
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field, etc
#include "stk_topology/topology.hpp"
#include "stk_mesh/base/CoordinateSystems.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include "stk_util/environment/ReportHandler.hpp"
#include "stk_mesh/base/TopologyDimensions.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"

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
  stk::mesh::put_field(coord_field, meta.universal_part(),
                       meta.spatial_dimension());
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

  stk::mesh::Field<double> & distribution_factors_field =
    meta.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "distribution_factors");
  stk::io::set_field_role(distribution_factors_field, Ioss::Field::MESH);

  /** \todo REFACTOR How to associate distribution_factors field
   * with the nodeset part if a node is a member of multiple
   * nodesets
   */

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
      it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      stk::mesh::Part* const part = meta.get_part(entity->name());

      STKIORequire(part != NULL);
      STKIORequire(entity->field_exists("distribution_factors"));

      stk::io::set_field_role(distribution_factors_field, Ioss::Field::MESH);
      stk::mesh::put_field(distribution_factors_field, *part);
    }
  }

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
      it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      stk::mesh::Part* const part = meta.get_part(entity->name());

      STKIORequire(part != NULL);
      STKIORequire(entity->field_exists("distribution_factors"));

      std::string nodesetName = part->name();
      std::string nodesetDistFieldName = "distribution_factors_" + nodesetName;

      stk::mesh::Field<double> & distribution_factors_field_per_nodeset =
        meta.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, nodesetDistFieldName);

      stk::io::set_field_role(distribution_factors_field_per_nodeset, Ioss::Field::MESH);
      stk::mesh::put_field(distribution_factors_field_per_nodeset, *part);
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
  STKIORequire(ss_part != NULL);

  stk::mesh::Field<double, stk::mesh::ElementNode> *distribution_factors_field = NULL;
  bool surface_df_defined = false; // Has the surface df field been defined yet?

  size_t block_count = sset->block_count();
  for (size_t i=0; i < block_count; i++) {
    Ioss::SideBlock *sb = sset->get_block(i);
    if (stk::io::include_entity(sb)) {
      stk::mesh::Part * const sb_part = meta.get_part(sb->name());
      STKIORequire(sb_part != NULL);
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
        stk::mesh::put_field(*distribution_factors_field,
                             *sb_part, side_node_count);
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

// ========================================================================


// ========================================================================
stk::mesh::EntityId get_side_entity_id(int64_t elem_id, int side_ordinal)
{
  // NOTE: This function uses a 1-based side ordinal
  int64_t ten = 10;
  stk::mesh::EntityId side_id = elem_id * ten + side_ordinal;
  return side_id;
}

template <typename INT>
void process_surface_entity(const Ioss::SideSet* sset, stk::mesh::BulkData & bulk, std::vector<ElemSidePartOrds> &sidesToMove, stk::io::StkMeshIoBroker::SideSetFaceCreationBehavior behavior)
{
    assert(sset->type() == Ioss::SIDESET);

    const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(bulk);

    size_t block_count = sset->block_count();
    for (size_t i=0; i < block_count; i++) {
        Ioss::SideBlock *block = sset->get_block(i);
        if (stk::io::include_entity(block)) {
            std::vector<INT> elem_side ;

            stk::mesh::Part * const sb_part = meta.get_part(block->name());
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
            stk::mesh::PartVector add_parts( 1 , sb_part );

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
                    stk::mesh::EntityId side_id = get_side_entity_id(elem_side[is*2], elem_side[is*2+1]);

                    if (par_dimen == 1) {
                        if(bulk.mesh_meta_data().spatial_dimension()==2 && behavior == stk::io::StkMeshIoBroker::STK_IO_SIDE_CREATION_USING_GRAPH_TEST)
                        {
                            stk::mesh::declare_element_side(bulk, elem, side_ordinal, add_parts);
                        }
                        else
                        {
                            stk::mesh::Entity side = stk::mesh::declare_element_edge(bulk, side_id, elem, side_ordinal);
                            bulk.change_entity_parts( side, add_parts );
                        }
                    }
                    else if (par_dimen == 2) {
                        if (behavior == stk::io::StkMeshIoBroker::STK_IO_SIDESET_FACE_CREATION_CLASSIC) {
                            stk::mesh::Entity side = stk::mesh::declare_element_side(bulk, side_id, elem, side_ordinal);
                            bulk.change_entity_parts( side, add_parts );
                        }
                        else if (behavior == stk::io::StkMeshIoBroker::STK_IO_SIDESET_FACE_CREATION_CURRENT) {
                            stk::mesh::Entity new_face = stk::mesh::impl::get_or_create_face_at_element_side(bulk,elem,side_ordinal,side_id,stk::mesh::PartVector(1,sb_part));
                            stk::mesh::impl::connect_face_to_other_elements(bulk,new_face,elem,side_ordinal);
                        }
                        else if (behavior == stk::io::StkMeshIoBroker::STK_IO_SIDE_CREATION_USING_GRAPH_TEST) {
                            if(bulk.bucket(elem).owned()) {
                                stk::mesh::declare_element_side(bulk, elem, side_ordinal, add_parts);
                            }
                            else
                            {
                                stk::mesh::OrdinalVector ords(add_parts.size());
                                for(size_t n=0; n<add_parts.size(); ++n) {
                                    ords[n] = add_parts[n]->mesh_meta_data_ordinal();
                                }
                                sidesToMove.push_back({bulk.identifier(elem), side_ordinal, ords});
                            }
                        }
                    }
                }
            }
        }
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
        comm.send_buffer(destProc).pack<stk::mesh::EntityId>(sideToMove.elem);
        comm.send_buffer(destProc).pack<unsigned>(sideToMove.sideOrdinal);
        pack_vector_to_proc(comm, sideToMove.partOrdinals, destProc);
    }
}

void unpack_and_declare_element_side(stk::CommSparse & comm, stk::mesh::BulkData & bulk, int procId)
{
    stk::mesh::EntityId elemId;
    unsigned sideOrdinal;
    comm.recv_buffer(procId).unpack<stk::mesh::EntityId>(elemId);
    comm.recv_buffer(procId).unpack<unsigned>(sideOrdinal);
    stk::mesh::OrdinalVector partOrdinals;
    unpack_vector_from_proc(comm, partOrdinals, procId);

    stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, elemId);
    ThrowRequireWithSierraHelpMsg(bulk.is_valid(elem));
    stk::mesh::PartVector add_parts;
    stk::mesh::impl::convert_part_ordinals_to_parts(bulk.mesh_meta_data(), partOrdinals, add_parts);
    stk::mesh::declare_element_side(bulk, elem, sideOrdinal, add_parts);
}

void move_sidset_to_follow_element(stk::mesh::BulkData &bulk, const stk::mesh::EntityIdProcMap &elemIdMovedToProc, const std::vector<ElemSidePartOrds> &sidesToMove)
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

    move_sidset_to_follow_element(bulk, elemIdMovedToProc, sidesToMove);
}

}} // namespace stk io
