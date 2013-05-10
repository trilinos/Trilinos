/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_FEMHelpers_hpp
#define stk_mesh_FEMHelpers_hpp

#include <stk_mesh/base/Types.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/CellTopology.hpp>
// This is needed for ElementNode class
#include <stk_mesh/base/TopologyDimensions.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_bulk_data_element
 *  \{
 */

//----------------------------------------------------------------------
/** \brief  Declare an element member of a Part with a CellTopology
 *          and nodes conformal to that topology.
 */
Entity declare_element( BulkData & mesh ,
                          Part & part ,
                          const EntityId elem_id ,
                          const EntityId node_id[] );


/** \brief  Create (or find) an element side.
 *
 *  The element must be a member of a Part with a CellTopology.
 */
Entity declare_element_side( BulkData & mesh ,
                               const stk::mesh::EntityId global_side_id ,
                               Entity elem ,
                               const unsigned local_side_id ,
                               Part * part = NULL );

/** \brief  Create (or find) an element edge.
 *
 *  The element must be a member of a Part with a CellTopology.
 */
Entity declare_element_edge( BulkData & mesh ,
                               const stk::mesh::EntityId global_side_id ,
                               Entity elem ,
                               const unsigned local_side_id ,
                               Part * part = NULL );

/** \brief  Create (or find) an element side.
 *
 *  The element must be a member of a Part with a CellTopology.
 */
Entity declare_element_side( BulkData & mesh ,
                               Entity elem ,
                               Entity side ,
                               const unsigned local_side_id ,
                               Part * part = NULL );



/** \brief  Create (or find) an element edge.
 *
 *  The element must be a member of a Part with a CellTopology.
 */
Entity declare_element_edge( BulkData & mesh ,
                               Entity elem ,
                               Entity edge ,
                               const unsigned local_edge_id ,
                               Part * part = NULL );



/** \brief  Declare a part with a given cell topology. This is just a convenient
            function that wraps MetaData's declare_part.
 */
template< class Top >
Part &declare_part(MetaData& meta_data, const std::string &name) {
  return meta_data.declare_part(name, shards::getCellTopologyData<Top>());
}

/**
 * Given an entity, subcell_rank, and subcell_id, return the nodes
 * that make up the subcell in a correct order for the given polarity.
 *
 * \param entity
 * \param subcell_rank
 * \param subcell_indentifier
 * \param subcell_nodes EntityVector output of the subcell nodes
 * \param use_reverse_polarity
 * \return CellTopologyData * of the requested subcell
 */
const CellTopologyData * get_subcell_nodes(const BulkData& mesh,
    const Entity entity ,
    EntityRank         subcell_rank ,
    unsigned           subcell_identifier ,
    EntityVector     & subcell_nodes
    );

/** \brief  Given an entity and collection of nodes, return the
 *          local id of the subcell that contains those nodes in the
 *          correct orientation.
 */
int get_entity_subcell_id( const BulkData& mesh, const Entity entity ,
                           const EntityRank          subcell_rank,
                           const CellTopologyData  * side_topology,
                           const EntityVector      & side_nodes );


template< class Traits >
void get_parts_with_topology(stk::mesh::BulkData& mesh,
                             stk::mesh::PartVector& parts,
                             bool skip_topology_root_parts=false)
{
  parts.clear();

  stk::mesh::MetaData & fem_meta = stk::mesh::MetaData::get(mesh);

  const stk::mesh::PartVector& all_parts = fem_meta.get_parts();

  stk::mesh::PartVector::const_iterator
    iter = all_parts.begin(),
    iter_end = all_parts.end();

  const CellTopologyData* topology = shards::getCellTopologyData<Traits>();

  for(; iter!=iter_end; ++iter) {
    stk::mesh::Part* part =  *iter;
    if (fem_meta.get_cell_topology(*part).getCellTopologyData() == topology) {
      if (skip_topology_root_parts && stk::mesh::is_cell_topology_root_part(*part)) {
        continue;
      }
      parts.push_back(part);
    }
  }
}

/** \} */

} //namespace mesh
} //namespace stk
#endif
