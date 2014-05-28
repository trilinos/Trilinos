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

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/CellTopology.hpp>
// This is needed for ElementNode class
#include <stk_mesh/fem/TopologyDimensions.hpp>

namespace stk {
namespace mesh {

class Bucket;
class Entity;

namespace fem {

/** \addtogroup stk_mesh_bulk_data_element
 *  \{
 */

//----------------------------------------------------------------------
/** \brief  Declare an element member of a Part with a CellTopology
 *          and nodes conformal to that topology.
 */
Entity & declare_element( BulkData & mesh ,
                          Part & part ,
                          const EntityId elem_id ,
                          const EntityId node_id[] );


/** \brief  Create (or find) an element side.
 *
 *  The element must be a member of a Part with a CellTopology.
 */
Entity & declare_element_side( BulkData & mesh ,
                               const stk::mesh::EntityId global_side_id ,
                               Entity & elem ,
                               const unsigned local_side_id ,
                               Part * part = NULL );

/** \brief  Create (or find) an element edge.
 *
 *  The element must be a member of a Part with a CellTopology.
 */
Entity & declare_element_edge( BulkData & mesh ,
                               const stk::mesh::EntityId global_side_id ,
                               Entity & elem ,
                               const unsigned local_side_id ,
                               Part * part = NULL );

/** \brief  Determine the polarity of the local side,
 *          more efficient if the local_side_id is known.
 */
bool element_side_polarity( const Entity & elem ,
                            const Entity & side , int local_side_id = -1 );

/** \brief  Create (or find) an element side.
 *
 *  The element must be a member of a Part with a CellTopology.
 */
Entity & declare_element_side( Entity & elem ,
                               Entity & side ,
                               const unsigned local_side_id ,
                               Part * part = NULL );



/** \brief  Create (or find) an element edge.
 *
 *  The element must be a member of a Part with a CellTopology.
 */
Entity & declare_element_edge( Entity & elem ,
                               Entity & edge ,
                               const unsigned local_edge_id ,
                               Part * part = NULL );



/** \brief  Declare a part with a given cell topology. This is just a convenient
            function that wraps FEMMetaData's declare_part.
 */
template< class Top >
Part &declare_part(FEMMetaData& meta_data, const std::string &name) {
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
const CellTopologyData * get_subcell_nodes(
    const Entity     & entity ,
    EntityRank         subcell_rank ,
    unsigned           subcell_identifier ,
    EntityVector     & subcell_nodes
    );

/** \brief  Given an entity and collection of nodes, return the
 *          local id of the subcell that contains those nodes in the
 *          correct orientation.
 */
int get_entity_subcell_id( const Entity            & entity ,
                           const EntityRank          subcell_rank,
                           const CellTopologyData  * side_topology,
                           const EntityVector      & side_nodes );

/** \brief Global counts for a mesh's entities. */
bool comm_mesh_counts( BulkData & ,
                       std::vector<size_t> & counts ,
                       bool = false );

typedef Field<double*,stk::mesh::ElementNode> ElementNodePointerField ;

/** \brief  Declare an element-to-node-data pointer field.
 */
template< class NodeField >
inline
ElementNodePointerField &
declare_element_node_pointer_field(
  FEMMetaData & fmd , const std::string & s ,
  NodeField & node_field )
{
  const unsigned num_states = node_field.number_of_states();

  ElementNodePointerField & f =
    fmd.template declare_field< ElementNodePointerField >( s, num_states );

  for ( unsigned i = 0 ; i < num_states ; ++i ) {
    FieldState state = (FieldState) i;
    fmd.declare_field_relation(
      f.field_of_state( state ) ,
      fem::get_element_node_stencil(fmd.spatial_dimension()) ,
      node_field.field_of_state( state ) );
  }

  return f ;
}

template< class Traits >
void get_parts_with_topology(stk::mesh::BulkData& mesh,
                             stk::mesh::PartVector& parts,
                             bool skip_topology_root_parts=false)
{
  parts.clear();

  stk::mesh::fem::FEMMetaData & fem_meta = stk::mesh::fem::FEMMetaData::get(mesh);

  const stk::mesh::PartVector& all_parts = fem_meta.get_parts();

  stk::mesh::PartVector::const_iterator
    iter = all_parts.begin(),
    iter_end = all_parts.end();

  const CellTopologyData* topology = shards::getCellTopologyData<Traits>();

  for(; iter!=iter_end; ++iter) {
    stk::mesh::Part* part =  *iter;
    if (fem_meta.get_cell_topology(*part).getCellTopologyData() == topology) {
      if (skip_topology_root_parts && stk::mesh::fem::is_cell_topology_root_part(*part)) {
        continue;
      }
      parts.push_back(part);
    }
  }
}

inline
unsigned get_spatial_dimension(const Entity& entity)
{
  // expose some dot-chain to ensure everything inlined
  return entity.bucket().mesh().mesh_meta_data().get_spatial_dimension();
}

/* The Fmwk uses an enum to identify nodes, edges, faces and elements. The
   toolkit is similar, but the the toolkit rank depends on the spatial
   dimension. For 3D parts they are identical. With 2D, the element rank is
   2 (not 3). etc.
*/

#ifdef SIERRA_MIGRATION

inline
unsigned convert_fmwk_rank_to_stk(unsigned fmwk_rank, unsigned spatial_dim)
{
  ThrowAssert(spatial_dim > 0);
  ThrowAssert(spatial_dim < 4);
  ThrowAssert(fmwk_rank <= 4); // up to four basic entities types and constraints
  static int map2Stk[4][5] = { {-1,-1,-1,-1, 4},
                               { 0,-1,-1, 1, 4},
                               { 0, 1,-1, 2, 4},
                               { 0, 1, 2, 3, 4}};
  int stk_rank = map2Stk[spatial_dim][fmwk_rank];

  return static_cast<unsigned>(stk_rank);
}

inline
unsigned convert_stk_rank_to_fmwk(unsigned stk_rank,  unsigned spatial_dim)
{
  ThrowAssert(spatial_dim > 0);
  ThrowAssert(spatial_dim < 4);
  ThrowAssert(stk_rank <= 4); // up to four basic entities types and constraints
  static int map2Fmwk[4][5]={ {-1,-1,-1,-1, 4},
                              { 0, 3,-1,-1, 4},
                              { 0, 1, 3,-1, 4},
                              { 0, 1, 2, 3, 4}};
  int fmwk_rank=map2Fmwk[spatial_dim][stk_rank];

  return static_cast<unsigned>(fmwk_rank);
}

inline
unsigned get_derived_type(const Entity& entity)
{
  return convert_stk_rank_to_fmwk(entity.entity_rank(), get_spatial_dimension(entity));
}

#endif

/** \} */

} //namespace fem
} //namespace mesh
} //namespace stk
#endif
