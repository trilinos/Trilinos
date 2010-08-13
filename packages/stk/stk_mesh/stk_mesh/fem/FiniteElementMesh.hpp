/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_FEMesh_hpp
#define stk_mesh_FEMesh_hpp

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_mesh/fem/FieldTraits.hpp>

#include <stk_mesh/femImpl/FiniteElementMeshImpl.hpp>
#include <stk_mesh/femImpl/PartCellTopologyMap.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  A Finite Element Mesh specializes a "generic" STK-Mesh with
 *          finite element topologies from the Shards library.
 *
 *  The spatial dimension of the mesh determines the entity-ranks
 *  and default cell topologies extracted from the Shards library.
 */

template< class ArgCoordinateFieldType = Field<double,Cartesian> >
class FiniteElementMesh {
public:
  //--------------------------------------------------------------------------

  typedef ArgCoordinateFieldType CoordinateFieldType ;

  MetaData metaData ;
  BulkData bulkData ;

  CoordinateFieldType & coordinate_field ;

  const unsigned spatial_dimension ; ///< Spatial dimension

  const unsigned node_rank ;    ///< Rank of nodes (always zero)
  const unsigned edge_rank ;    ///< Rank of edges (1 for 2D and 3D)
  const unsigned side_rank ;    ///< Rank of sides (1 for 2D, 2 for 3D)
  const unsigned element_rank ; ///< Rank of elements (spatial_dimension)
  const unsigned patch_rank ;   ///< Rank of arbitrary patches (element+1)

  //--------------------------------------------------------------------------
  /** \brief  Declare part of a given name and cell topology */

  Part & declare_part( const std::string & name ,
                       const CellTopologyData & top )
  {
    static const char method[] = "stk::mesh::FiniteElementMesh::declare_part" ;

    Part & top_part = * m_part_cell_topology_map.get_part( top , method );

    Part & part = metaData.declare_part( name, top_part.primary_entity_rank() );

    metaData.declare_part_subset( top_part , part );

    return part ;
  }

  template< class Top >
  Part & declare_part( const std::string & name )
  { return declare_part( name , * shards::getCellTopologyData<Top>() ); }

  //--------------------------------------------------------------------------
  /** \brief  Query the unique part declared for the given cell topology.
   *          Returns NULL if there is no associated part. 
   */
  Part * query_part( const CellTopologyData & top ) const
    { return m_part_cell_topology_map.get_part( top , NULL ); }

  //--------------------------------------------------------------------------
  /** \brief  Construct a finite element mesh and populate its
   *          meta data with the coordinate field on every node
   *          and shards' predefined cell topologies appropriate for
   *          the coordinate system and dimension.
   *
   *          The mesh meta data is not committed.
   */
  FiniteElementMesh(
    unsigned spatial_dimension , ParallelMachine parallel_machine )
  : metaData( impl::finite_element_mesh_entity_rank_names(spatial_dimension) )
  , bulkData( metaData , parallel_machine )
  , coordinate_field(
      metaData.declare_field<CoordinateFieldType>("coordinates",1) )
  , spatial_dimension( spatial_dimension )
  , node_rank( 0 )
  , edge_rank( 1 < spatial_dimension ? 1 : 0 )
  , side_rank( 2 < spatial_dimension ? 2 : edge_rank )
  , element_rank( spatial_dimension )
  , patch_rank( spatial_dimension + 1 )
  , m_part_cell_topology_map( metaData , spatial_dimension )
  {
    put_field(
       coordinate_field, 1u, metaData.universal_part(), spatial_dimension );
  }

  ~FiniteElementMesh() {}

  //--------------------------------------------------------------------------
  /** \brief  Query or declare the unique part with the given cell topology.
   *
   *  Parts for the conventional cell topologies are declared by default.
   *  This method provides an "extension point" in the design for
   *  applications to introduce specialized cell topologies as needed.
   */
  Part & declare_part( const CellTopologyData & top )
    { return m_part_cell_topology_map.declare_part( top , top.dimension ); }

private: 

  impl::PartCellTopologyMap m_part_cell_topology_map ;

  FiniteElementMesh();
  FiniteElementMesh( const FiniteElementMesh & );
  FiniteElementMesh & operator = ( const FiniteElementMesh & );
};

//----------------------------------------------------------------------------
/** \brief  Query topology associated with a part.
 *
 *          If the input part does not have an associated topology
 *          then the supersets of the part are queried.
 *          If there is exactly one cell topology then it is returned,
 *          otherwise NULL is returned.
 */
inline
const CellTopologyData * get_cell_topology(
  const Part & part , const char * const required_by = NULL )
{ return impl::PartCellTopologyMap::get_cell_topology( part , required_by ); }

/** \brief  Query topology associated with a set of homogeneous entities. */
inline
const CellTopologyData * get_cell_topology(
  const Bucket & bucket , const char * const required_by = NULL )
{ return impl::PartCellTopologyMap::get_cell_topology( bucket , required_by ); }

/** \brief  Query topology associated with an entity. */
inline
const CellTopologyData * get_cell_topology(
  const Entity & entity , const char * const required_by = NULL )
{ return impl::PartCellTopologyMap::get_cell_topology( entity.bucket() , required_by ); }

//----------------------------------------------------------------------------
/** \brief  Given an element and collection of side nodes
 *          determine the subcell and permutation indices for the
 *          local side topology that matches the given nodes.
 *          ( local_side_index , local_side_permutation )
 *
 *          If their is no matching side then return (-1,-1).
 */
std::pair<int,int>
  element_side_local_info( const Entity & element ,
                           const std::vector<Entity*> & side_nodes );

//----------------------------------------------------------------------------
/** \brief  Declare an element with a defined cell topology.
 *
 *  The 'element_part' must have an associated topology or be the subset
 *  of a part that has an associated topology.  The length of the 
 *  'node_ids' array must be the number of nodes of that topology.
 *  The returned entity will be of rank 'element_rank'.
 *  Nodes of the given identifiers will be declared or queried.
 */
template< typename IdType >
Entity & declare_element( BulkData     & bulk_data ,
                          Part         & element_part ,
                          const IdType & element_id ,
                          const IdType   node_ids[] )
{
  static const char method[] = "stk::mesh::declare_element" ;

  const unsigned node_rank = 0 ; // This is hard-wired for any FEMesh
  const unsigned elem_rank = element_part.primary_entity_rank();

  const CellTopologyData * const top =
    get_cell_topology( element_part , method );

  PartVector empty ;
  PartVector add( 1 ); add[0] = & element_part ;

  Entity & element = bulk_data.declare_entity( elem_rank , element_id, add );

  for ( unsigned i = 0 ; i < top->node_count ; ++i ) {
    Entity & node = bulk_data.declare_entity( node_rank , node_ids[i], empty );
    bulk_data.declare_relation( element , node , i );
  }
  return element ;
}

//--------------------------------------------------------------------------
/** \brief  Declare the side of an element.
 *
 *  The element must have an associated topology.
 *  No attempt is made to determine if coincident side
 *  already exists attached to a coincident or adjacent element.
 */
template< class FEMesh , typename IdType >
Entity & declare_element_side( FEMesh & mesh ,
                               Part   & side_part ,
                               Entity & element ,
                               const IdType & side_global_id ,
                               const IdType & side_local_ordinal );

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

