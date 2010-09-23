/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_TopologicalMetaData_hpp
#define stk_mesh_TopologicalMetaData_hpp

#include <Shards_CellTopologyData.h>
#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/fem/FEMTypes.hpp>

namespace stk {
namespace mesh {


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  Augment mesh meta data with topological meta data including
 *          entity ranks, shards' cell topologies appropriate for
 *          the given spatial dimension, and part to cell topology mapping.
 */
class TopologicalMetaData {
public:
  //--------------------------------------------------------------------------

  static void verify_spatial_dimension( unsigned spatial_dimension ,
                                        const char * method );

  /** \brief  Entity rank text names that match the entity rank values in this class */
  static std::vector<std::string> entity_rank_names( unsigned spatial_dimension );

  /** \brief  Construct topological meta data including the map
   *          from Parts to shards' cell topologies.
   *          Required that 1 <= spatial_dimension <= 3 .
   */
  TopologicalMetaData( MetaData & , unsigned spatial_dimension );

  ~TopologicalMetaData();

  //--------------------------------------------------------------------------

  const unsigned spatial_dimension ; ///< Spatial dimension

  const EntityRank node_rank ;    ///< Rank of nodes (always zero)
  const EntityRank edge_rank ;    ///< Rank of edges (1 for 2D and 3D)
  const EntityRank side_rank ;    ///< Rank of sides (1 for 2D, 2 for 3D)
  const EntityRank element_rank ; ///< Rank of elements (spatial_dimension)
  const EntityRank patch_rank ;   ///< Rank of arbitrary patches (element+1)

  //--------------------------------------------------------------------------
  /** \brief  Get the cell topology associated with a mesh part.
   *          The supersets of the mesh part are not checked.
   */
  static const CellTopologyData * get_cell_topology( const Part & part , const char * required_by = 0 );

/** \brief  Query the cell topology associated with a bucket.
 *
 *          The bucket's superset mesh parts of the bucket's entity rank
 *          are queried for their associated cell topology.
 *
 *          If no cell topology is found and required_by is zero then
 *          zero is returned, otherwise an exception is thrown with the
 *          required_by string in the message.
 *
 *          If there is exactly one cell topology then it is returned.
 *
 *          If more than one cell topology is found then an exception
 *          is thrown.
 */
  static const CellTopologyData * get_cell_topology( const Bucket & bucket , const char * required_by = 0 );
  static const CellTopologyData * get_cell_topology( const Entity & entity , const char * required_by = 0 );

  //--------------------------------------------------------------------------
  /** \brief  Declare part of a given name and associated cell topology.
   *
   *          If the cell topology is not already defined
   *          then the entity rank is defined to be the
   *          dimension of the cell topology.
   */
  Part & declare_part( const std::string & name ,
                       const CellTopologyData * top );

  /** \brief  Declare part of a given name and associated cell topology */
  template< class Top >
  Part & declare_part( const std::string & name )
    { return declare_part( name , shards::getCellTopologyData<Top>() ); }

  //--------------------------------------------------------------------------
  /** \brief  Extend the list of defined cell topologies for
   *          mesh entities of the given rank.
   */
  void declare_cell_topology( const CellTopologyData * , EntityRank entity_rank );

  //throws if not found
  EntityRank get_entity_rank( const CellTopologyData * ) const ;

private:
 
  // Base meta data to which this topological meta data is attached.
  MetaData & m_meta_data ;

  // Defined cell topologies and associated mesh entity rank.
  std::vector< std::pair< const CellTopologyData * , EntityRank > > m_top_rank ;

  // Map part meta data ordinals to cell topologies.
  std::vector< std::pair< unsigned , const CellTopologyData * > > m_part_top_map ;

  static const TopologicalMetaData * internal_get( const MetaData & );

  const CellTopologyData * internal_get_cell_topology( unsigned part_ordinal) const ;

  void internal_set_entity_rank( const CellTopologyData * , EntityRank );
  
  // Set sell topology internally, to be used only from declare_part
  void internal_set_cell_topology( 
      Part & part, 
      unsigned entity_rank, 
      const CellTopologyData * top
      );


  void throw_ambiguous( const Part & ) const ;
  void throw_ambiguous( const Bucket & ) const ;
  static void throw_required( const Part & , const char * );
  static void throw_required( const Bucket & , const char * );


  TopologicalMetaData();
  TopologicalMetaData( const TopologicalMetaData & );
  TopologicalMetaData & operator = ( const TopologicalMetaData & );
};

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

