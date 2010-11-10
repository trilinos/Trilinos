/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_TopologyHelpers_hpp
#define stk_mesh_TopologyHelpers_hpp

#include <sstream>
#include <stdexcept>
#include <Shards_CellTopologyTraits.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/FEMTypes.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/TopologicalMetaData.hpp>

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
#include <stk_mesh/fem/TopologyHelpersDeprecated.hpp>
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_bulk_data_element
 *  \{
 */

/// \todo REFACTOR: The functions in this file represent a "bridge"
///between the mesh and the Shards_CellTopologyData stuff. Does it
///belong here?

//----------------------------------------------------------------------
template< class Traits >
void get_parts_with_topology(stk::mesh::BulkData& mesh,
                             stk::mesh::PartVector& parts)
{
  parts.clear();

  const stk::mesh::PartVector& all_parts = mesh.mesh_meta_data().get_parts();

  stk::mesh::PartVector::const_iterator
    iter = all_parts.begin(),
    iter_end = all_parts.end();

  const CellTopologyData* topology = shards::getCellTopologyData<Traits>();

  for(; iter!=iter_end; ++iter) {
    stk::mesh::Part* part =  *iter;
#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
    if (get_cell_topology(*part) == topology) {
      parts.push_back(part);
    }
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
    if (TopologicalMetaData::get_cell_topology(*part) == topology) {
      parts.push_back(part);
    }
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  }
}

//----------------------------------------------------------------------
/** \brief  Declare an element member of a Part with a CellTopology
 *          and nodes conformal to that topology.
 */
template< typename IdType >
inline
Entity & declare_element( BulkData & mesh ,
                          Part & part ,
                          const IdType elem_id ,
                          const IdType node_id[] )
{
#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData * const top = get_cell_topology( part );
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData * const top = TopologicalMetaData::get_cell_topology( part );
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

  if ( top == NULL ) {
    std::ostringstream msg ;
    msg << "stk::mesh::declare_element( mesh , " ;
    msg << part.name();
    msg << " , " ;
    msg << elem_id ;
    msg << " , node_id[] ) ERROR, Part does not have a local topology" ;
    throw std::runtime_error( msg.str() );
  }

  PartVector empty ;
  PartVector add( 1 ); add[0] = & part ;

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const EntityRank entity_rank = element_rank_deprecated(part.mesh_meta_data());
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const EntityRank entity_rank = top->dimension;
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

  Entity & elem = mesh.declare_entity( entity_rank, elem_id, add );

  for ( unsigned i = 0 ; i < top->node_count ; ++i ) {
    //declare node if it doesn't already exist
    Entity * node = mesh.get_entity( NodeRank , node_id[i]);
    if ( NULL == node) {
      node = & mesh.declare_entity( NodeRank , node_id[i], empty );
    }

    mesh.declare_relation( elem , *node , i );
  }
  return elem ;
}

/** \brief  Declare an element member of a Part with a CellTopology
 *          and nodes conformal to that topology.
 */
template< typename IdType >
Entity & declare_element( BulkData & mesh ,
                          Part & part ,
                          const IdType elem_id ,
                          Entity * node[] )
{
#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData * const top = get_cell_topology( part );
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData * const top = TopologicalMetaData::get_cell_topology( part );
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

  if ( top == NULL ) {
    std::ostringstream msg ;
    msg << "stk::mesh::declare_element( mesh , " ;
    msg << part.name();
    msg << " , " ;
    msg << elem_id ;
    msg << " , node[] ) ERROR, Part does not have a local topology" ;
    throw std::runtime_error( msg.str() );
  }

  PartVector add( 1 ); add[0] = & part ;

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const EntityRank entity_rank = element_rank_deprecated(part.mesh_meta_data());
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const EntityRank entity_rank = top->dimension;
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

  Entity & elem = mesh.declare_entity( entity_rank, elem_id, add );

  for ( unsigned i = 0 ; i < top->node_count ; ++i ) {
    mesh.declare_relation( elem , *node[i] , i );
  }
  return elem ;
}

//----------------------------------------------------------------------
/** \brief  Create (or find) an element side.
 *
 *  The element must be a member of a Part with a CellTopology.
 */
Entity & declare_element_side( BulkData & mesh ,
                               const stk::mesh::EntityId global_side_id ,
                               Entity & elem ,
                               const unsigned local_side_id ,
                               Part * part = NULL );

Entity & declare_element_side( Entity & elem ,
                               Entity & side ,
                               const unsigned local_side_id ,
                               Part * part = NULL );

/** \brief  Determine the polarity of the local side,
 *          more efficient if the local_side_id is known.
 */
bool element_side_polarity( const Entity & elem ,
                            const Entity & side , int local_side_id = -1 );


/** \brief  Given an element and collection of nodes, return the
 *          local id of the side that contains those nodes in the
 *          correct orientation.
 */
int element_local_side_id( const Entity & elem ,
                           const CellTopologyData * side_topology,
                           const EntityVector & side_nodes );

//----------------------------------------------------------------------

/** \brief Given an element and a side ordinal, populate a vector of nodes that make up the side.
 * The nodes are ordered such that the correct node ordering for the side is preserved and the node
 * with the lowest identifier comes first
 *
 * return the CellTopologyData * for the given side if it exist, else return null
 */
const CellTopologyData * get_elem_side_nodes( const Entity & elem,
                          RelationIdentifier side_ordinal,
                          EntityVector & side_key_nodes
                        );
/** \} */

}//namespace mesh
}//namespace stk

/** 11/8/10 Design meeting to decide on FEM layer interfaces
 * Attendance:  Carter, Alan, Greg, Jim, Dave, Todd, Dan
 * Issues to resolve:
 *   * EntityRanks defined by discretization names vs topological names
 *   * discretization names:  vertex, side, element
 *   * topological names:  node, edge, face, solid (3D)
 *   Element Rank depends on (equals) spatial dimension  [ essential complexity ]
 *   Entity Rank is not 1:1 with respect to topological dimension
 *
 * Three fundamental components of FEM: (agreed on by all in attendance)
 *   1.  Cell Topology
 *   2.  Entity Ranks
 *   3.  Spatial Dimension
 *
 * Where do coordinates go?  In Fields...
 * Where does axi-symmetric information go?
 *
 * Issue of attribute mechanism on MetaData giving const only access:
 * Carter said this is an arbitrary restriction and it should be removed now that we have a need.
 * All in agreement that attributes should be const/nonconst accessible.
 *
 * Issue of Pure Abstract Class design for FEM layer.
 *   * Initial idea was to allow applications store mapping of parts to cell
 *     topologies in different ways.  As a map on MetaData or as attributes on
 *     parts.  Alan suggested we store them in an std::vector<cell_topology*>
 *     with one index into vector for every Part (ordinal).  This would make
 *     lookup constant which is about as fast as you can make it. 
 *   * What about a template for this variation?
 *     Concern that this would move implementation into headers and require
 *     templatization of stk_io.  
 *   * How does refinement topology get added into the mix?  Is this an
 *     extension of FEM or a different plugin entirely?  Carter said a
 *     complaint in the framework has been "Why do I have to see refinement
 *     topology when I don't need it?"  This is an argument for a different
 *     plug-in.
 *   * Pro/Con for abstract base class:
 *     Pro:  specialization allowed easily (axi-symmetric information,
 *       refinement topology, alternative implementation of cell_topology
 *       lookups, etc.), user-defined ranks, Java style "interface" API
 *     Con:  unessential complexity added, "abstract", no other examples in
 *       stk_mesh, no need for it now, virtual symbol table, no inline calls,
 *       non-obvious performance concerns (e.g. fem.element_rank() looks
 *       light-weight, but may not be).
 *   * Favor composition over inheritance.  Where does this come into play here?
 *   * What is lesson learned from Framework use of inheritance?  Dave & Carter
 *     discussed this and it appeared that in the framework the base classes
 *     were not pure-virtual "interfaces" and the resulting extensive
 *     inheritance caused great problems.  The idea here is that a pure-virtual
 *     base class is safer.  There is also the idea that this class should only
 *     be inherited from once.
 *   * Voting on abstract vs. concrete implementation:
 *     Abstract:  Todd, Alan, Greg, Dan.
 *     Concrete:  Carter, Jim.
 *     Abstain:  Dave.
 *     Greg comment:  concrete classes that need specializations lead to work-arounds.
 *     Dan comment:  inner loop performance issues are easy to resolve.
 *   * We also talked about how to access the element_rank EntityRank
 *     information.  Currently Dave has some free functions that use MetaData
 *     to get the attribute for the FEM interface and it gets spatial_dimension
 *     off of that class and computes the element rank.  Alan suggested that we
 *     move the element_rank function to the FEM interface class instead and
 *     then we will use a block of code like:
 *       FEMInterface & fem = get_fem(part);
 *       fem.element_rank();
 *     In this way, the DefaultFEM class would set the internal element_rank
 *     value at construction with the spatial dimension and then the accessor
 *     would simply return the value.
 *
 *  Issue of constructor for MetaData:
 *    The main issue is that MetaData requires a list of entity rank names at
 *      construction and it uses this internally to check invariants related to
 *      the maximum entity rank and it is used in BulkData for both the maximum
 *      entity rank and the actual entity rank names in error messages.  When
 *      reading from an exodus file, you need to pass MetaData into stk_io but
 *      you don't know if the problem is 2D or 3D yet and this is necessary to
 *      create the correct list of entity rank names.
 *    Dave added a default MetaData constructor and an additional state for
 *      before/after setting entity rank names.  This allowed him to create a
 *      MetaData and pass it to stk_io to fill up with MetaData and entity rank
 *      names based on spatial dimension from the exodus file.
 *    Dan/Greg agree that stk_io needs better code for interacting with exodus,
 *      MetaData, and BulkData.  Currenlty, there is a use_case_mesh class that
 *      has a nice set of features but because its in a use case its not really
 *      part of the library, yet it is being used heavily clients of stk_mesh to
 *      read in exodus files.
 *    Alan suggested that we should have a factory in stk_io which creates
 *      MetaData and passes it back through a shared pointer.  
 *    We had a discussion about whether entity ranks must be contiguous or not.
 *      It appears that base is okay with non-contiguous entty ranks, but FEM
 *      requires them to be contiguous.
 *    We then had a discussion about who is a client of base and not FEM.  It
 *      appears that arbitrary polyhedra is building on top of FEM and
 *      peridynamics is also building on top of FEM.  Those were the two big
 *      examples that were supposed to build on top of base directly and not
 *      use the FEM layer.  So at this point we don't have even a single
 *      use-case for base all by itself.  This raised the question of whether
 *      we should just move the FEM interface calls into MetaData directly and
 *      access an internal FEM class to store all the mappings.  This part of
 *      the conversation occured at the end of 2.5 hours of discussion and
 *      trailed off with no conclusion.
 *    Dave is going to implement an additional state in MetaData for the
 *      before/after setting entity rank names and implement checks for this
 *      state in the appropriate member functions in MetaData.  This will allow
 *      his current implementation with stk_io to move forward easily.
 *    
 **/
#endif

