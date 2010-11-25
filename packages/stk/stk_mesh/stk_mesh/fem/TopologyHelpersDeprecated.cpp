/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards
 */

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <cassert>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologyHelpersDeprecated.hpp>

#include <stk_mesh/fem/DefaultFEM.hpp>

#include <stk_util/util/StaticAssert.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

// DEPRECATED: 09/15/10 FEM refactor
const CellTopologyData * get_cell_topology_deprecated( const Part & p )
{
  return p.attribute<CellTopologyData>();
}
// DEPRECATED: 09/15/10 FEM refactor
const CellTopologyData * get_cell_topology( const Part & p )
{
  const fem::FEMInterface * fem = p.mesh_meta_data().get_attribute<fem::FEMInterface>();
  const CellTopologyData * cell_topology_data;
  if (fem) {
    cell_topology_data = fem::get_cell_topology(p).getTopologyData();
  } else {
    cell_topology_data = get_cell_topology_deprecated(p);
  }
  return cell_topology_data;
}


// DEPRECATED: 09/15/10 FEM refactor
void set_cell_topology( Part & p , const CellTopologyData * singleton )
{
  set_cell_topology_deprecated( p, singleton );
}


// DEPRECATED: 09/15/10 FEM refactor
void set_cell_topology_deprecated( Part & p , const CellTopologyData * singleton )
{
  static const char method[] = "stk::mesh::set_cell_topology" ;

  MetaData & m = p.mesh_meta_data();

  const CellTopologyData * t = NULL ;

  if ( p.mesh_meta_data().entity_rank_count() <= p.primary_entity_rank() ||
       singleton == NULL ||
       singleton != ( t = m.declare_attribute_no_delete(p,singleton) ) ) {
    std::ostringstream msg ;
    msg << method << "( " << p.name();
    msg << " entity_rank(" << p.primary_entity_rank() << ") , " ;
    if ( singleton ) { msg << singleton->name ; }
    else             { msg << "NULL" ; }
    msg << " ) ERROR" ;
    if ( t ) { msg << "Existing topology = " << t->name ; }
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------


// DEPRECATED: 09/15/10 FEM refactor
const CellTopologyData * get_cell_topology_deprecated( const Bucket & bucket )
{
  const CellTopologyData * top = NULL ;
  PartVector parts ;
  bucket.supersets( parts );

  PartVector::iterator i = parts.begin() ;

  for ( ; NULL == top && i != parts.end() ; ++i ) {
    if ( bucket.entity_rank() == (**i).primary_entity_rank() ) {
      top = get_cell_topology( **i );
    }
  }

  bool ok = true ;

  for ( ; ok && i != parts.end() ; ++i ) {
    if ( bucket.entity_rank() == (**i).primary_entity_rank() ) {
      const CellTopologyData * const tmp = get_cell_topology( **i );
      ok = ((tmp == NULL) || (tmp == top)) ;
    }
  }

  if ( ! ok ) {
    std::ostringstream msg ;
    msg << "stk::mesh::get_cell_topology( Bucket[" ;
    for ( i = parts.begin() ; i != parts.end() ; ++i ) {
      if ( bucket.entity_rank() == (**i).primary_entity_rank() ) {
        const CellTopologyData * const tmp = get_cell_topology( **i );
        msg << " " << (*i)->name();
        if ( tmp ) { msg << "->" << tmp->name ; }
        msg << " ] ).getTopologyData() FAILED WITH MULTIPLE LOCAL TOPOLOGIES" ;
        throw std::runtime_error( msg.str() );
      }
    }
  }
  return top ;
}


// DEPRECATED: 09/15/10 FEM refactor
const CellTopologyData * get_cell_topology( const Bucket & bucket )
{
  const fem::FEMInterface * fem = bucket.mesh().mesh_meta_data().get_attribute< fem::FEMInterface >();
  const CellTopologyData * cell_topology_data;
  if (fem) {
    cell_topology_data = fem::get_cell_topology(bucket).getTopologyData();
  } else {
    cell_topology_data = get_cell_topology_deprecated(bucket);
  }
  return cell_topology_data;
}


// DEPRECATED: 09/15/10 FEM refactor
const CellTopologyData * get_cell_topology_deprecated( const Entity & entity )
{ return get_cell_topology_deprecated( entity.bucket() ); }


// DEPRECATED: 09/15/10 FEM refactor
const CellTopologyData * get_cell_topology( const Entity & entity )
{ return get_cell_topology(entity.bucket()); }

//----------------------------------------------------------------------

}// namespace mesh
}// namespace stk

#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
