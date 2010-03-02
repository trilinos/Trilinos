/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_EntityComm_hpp
#define stk_mesh_EntityComm_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <vector>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_mesh/base/Types.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

void pack_entity_info( CommBuffer & buf , const Entity & entity );
 
void unpack_entity_info(
  CommBuffer     & buf,
  const BulkData & mesh ,
  EntityKey      & key ,
  unsigned       & owner , 
  PartVector     & parts , 
  std::vector<Relation> & relations );
 

/** \brief  Pack an entity's field values into a buffer */
void pack_field_values( CommBuffer & , Entity & );

/** \brief  Unpack an entity's field values from a buffer */
bool unpack_field_values( CommBuffer & , Entity & , std::ostream & error_msg );

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

